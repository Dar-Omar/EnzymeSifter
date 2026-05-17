#!/usr/bin/env python3
"""
Inputs:
  --clades     TSV: tip_name, clade_id
  --hit_list   text file with one pdb basename per line (used to
               disambiguate single-chain PDBs whose filename ends in
               _<digits>, e.g. 'zt2_1')
  --filtered   TSV(s) of threshold-passing enzymes (no score column)
  --output     Output TSV path
  Filter specs (same syntax as filter_predictions.py; used for scoring):
    --solubility, --usability, --phopt, --topt, --tm
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

COLUMN_FLAG_MAP = {
    "solubility": "predicted_solubility",
    "usability":  "predicted_usability",
    "phopt":      "predicted_ph_opt",
    "topt":       "predicted_topt_C",
    "tm":         "predicted_tm_C",
}
CUTOFF_ONLY   = {"solubility", "usability"}
INTERVAL_ONLY = {"phopt", "topt"}


class Spec:
    __slots__ = ("flag", "column", "mode", "lo", "hi")

    def __init__(self, flag, column, mode, lo, hi):
        self.flag   = flag
        self.column = column
        self.mode   = mode
        self.lo     = lo
        self.hi     = hi

    def midpoint(self):
        return (self.lo + self.hi) / 2.0

    def half_width(self):
        return (self.hi - self.lo) / 2.0

    def score(self, value_str, pool_min, pool_max):
        """Per-property score in [0, 1]; None if value is NA/unparseable."""
        if value_str in ("NA", "", None):
            return None
        try:
            v = float(value_str)
        except ValueError:
            return None
        if self.mode == "cutoff":
            lo = max(self.lo, pool_min)
            hi = pool_max
            if hi <= lo:
                return 1.0
            return max(0.0, min(1.0, (v - lo) / (hi - lo)))                #imp
        hw = self.half_width()
        if hw <= 0:
            return 1.0
        return max(0.0, min(1.0, 1.0 - abs(v - self.midpoint()) / hw))     #imp


def parse_spec(flag, raw):
    if raw is None:
        return None
    column = COLUMN_FLAG_MAP[flag]
    has_colon = ":" in raw
    if has_colon:
        if flag in CUTOFF_ONLY:
            sys.exit(f"[ERROR] --{flag} does not accept an interval "
                     f"(got {raw!r}); provide a single number.")
        try:
            lo_str, hi_str = raw.split(":", 1)
            lo, hi = float(lo_str), float(hi_str)
        except ValueError:
            sys.exit(f"[ERROR] --{flag} interval must be 'min:max' "
                     f"(got {raw!r}).")
        if hi <= lo:
            sys.exit(f"[ERROR] --{flag} interval requires max > min "
                     f"(got {lo}:{hi}).")
        return Spec(flag, column, "interval", lo, hi)
    else:
        if flag in INTERVAL_ONLY:
            sys.exit(f"[ERROR] --{flag} requires an interval 'min:max' "
                     f"(got {raw!r}).")
        try:
            lo = float(raw)
        except ValueError:
            sys.exit(f"[ERROR] --{flag} must be a number (got {raw!r}).")
        return Spec(flag, column, "cutoff", lo, None)


def tip_to_pdb(tip_name, known_pdbs=None):
    # Handles the presence of _ in IDs
    if known_pdbs is not None and tip_name in known_pdbs:
        return tip_name
    parts = tip_name.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return tip_name


def load_clades(path):
    tip_to_clade = {}
    members = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tip = row.get("tip_name", "").strip()
            cid = row.get("clade_id", "").strip()
            if not tip or not cid:
                continue
            cid_int = int(cid)
            tip_to_clade[tip] = cid_int
            members[cid_int].append(tip)
    return tip_to_clade, members, sorted(members.keys())


def load_filtered(paths):
    rows_by_pdb = defaultdict(list)
    extra_cols = []

    RESERVED = {"ID", "score", "chain_id", "clade_id",
                "n_eligible_in_clade", "n_members_in_clade"}

    for p in paths:
        if not Path(p).exists():
            print(f"[WARN] Filtered TSV not found: {p}", file=sys.stderr)
            continue
        with open(p, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            headers = reader.fieldnames or []
            for h in headers:
                if h in RESERVED:
                    continue
                if h not in extra_cols:
                    extra_cols.append(h)
            for row in reader:
                pid = row.get("ID", "").strip()
                if not pid:
                    continue
                rows_by_pdb[pid].append(row)

    eligible = {}
    excluded = []
    for pid, entries in rows_by_pdb.items():
        if len(entries) == 1:
            eligible[pid] = entries[0]
        else:
            excluded.append(pid)

    if excluded:
        shown = ", ".join(sorted(excluded)[:10])
        more = "" if len(excluded) <= 10 else f" (+{len(excluded) - 10} more)"
        print(f"[INFO] Excluded {len(excluded)} multi-sequence PDB(s) from "
              f"representative selection: {shown}{more}", file=sys.stderr)

    return eligible, extra_cols, len(excluded)



def compute_pool_extremes(rows, specs):
    """For each cutoff-mode spec, observed (min, max) across the pool."""
    extremes = {}
    for s in specs:
        if s.mode != "cutoff":
            continue
        vals = []
        for row in rows:
            raw = row.get(s.column, "NA")
            if raw in ("NA", "", None):
                continue
            try:
                vals.append(float(raw))
            except ValueError:
                continue
        extremes[s.flag] = (min(vals), max(vals)) if vals else (s.lo, s.lo)
    return extremes


def row_score(row, specs, pool_extremes):
    """Combined score for one row; None if nothing scorable."""
    parts = []
    for s in specs:
        if s.column not in row:
            continue
        if s.mode == "cutoff":
            pmin, pmax = pool_extremes.get(s.flag, (s.lo, s.lo))
            sc = s.score(row.get(s.column, "NA"), pmin, pmax)
        else:
            sc = s.score(row.get(s.column, "NA"), None, None)
        if sc is not None:
            parts.append(sc)
    if not parts:
        return None
    return sum(parts) / len(parts)



def main():
    p = argparse.ArgumentParser(
        description="Select best-scoring representative per clade "
                    "(single-sequence enzymes only)."
    )
    p.add_argument("--clades", required=True,
                   help="clade_assignments.tsv")
    p.add_argument("--hit_list", required=False, default=None,
                   help="data/enzymm/hit_pdbs.txt — used to disambiguate "
                        "single-chain PDBs whose filename ends in _<digits>")
    p.add_argument("--filtered", required=True, nargs="+",
                   help="One or more *_filtered.tsv files")
    p.add_argument("--output", required=True,
                   help="Output TSV (clade_representatives.tsv)")
    for flag in COLUMN_FLAG_MAP:
        p.add_argument(f"--{flag}", type=str, default=None,
                       help=f"Scoring spec for {COLUMN_FLAG_MAP[flag]} "
                            f"(number or min:max).")
    args = p.parse_args()

    specs = []
    for flag in COLUMN_FLAG_MAP:
        raw = getattr(args, flag)
        s = parse_spec(flag, raw)
        if s is not None:
            specs.append(s)

    if not specs:
        sys.exit("[ERROR] At least one scoring spec is required "
                 "(e.g. --tm 50:80). Without specs there's no way to "
                 "rank enzymes.")

    known_pdbs = None
    if args.hit_list:
        with open(args.hit_list) as fh:
            known_pdbs = {ln.strip() for ln in fh if ln.strip()}
        print(f"[INFO] Loaded {len(known_pdbs)} known pdb_id(s) from "
              f"{args.hit_list}", file=sys.stderr)

    tip_to_clade, members_by_clade, all_clades = load_clades(args.clades)
    eligible, extra_cols, n_excluded = load_filtered(args.filtered)

    print(f"[INFO] {len(all_clades)} clade(s); "
          f"{len(eligible)} eligible single-sequence enzyme(s)",
          file=sys.stderr)

    sample_cols = set()
    for row in eligible.values():
        sample_cols.update(row.keys())
    active_specs = [s for s in specs if s.column in sample_cols]
    if not active_specs:
        sys.exit("[ERROR] None of the scoring specs match columns present "
                 "in the filtered TSVs.")

    pool_extremes = compute_pool_extremes(eligible.values(), active_specs)

    eligible_by_clade = defaultdict(list)
    for pid, row in eligible.items():
        assigned = None
        for tip, cid in tip_to_clade.items():
            if tip_to_pdb(tip, known_pdbs) == pid:
                assigned = cid
                break
        if assigned is None:
            continue
        sc = row_score(row, active_specs, pool_extremes)
        sc_num = sc if sc is not None else float("-inf")
        eligible_by_clade[assigned].append((sc_num, sc, row))

    n_members = {cid: len(members_by_clade[cid]) for cid in all_clades}

    out_headers = (["clade_id", "ID", "score",
                    "n_eligible_in_clade", "n_members_in_clade"]
                   + extra_cols)

    rows_out = []
    for cid in all_clades:
        entries = eligible_by_clade.get(cid, [])
        n_elig = len(entries)
        if n_elig == 0:
            rows_out.append({
                "clade_id": cid, "ID": "NA", "score": "NA",
                "n_eligible_in_clade": 0,
                "n_members_in_clade": n_members[cid],
            })
            continue
        entries.sort(key=lambda e: (-e[0], e[2].get("ID", "")))
        _, sc_real, winner_row = entries[0]
        row = {
            "clade_id": cid,
            "ID": winner_row.get("ID", ""),
            "score":  f"{sc_real:.4f}" if sc_real is not None else "NA",
            "n_eligible_in_clade": n_elig,
            "n_members_in_clade":  n_members[cid],
        }
        for col in extra_cols:
            row[col] = winner_row.get(col, "NA")
        rows_out.append(row)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_headers, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    n_winners = sum(1 for r in rows_out if r["ID"] != "NA")
    print(f"[INFO] Wrote {len(rows_out)} clade row(s); "
          f"{n_winners} winner(s), "
          f"{len(rows_out) - n_winners} empty clade(s) → {out}",
          file=sys.stderr)
    if n_excluded > 0:
        print(f"[INFO] {n_excluded} multi-sequence PDB(s) were not "
              f"considered.", file=sys.stderr)


if __name__ == "__main__":
    main()
