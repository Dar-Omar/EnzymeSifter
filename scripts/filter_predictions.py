#!/usr/bin/env python3

import argparse
import csv
import sys
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
EITHER        = {"tm"}

FILTERABLE_STEMS = {
    "all_predictions",
    "all_predictions_structure",
    "all_predictions_chains",
}


class Spec:
    """Parsed filter spec for one flag."""
    __slots__ = ("flag", "column", "mode", "lo", "hi")

    def __init__(self, flag, column, mode, lo, hi):
        self.flag   = flag
        self.column = column
        self.mode   = mode
        self.lo     = lo
        self.hi     = hi

    def passes(self, value_str):
        if value_str in ("NA", "", None):
            return True
        try:
            v = float(value_str)
        except ValueError:
            return True
        if self.mode == "cutoff":
            return v >= self.lo
        return self.lo <= v <= self.hi


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


def filter_file(input_path, output_path, specs):
    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        headers = list(reader.fieldnames or [])
        rows = list(reader)

    active = [s for s in specs if s.column in headers]
    if not active:
        print(f"[WARN] No filterable columns in {input_path.name} — "
              f"copying as-is.", file=sys.stderr)
        filtered = rows
    else:
        filtered = [r for r in rows if all(s.passes(r.get(s.column, "NA"))
                                           for s in active)]

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(filtered)

    print(f"[INFO] {input_path.name}: {len(filtered)}/{len(rows)} rows "
          f"passed → {output_path.name}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Filter prediction TSVs by threshold/interval flags."
    )
    parser.add_argument("--output_dir", required=True,
                        help="predictions_output directory")
    for flag in COLUMN_FLAG_MAP:
        parser.add_argument(f"--{flag}", type=str, default=None,
                            help=f"Spec for {COLUMN_FLAG_MAP[flag]} "
                                 f"(number or min:max).")
    args = parser.parse_args()

    specs = []
    for flag in COLUMN_FLAG_MAP:
        raw = getattr(args, flag)
        s = parse_spec(flag, raw)
        if s is not None:
            specs.append(s)

    if not specs:
        print("[INFO] No thresholds provided — nothing to filter.",
              file=sys.stderr)
        return

    summary = []
    for s in specs:
        if s.mode == "cutoff":
            summary.append(f"{s.column} >= {s.lo}")
        else:
            summary.append(f"{s.column} in [{s.lo}, {s.hi}]")
    print(f"[INFO] Active filters: {', '.join(summary)}", file=sys.stderr)

    out_dir = Path(args.output_dir)
    for tsv in sorted(out_dir.glob("*.tsv")):
        if tsv.stem.endswith("_filtered"):
            continue
        if tsv.stem not in FILTERABLE_STEMS:
            continue
        filtered_path = tsv.with_name(tsv.stem + "_filtered.tsv")
        filter_file(tsv, filtered_path, specs)


if __name__ == "__main__":
    main()
