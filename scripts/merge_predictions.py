#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def base_id(seq_id: str) -> str:
    return seq_id.split("|")[0]


def seq_id_to_pdb_id(seq_id: str, known_pdbs=None) -> str:
    b = base_id(seq_id)
    if known_pdbs is not None and b in known_pdbs:
        return b
    parts = b.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return b


def load_tsv(path: str, label: str):
    p = Path(path)
    if not p.exists():
        print(f"[WARN] {label} file not found: {path}", file=sys.stderr)
        return [], []
    with open(p, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    headers = reader.fieldnames or []
    print(f"[INFO] Loaded {label}: {len(rows)} rows, cols={headers}", file=sys.stderr)
    return headers, rows


def main():
    parser = argparse.ArgumentParser(
        description="Merge all predictions into one or two TSVs."
    )
    parser.add_argument("--hit_list",   required=True)
    parser.add_argument("--netsolp",    required=True)
    parser.add_argument("--phoptnn",    required=True)
    parser.add_argument("--seq2topt",   required=True)
    parser.add_argument("--seq2tm",     required=True)
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write output file(s) into")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.hit_list) as fh:
        hit_pdbs = [line.strip() for line in fh if line.strip()]

    if not hit_pdbs:
        print("[WARN] No hits — writing empty output.", file=sys.stderr)
        with open(out_dir / "all_predictions.tsv", "w") as f:
            f.write("ID\n")
        return

    print(f"[INFO] {len(hit_pdbs)} hit PDB(s)", file=sys.stderr)
    known_pdbs = set(hit_pdbs)

    _, netsolp_rows = load_tsv(args.netsolp,  "NetSolP")
    _, phoptnn_rows = load_tsv(args.phoptnn,  "pHoptNN")
    _, topt_rows    = load_tsv(args.seq2topt, "Seq2Topt")
    _, tm_rows      = load_tsv(args.seq2tm,   "Seq2Tm")


    ns_by_chain     = {}
    topt_by_chain   = {}
    tm_by_chain     = {}
    chains_by_pdb   = defaultdict(list)

    def _register(cid, pid):
        if cid not in chains_by_pdb[pid]:
            chains_by_pdb[pid].append(cid)

    for row in netsolp_rows:
        sid = row.get("sid", "")
        cid, pid = base_id(sid), seq_id_to_pdb_id(sid, known_pdbs)
        ns_by_chain[cid] = {
            "predicted_solubility": row.get("predicted_solubility", ""),
            "predicted_usability":  row.get("predicted_usability",  ""),
        }
        _register(cid, pid)

    for row in topt_rows:
        sid = row.get("seq_id", "")
        cid, pid = base_id(sid), seq_id_to_pdb_id(sid, known_pdbs)
        topt_by_chain[cid] = row.get("predicted_topt_C", "")
        _register(cid, pid)

    for row in tm_rows:
        sid = row.get("seq_id", "")
        cid, pid = base_id(sid), seq_id_to_pdb_id(sid, known_pdbs)
        tm_by_chain[cid] = row.get("predicted_tm_C", "")
        _register(cid, pid)

    ph_by_pdb = {}
    for row in phoptnn_rows:
        pid = row.get("ID", "")
        ph_by_pdb[pid] = row.get("predicted_ph_opt", "")

    multi_chain = {pid for pid in hit_pdbs if len(chains_by_pdb.get(pid, [])) > 1}
    has_multi   = bool(multi_chain)

    print(f"[INFO] {len(multi_chain)} multi-chain PDB(s): "
          f"{sorted(multi_chain) if multi_chain else 'none'}", file=sys.stderr)

    def chain_seq_preds(cid):
        return {
            "predicted_solubility": ns_by_chain.get(cid, {}).get("predicted_solubility", "NA"),
            "predicted_usability":  ns_by_chain.get(cid, {}).get("predicted_usability",  "NA"),
            "predicted_topt_C":     topt_by_chain.get(cid, "NA"),
            "predicted_tm_C":       tm_by_chain.get(cid,  "NA"),
        }

    if not has_multi:
        headers = ["ID", "predicted_solubility", "predicted_usability",
                   "predicted_ph_opt", "predicted_topt_C", "predicted_tm_C"]
        out_path = out_dir / "all_predictions.tsv"

        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
            w.writeheader()
            for pid in hit_pdbs:
                cids = chains_by_pdb.get(pid, [])
                cid  = cids[0] if cids else None
                row  = {"ID": pid, "predicted_ph_opt": ph_by_pdb.get(pid, "NA")}
                row.update(chain_seq_preds(cid) if cid else {
                    "predicted_solubility": "NA", "predicted_usability": "NA",
                    "predicted_topt_C": "NA", "predicted_tm_C": "NA",
                })
                w.writerow(row)

        print(f"[INFO] Wrote {len(hit_pdbs)} rows → {out_path}", file=sys.stderr)
        _log_matches(hit_pdbs, ph_by_pdb, ns_by_chain, topt_by_chain,
                     tm_by_chain, chains_by_pdb)
        return

    struct_headers = ["ID", "predicted_solubility", "predicted_usability",
                      "predicted_ph_opt", "predicted_topt_C", "predicted_tm_C"]
    chain_headers  = ["chain_id", "ID", "predicted_solubility",
                      "predicted_usability", "predicted_topt_C", "predicted_tm_C"]

    struct_rows = []
    chain_rows  = []

    for pid in hit_pdbs:
        cids     = chains_by_pdb.get(pid, [])
        is_multi = pid in multi_chain

        sr = {"ID": pid, "predicted_ph_opt": ph_by_pdb.get(pid, "NA")}
        if is_multi:
            sr.update({
                "predicted_solubility": "NA", "predicted_usability": "NA",
                "predicted_topt_C": "NA",     "predicted_tm_C": "NA",
            })
        else:
            cid = cids[0] if cids else None
            sr.update(chain_seq_preds(cid) if cid else {
                "predicted_solubility": "NA", "predicted_usability": "NA",
                "predicted_topt_C": "NA", "predicted_tm_C": "NA",
            })
        struct_rows.append(sr)

        if is_multi:
            for cid in cids:
                cr = {"chain_id": cid, "ID": pid}
                cr.update(chain_seq_preds(cid))
                chain_rows.append(cr)

    struct_path = out_dir / "all_predictions_structure.tsv"
    chain_path  = out_dir / "all_predictions_chains.tsv"

    with open(struct_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=struct_headers, delimiter="\t")
        w.writeheader()
        w.writerows(struct_rows)
    print(f"[INFO] Wrote {len(struct_rows)} structure rows → {struct_path}",
          file=sys.stderr)

    with open(chain_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=chain_headers, delimiter="\t")
        w.writeheader()
        w.writerows(chain_rows)
    print(f"[INFO] Wrote {len(chain_rows)} chain rows → {chain_path}",
          file=sys.stderr)

    _log_matches(hit_pdbs, ph_by_pdb, ns_by_chain, topt_by_chain,
                 tm_by_chain, chains_by_pdb)


def _log_matches(hit_pdbs, ph_by_pdb, ns_by_chain, topt_by_chain,
                 tm_by_chain, chains_by_pdb):
    """Print a brief match summary to stderr."""
    n = len(hit_pdbs)
    ph_match = sum(1 for p in hit_pdbs if p in ph_by_pdb)
    ns_match = sum(1 for p in hit_pdbs
                   if any(c in ns_by_chain for c in chains_by_pdb.get(p, [])))
    tp_match = sum(1 for p in hit_pdbs
                   if any(c in topt_by_chain for c in chains_by_pdb.get(p, [])))
    tm_match = sum(1 for p in hit_pdbs
                   if any(c in tm_by_chain for c in chains_by_pdb.get(p, [])))
    print(f"[INFO] Match summary  pHoptNN {ph_match}/{n}  "
          f"NetSolP {ns_match}/{n}  Seq2Topt {tp_match}/{n}  Seq2Tm {tm_match}/{n}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
