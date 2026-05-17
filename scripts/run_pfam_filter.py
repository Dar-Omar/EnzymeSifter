#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from Bio import SeqIO


PFAM_ACCESSION_RE = re.compile(r"^PF\d{5}$")


def validate_accessions(accession_str):
    accessions = [a.strip() for a in accession_str.split(",") if a.strip()]
    if not accessions:
        sys.exit("[ERROR] --pfam must be a non-empty list of accessions.")
    for acc in accessions:
        if not PFAM_ACCESSION_RE.match(acc):
            sys.exit(f"[ERROR] Invalid Pfam accession {acc!r}. "
                     f"Expected format PFXXXXX (five digits).")
    return accessions

def fetch_hmm_subset(pfam_db, accessions, subset_hmm):
    """Extract only the wanted HMMs into a small temp file."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".keys", delete=False, prefix="hmmfetch_keys_"
    ) as kf:
        for acc in accessions:
            kf.write(acc + "\n")
        keys_path = Path(kf.name)
    try:
        cmd = ["hmmfetch", "-f", "-o", str(subset_hmm),
               str(pfam_db), str(keys_path)]
        print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        if result.returncode != 0:
            sys.exit(f"[ERROR] hmmfetch exited with code {result.returncode}")
    finally:
        keys_path.unlink(missing_ok=True)


def run_hmmsearch(fasta, hmm_file, out_tsv, cpu):
    cmd = [
        "hmmsearch",
        "--cpu", str(cpu),
        "--cut_ga",
        "--tblout", str(out_tsv),
        "--noali",
        str(hmm_file),
        str(fasta),
    ]
    print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        sys.exit(f"[ERROR] hmmsearch exited with code {result.returncode}")

def parse_tblout_hmmsearch(tblout_path, wanted_accessions):
    hits_by_seq = {}
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            seq_id      = parts[0]
            hmm_name    = parts[2]
            hmm_acc     = parts[3]
            try:
                evalue = float(parts[4])
            except ValueError:
                continue
            bare_acc = hmm_acc.split(".")[0]
            if bare_acc in wanted_accessions:
                hits_by_seq.setdefault(seq_id, []).append(
                    (bare_acc, hmm_name, evalue)
                )
    return hits_by_seq

def main():
    p = argparse.ArgumentParser(description="Filter sequences by Pfam accession")
    p.add_argument("--input",    required=True, help="Input FASTA")
    p.add_argument("--output",   required=True, help="Output FASTA (matched sequences)")
    p.add_argument("--report",   required=True, help="Output Pfam report TSV")
    p.add_argument("--pfam",     required=True, help="Comma-separated Pfam accessions (PFXXXXX)")
    p.add_argument("--pfam_db",  required=True, help="Path to Pfam-A.hmm")

    default_cpu = min(16, os.cpu_count() or 1)
    p.add_argument("--cpu", type=int, default=default_cpu,
                   help=f"hmmsearch worker threads "
                        f"(default: {default_cpu}, auto-detected)")

    args = p.parse_args()

    accessions = validate_accessions(args.pfam)
    wanted = set(accessions)
    print(f"[INFO] Filtering for Pfam accessions: {', '.join(accessions)}",
          file=sys.stderr)

    input_fasta  = Path(args.input).resolve()
    output_fasta = Path(args.output).resolve()
    report_tsv   = Path(args.report).resolve()
    pfam_db      = Path(args.pfam_db).resolve()

    for pth in (output_fasta, report_tsv):
        pth.parent.mkdir(parents=True, exist_ok=True)

    if not pfam_db.exists():
        sys.exit(f"[ERROR] Pfam HMM file not found: {pfam_db}")

    ssi = Path(f"{pfam_db}.ssi")
    if not ssi.exists():
        sys.exit(f"[ERROR] hmmfetch index missing: {ssi}\n"
             f"        Re-run setup_pfam.sh.")

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tbl", delete=False, prefix="hmmsearch_"
    ) as tmp:
        tblout_path = Path(tmp.name)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".hmm", delete=False, prefix="pfam_subset_"
    ) as tmp:
        subset_hmm = Path(tmp.name)

    try:
        fetch_hmm_subset(pfam_db, accessions, subset_hmm)
        run_hmmsearch(input_fasta, subset_hmm, tblout_path, args.cpu)
        hits_by_seq = parse_tblout_hmmsearch(tblout_path, wanted)
    finally:
        tblout_path.unlink(missing_ok=True)
        subset_hmm.unlink(missing_ok=True)

    n_in = 0
    n_hit = 0
    with open(output_fasta, "w") as out_fh, open(report_tsv, "w") as rpt_fh:
        rpt_fh.write("seq_id\tpfam_accession\tpfam_name\n")
        for rec in SeqIO.parse(str(input_fasta), "fasta"):
            n_in += 1
            hits = hits_by_seq.get(rec.id, [])
            if hits:
                n_hit += 1
                SeqIO.write([rec], out_fh, "fasta")
                for acc, name, ev in hits:
                    rpt_fh.write(f"{rec.id}\t{acc}\t{name}\n")
    print(f"[INFO] Pfam filter: {n_hit}/{n_in} sequences matched",
          file=sys.stderr)


if __name__ == "__main__":
    main()
