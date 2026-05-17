#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from Bio import SeqIO

MAX_ESM1B_LEN = 1022

def parse_ec_specs(ec_str):
    specs = []
    for raw in ec_str.split(","):
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split(".")
        while parts and parts[-1] == "-":
            parts.pop()
        if not parts:
            sys.exit(f"[ERROR] EC spec is empty after stripping placeholders: {raw!r}")
        if len(parts) > 4:
            sys.exit(f"[ERROR] EC spec {raw!r} has more than 4 levels.")
        for seg in parts:
            if not seg.isdigit():
                sys.exit(f"[ERROR] Invalid EC segment {seg!r} in {raw!r}; "
                         f"expected digits or '-'.")
        specs.append(tuple(parts))
    if not specs:
        sys.exit("[ERROR] --ec must contain at least one EC number.")
    return specs


def ec_matches(predicted_ec, specs):
    pred_parts = predicted_ec.split(".")
    for spec in specs:
        if len(pred_parts) < len(spec):
            continue
        if all(pred_parts[i] == spec[i] for i in range(len(spec))):
            return True
    return False


def parse_clean_output(csv_path):
    preds = {}
    with open(csv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            seq_id = parts[0]
            ec_pairs = []
            for tok in parts[1:]:
                tok = tok.strip()
                if not tok or not tok.startswith("EC:"):
                    continue
                body = tok[3:]
                if "/" in body:
                    ec, dist_str = body.split("/", 1)
                    try:
                        dist = float(dist_str)
                    except ValueError:
                        dist = float("nan")
                else:
                    ec, dist = body, float("nan")
                ec_pairs.append((ec, dist))
            preds[seq_id] = ec_pairs
    return preds


def run_clean(clean_dir, fasta_basename):
    app_dir = Path(clean_dir) / "app"
    infer_py = app_dir / "CLEAN_infer_fasta.py"
    if not infer_py.exists():
        sys.exit(f"[ERROR] CLEAN_infer_fasta.py not found at {infer_py}")
    cmd = [sys.executable, str(infer_py), "--fasta_data", fasta_basename]
    print(f"[INFO] Running CLEAN: {' '.join(cmd)} (cwd={app_dir})",
          file=sys.stderr)
    result = subprocess.run(cmd, cwd=str(app_dir),
                            capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            sys.stderr.write(result.stdout)
        if result.stderr:
            sys.stderr.write(result.stderr)
        sys.exit(f"[ERROR] CLEAN_infer_fasta.py exited with code "
                 f"{result.returncode}")

def main():
    p = argparse.ArgumentParser(
        description="Filter sequences by predicted EC number using CLEAN."
    )
    p.add_argument("--input",     required=True, help="Input FASTA")
    p.add_argument("--output",    required=True,
                   help="Output FASTA (sequences whose predicted EC matches)")
    p.add_argument("--report",    required=True, help="Output EC report TSV")
    p.add_argument("--ec",        required=True,
                   help="Comma-separated EC numbers; partial allowed "
                        "(e.g. 3.13.1.8 or 3.13.-.- or 3.13)")
    p.add_argument("--clean_dir", required=True,
                   help="Path to cloned CLEAN repo (must contain app/ and "
                        "have been set up via setup_clean.sh)")
    args = p.parse_args()

    specs = parse_ec_specs(args.ec)
    pretty_specs = ", ".join(".".join(s) for s in specs)
    print(f"[INFO] Filtering for EC specs: {pretty_specs}", file=sys.stderr)

    input_fasta  = Path(args.input).resolve()
    output_fasta = Path(args.output).resolve()
    report_tsv   = Path(args.report).resolve()
    clean_dir    = Path(args.clean_dir).resolve()
    app_dir      = clean_dir / "app"

    if not input_fasta.exists():
        sys.exit(f"[ERROR] Input FASTA not found: {input_fasta}")
    if not (clean_dir / ".setup_done").exists():
        sys.exit(f"[ERROR] CLEAN not set up at {clean_dir} "
                 f"(run setup_clean.sh first).")

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    staged_basename = f"stage1_clean_{os.getpid()}_{int(time.time() * 1000) % 10**9}"
    staged_fasta    = app_dir / "data" / "inputs" / f"{staged_basename}.fasta"
    result_csv      = app_dir / "results" / "inputs" / f"{staged_basename}_maxsep.csv"

    staged_fasta.parent.mkdir(parents=True, exist_ok=True)
    n_staged = 0
    skipped_long = []
    with open(staged_fasta, "w") as out_fh:
        for rec in SeqIO.parse(str(input_fasta), "fasta"):
            seq_str = str(rec.seq)
            if len(seq_str) > MAX_ESM1B_LEN:
                skipped_long.append((rec.id, len(seq_str)))
                continue
            out_fh.write(f">{rec.id}\n{seq_str}\n")
            n_staged += 1
    print(f"[INFO] Staged {n_staged} sequence(s) to {staged_fasta} "
          f"(headers normalised to id-only).", file=sys.stderr)
    if skipped_long:
        print(f"[WARN] Skipped {len(skipped_long)} sequence(s) longer than "
              f"{MAX_ESM1B_LEN} aa (ESM-1b's maximum input length); these "
              f"cannot receive a CLEAN prediction:", file=sys.stderr)
        for sid, slen in skipped_long:
            print(f"[WARN]   {sid}  ({slen} aa)", file=sys.stderr)
    if n_staged == 0:
        sys.exit("[ERROR] No sequences remain after length filtering; "
                 "all input sequences exceed ESM-1b's 1022 aa limit.")
    try:
        run_clean(clean_dir, staged_basename)

        if not result_csv.exists():
            sys.exit(f"[ERROR] CLEAN result not found at {result_csv}")

        preds = parse_clean_output(result_csv)
        print(f"[INFO] CLEAN produced predictions for {len(preds)} sequence(s).",
              file=sys.stderr)

        skipped_ids = {sid for sid, _ in skipped_long}
        n_in = 0
        n_hit = 0
        with open(output_fasta, "w") as out_fh, \
             open(report_tsv,   "w") as rpt_fh:
            rpt_fh.write("seq_id\tmatched_ec\tall_predicted_ecs\n")
            for rec in SeqIO.parse(str(input_fasta), "fasta"):
                n_in += 1
                if rec.id in skipped_ids:
                    rpt_fh.write(f"{rec.id}\tSKIPPED_TOO_LONG\tNA\n")
                    continue
                rec_preds = preds.get(rec.id, [])
                matches = [(ec, d) for ec, d in rec_preds
                           if ec_matches(ec, specs)]
                if matches:
                    n_hit += 1
                    SeqIO.write([rec], out_fh, "fasta")
                    all_ecs = ";".join(
                        f"{ec}/{d:.4f}" if d == d else f"{ec}/NA"
                        for ec, d in rec_preds
                    )
                    for ec, d in matches:
                        d_str = f"{d:.4f}" if d == d else "NA"
                        rpt_fh.write(f"{rec.id}\t{ec}/{d_str}\t{all_ecs}\n")
        print(f"[INFO] CLEAN EC filter: {n_hit}/{n_in} sequences matched "
              f"(specs: {pretty_specs})", file=sys.stderr)
    finally:
        for path in (staged_fasta, result_csv):
            try:
                path.unlink(missing_ok=True)
            except Exception:
                pass


if __name__ == "__main__":
    main()
