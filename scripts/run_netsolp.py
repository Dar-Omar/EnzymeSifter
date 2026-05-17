#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import tempfile
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Run NetSolP for Solubility & Usability, output clean TSV."
    )
    parser.add_argument("--fasta", required=True, help="Input multi-FASTA")
    parser.add_argument("--output", required=True, help="Final TSV")
    parser.add_argument(
        "--netsolp_dir", required=True,
        help="Path to NetSolP-1.0/PredictionServer directory",
    )
    parser.add_argument(
        "--model_type", default="Distilled",
        choices=["ESM12", "ESM1b", "Distilled", "Both"],
        help="Transformer backbone (default: Distilled)",
    )
    args = parser.parse_args()

    fasta = os.path.abspath(args.fasta)
    output = os.path.abspath(args.output)
    netsolp_dir = os.path.abspath(args.netsolp_dir)
    predict_script = os.path.join(netsolp_dir, "predict.py")

    if not os.path.isfile(predict_script):
        sys.exit(f"[ERROR] predict.py not found at {predict_script}")

    os.makedirs(os.path.dirname(output), exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        tmp_csv = tmp.name

    cmd = [
        sys.executable, predict_script,
        "--FASTA_PATH", fasta,
        "--OUTPUT_PATH", tmp_csv,
        "--MODEL_TYPE", args.model_type,
        "--PREDICTION_TYPE", "SU",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, cwd=netsolp_dir, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout, file=sys.stderr)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        sys.exit(f"[ERROR] predict.py exited with code {result.returncode}")

    df = pd.read_csv(tmp_csv)
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")
    df.to_csv(output, sep="\t", index=False)

    print(f"[INFO] Wrote {len(df)} rows to {output}", file=sys.stderr)

    try:
        os.unlink(tmp_csv)
    except OSError:
        pass


if __name__ == "__main__":
    main()
