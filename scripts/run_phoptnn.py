#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import tempfile
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser(
        description="Run pHoptNN on EnzyMM-hit PDBs and output a TSV."
    )
    p.add_argument("--pdb_dir", required=True,
                   help="Directory containing all input PDB files")
    p.add_argument("--hit_list", required=True,
                   help="Text file with one PDB basename per line (no extension)")
    p.add_argument("--phoptnn_dir", required=True,
                   help="Path to cloned pHoptNN repository")
    p.add_argument("--output", required=True,
                   help="Output TSV path")
    p.add_argument("--quiet", action="store_true",
                   help="Silence pHoptNN output")
    args = p.parse_args()

    pdb_dir = Path(args.pdb_dir).resolve()
    phoptnn_dir = Path(args.phoptnn_dir).resolve()
    output = Path(args.output).resolve()
    interface_py = phoptnn_dir / "phoptnn_interface.py"

    if not interface_py.exists():
        sys.exit(f"[ERROR] phoptnn_interface.py not found at {interface_py}")

    with open(args.hit_list) as fh:
        hits = [line.strip() for line in fh if line.strip()]

    if not hits:
        print("[WARN] No hits in hit_list — writing empty TSV.", file=sys.stderr)
        pd.DataFrame(columns=["ID", "predicted_ph_opt"]).to_csv(
            output, sep="\t", index=False
        )
        return

    with tempfile.TemporaryDirectory(prefix="phoptnn_input_") as tmpdir:
        tmp_input = Path(tmpdir) / "pdbs"
        tmp_input.mkdir()
        tmp_output = Path(tmpdir) / "pred_out"
        tmp_pqr = Path(tmpdir) / "pqr_files"

        linked = 0
        for pdb_name in hits:
            src = pdb_dir / f"{pdb_name}.pdb"
            if src.exists():
                dst = tmp_input / f"{pdb_name}.pdb"
                dst.symlink_to(src)
                linked += 1
            else:
                print(f"[WARN] PDB not found, skipping: {src}", file=sys.stderr)

        if linked == 0:
            sys.exit("[ERROR] No PDB files could be linked for pHoptNN.")

        print(f"[INFO] Running pHoptNN on {linked} PDB(s)...", file=sys.stderr)

        cmd = [
            sys.executable, str(interface_py),
            str(tmp_input),
            "--save_dir", str(tmp_output),
            "--pqr_dir", str(tmp_pqr),
        ]
        if args.quiet:
            cmd.append("--quiet")

        subprocess.run(cmd, check=True, cwd=str(phoptnn_dir))

        pred_csv = tmp_output / "predictions.csv"
        if not pred_csv.exists():
            candidates = list(tmp_output.glob("*.csv"))
            if candidates:
                pred_csv = candidates[0]
            else:
                sys.exit(f"[ERROR] No prediction CSV found in {tmp_output}")

        df = pd.read_csv(pred_csv)
        df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")

        id_col = None
        for candidate in ["protein", "name", "file", "pdb", "input", "id"]:
            if candidate in df.columns:
                id_col = candidate
                break

        if id_col is None:
            id_col = df.columns[0]

        df["ID"] = df[id_col].apply(
            lambda x: Path(str(x)).stem
        )

        ph_col = None
        for candidate in ["predicted_ph", "ph_opt", "ph_optimum",
                          "prediction", "predicted_ph_opt", "ph_optimum_pred",
                          "pred", "ph"]:
            if candidate in df.columns:
                ph_col = candidate
                break

        if ph_col is None:
            numeric_cols = df.select_dtypes(include="number").columns.tolist()
            if numeric_cols:
                ph_col = numeric_cols[-1]
            else:
                sys.exit("[ERROR] Cannot identify pH prediction column in output.")

        out_df = df[["ID", ph_col]].copy()
        out_df.columns = ["ID", "predicted_ph_opt"]

    output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output, sep="\t", index=False)
    print(f"[INFO] Wrote {len(out_df)} predictions to {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
