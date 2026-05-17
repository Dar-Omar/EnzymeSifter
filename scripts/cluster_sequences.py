#!/usr/bin/env python3

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from Bio import SeqIO


def count_sequences(fasta_path):
    return sum(1 for _ in SeqIO.parse(str(fasta_path), "fasta"))


def resolve_input(input_path):
    p = Path(input_path).resolve()
    if not p.exists():
        sys.exit(f"[ERROR] Input path does not exist: {p}")
    if p.is_file():
        return p, (lambda: None)
    if p.is_dir():
        fasta_files = sorted(
            f for ext in ("*.fasta", "*.fa", "*.faa")
            for f in p.glob(ext)
        )
        if not fasta_files:
            sys.exit(f"[ERROR] No .fasta/.fa/.faa files found in directory: {p}")
        print(f"[INFO] Found {len(fasta_files)} FASTA file(s) in {p}; "
              f"concatenating.", file=sys.stderr)
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, prefix="stage1_merged_"
        )
        try:
            for f in fasta_files:
                for rec in SeqIO.parse(str(f), "fasta"):
                    SeqIO.write([rec], tmp, "fasta")
        finally:
            tmp.close()
        merged_path = Path(tmp.name)
        return merged_path, (lambda: merged_path.unlink(missing_ok=True))
    sys.exit(f"[ERROR] Input path is neither a file nor a directory: {p}")


def compile_motif(motif_str):
    try:
        return re.compile(motif_str, re.IGNORECASE)
    except re.error as e:
        sys.exit(f"[ERROR] Invalid motif pattern {motif_str!r}: {e}")


def filter_by_motif(input_fasta, output_fasta, motif_report, pattern, motif_str):
    n_in = 0
    n_hit = 0
    with open(output_fasta, "w") as out_fh, open(motif_report, "w") as rpt_fh:
        rpt_fh.write("seq_id\tmatch_start\tmatch_end\tmatch_residues\n")
        for rec in SeqIO.parse(str(input_fasta), "fasta"):
            n_in += 1
            seq_upper = str(rec.seq).upper()
            m = pattern.search(seq_upper)
            if m:
                n_hit += 1
                SeqIO.write([rec], out_fh, "fasta")
                rpt_fh.write(f"{rec.id}\t{m.start() + 1}\t{m.end()}\t{m.group()}\n")
    print(f"[INFO] Motif {motif_str!r}: {n_hit}/{n_in} sequences matched",
          file=sys.stderr)
    return n_hit

def filter_by_pfam(input_fasta, output_fasta, pfam_report, pfam_arg,
                   pfam_db):
    script = Path(__file__).parent / "run_pfam_filter.py"
    cmd = [
        sys.executable, str(script),
        "--input",   str(input_fasta),
        "--output",  str(output_fasta),
        "--report",  str(pfam_report),
        "--pfam",    pfam_arg,
        "--pfam_db", str(pfam_db),
    ]
    print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"[ERROR] Pfam filter exited with code {result.returncode}")
    return count_sequences(output_fasta)

def filter_by_ec(input_fasta, output_fasta, ec_report, ec_arg, clean_dir):
    script = Path(__file__).parent / "run_clean_filter.py"
    cmd = [
        sys.executable, str(script),
        "--input",     str(input_fasta),
        "--output",    str(output_fasta),
        "--report",    str(ec_report),
        "--ec",        ec_arg,
        "--clean_dir", str(clean_dir),
    ]
    print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"[ERROR] CLEAN EC filter exited with code "
                 f"{result.returncode}")
    return count_sequences(output_fasta)


def run_mmseqs2(input_fasta, output_prefix, identity_fraction, tmp_dir):
    cmd = [
        "mmseqs", "easy-cluster",
        str(input_fasta),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id", str(identity_fraction),
        "-c", "0.8",
        "--cov-mode", "0",
        "--cluster-mode", "2",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, file=sys.stderr)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        sys.exit(f"[ERROR] mmseqs2 exited with code {result.returncode}")


def write_clustering_report(mmseqs_cluster_tsv, report_path):
    with open(mmseqs_cluster_tsv) as fin, open(report_path, "w") as fout:
        fout.write("representative\tmember\n")
        for line in fin:
            fout.write(line)


def print_summary(n_in, n_after_motif, n_after_pfam, n_after_ec, n_out,
                  identity, motif, pfam, ec,
                  nr_fasta, cluster_report, motif_report, pfam_report,
                  ec_report):
    bar = "=" * 60
    print(f"\n{bar}")
    print("  Stage 1 complete")
    print(bar)
    print(f"  Input sequences:         {n_in}")
    if motif is not None:
        print(f"  Motif filter ({motif!r}): {n_after_motif} matched")
    if pfam is not None:
        print(f"  Pfam filter ({pfam}): {n_after_pfam} matched")
    if ec is not None:
        print(f"  EC filter ({ec}): {n_after_ec} matched")
    print(f"  Non-redundant sequences: {n_out}")
    if identity is not None:
        print(f"  Identity threshold:      {identity}%")
    if (identity is None and motif is None and pfam is None and ec is None):
        print(f"  (no filters given, sequences passed through unchanged)")
    print("")
    print(f"  Non-redundant FASTA:  {nr_fasta}")
    print(f"  Clustering report:    {cluster_report}")
    if motif is not None:
        print(f"  Motif report:         {motif_report}")
    if pfam is not None:
        print(f"  Pfam report:          {pfam_report}")
    if ec is not None:
        print(f"  EC report:            {ec_report}")
    print("")
    if n_out == 0:
        print(f"  No sequences remain after filtering. Check your parameters.")
    else:
        print(f"  Next step: predict structures for the {n_out} non-redundant")
        print(f"  sequences, then run stage 2:")
        print(f"    ./run_stage2.sh /path/to/predicted_pdbs")
    print(f"{bar}\n")

def main():
    p = argparse.ArgumentParser(description="Stage 1: sequence filtering & clustering")
    p.add_argument("--input", required=True,
                   help="Input FASTA file or directory of FASTA files")
    p.add_argument("--nonredundant", required=True)
    p.add_argument("--cluster_report", required=True)
    p.add_argument("--motif_report", required=True)
    p.add_argument("--pfam_report", required=True)
    p.add_argument("--ec_report", required=True)
    p.add_argument("--identity", type=float, default=None)
    p.add_argument("--residues", type=str, default=None)
    p.add_argument("--pfam", type=str, default=None,
                   help="Comma-separated Pfam accessions (PFXXXXX)")
    p.add_argument("--pfam_db", type=str, default=None,
                   help="Path to Pfam-A.hmm (required if --pfam is given)")
    p.add_argument("--ec", type=str, default=None,
                   help="Comma-separated EC numbers; partial allowed "
                        "(e.g. 3.13.1.8 or 3.13.-.-)")
    p.add_argument("--clean_dir", type=str, default=None,
                   help="Path to cloned CLEAN repo (required if --ec is given)")
    args = p.parse_args()

    if args.pfam is not None and args.pfam_db is None:
        sys.exit("[ERROR] --pfam_db is required when --pfam is given.")
    if args.ec is not None and args.clean_dir is None:
        sys.exit("[ERROR] --clean_dir is required when --ec is given.")

    if args.residues is None and args.pfam is None and args.ec is None and args.identity is None:
        print("[WARN] No filters specified (--residues, --pfam, --ec, --identity). "
              "Stage 1 will pass sequences through unchanged. ",
              file=sys.stderr)

    input_fasta, cleanup = resolve_input(args.input)

    nr_fasta       = Path(args.nonredundant).resolve()
    cluster_report = Path(args.cluster_report).resolve()
    motif_report   = Path(args.motif_report).resolve()
    pfam_report    = Path(args.pfam_report).resolve()
    ec_report      = Path(args.ec_report).resolve()
    for pth in (nr_fasta, cluster_report, motif_report, pfam_report, ec_report):
        pth.parent.mkdir(parents=True, exist_ok=True)

    cleanups = [cleanup]

    try:
        n_in = count_sequences(input_fasta)
        if n_in == 0:
            sys.exit(f"[ERROR] No sequences found in {input_fasta}")
        print(f"[INFO] Input sequences: {n_in}", file=sys.stderr)

        working_fasta = input_fasta

        # Stage 1a: motif filter
        if args.residues is not None:
            pattern = compile_motif(args.residues)
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False, prefix="motif_filtered_"
            ) as tmp:
                motif_filtered = Path(tmp.name)
            cleanups.append(lambda p=motif_filtered: p.unlink(missing_ok=True))
            n_after_motif = filter_by_motif(
                working_fasta, motif_filtered, motif_report,
                pattern, args.residues,
            )
            working_fasta = motif_filtered
        else:
            n_after_motif = n_in
            motif_report.write_text("seq_id\tmatch_start\tmatch_end\tmatch_residues\n")

        if n_after_motif == 0 and args.residues is not None:
            nr_fasta.write_text("")
            cluster_report.write_text("representative\tmember\n")
            pfam_report.write_text("seq_id\tpfam_accession\tpfam_name\n")
            ec_report.write_text("seq_id\tmatched_ec\tall_predicted_ecs\n")
            print_summary(n_in, 0, 0, 0, 0,
                          args.identity, args.residues, args.pfam, args.ec,
                          nr_fasta, cluster_report, motif_report,
                          pfam_report, ec_report)
            return
        # Stage 1b: Pfam filter
        if args.pfam is not None:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False, prefix="pfam_filtered_"
            ) as tmp:
                pfam_filtered = Path(tmp.name)
            cleanups.append(lambda p=pfam_filtered: p.unlink(missing_ok=True))
            n_after_pfam = filter_by_pfam(
                working_fasta, pfam_filtered, pfam_report,
                args.pfam, args.pfam_db,
            )
            working_fasta = pfam_filtered
        else:
            n_after_pfam = n_after_motif
            pfam_report.write_text("seq_id\tpfam_accession\tpfam_name\n")

        if n_after_pfam == 0 and args.pfam is not None:
            nr_fasta.write_text("")
            cluster_report.write_text("representative\tmember\n")
            ec_report.write_text("seq_id\tmatched_ec\tall_predicted_ecs\n")
            print_summary(n_in, n_after_motif, 0, 0, 0,
                          args.identity, args.residues, args.pfam, args.ec,
                          nr_fasta, cluster_report, motif_report,
                          pfam_report, ec_report)
            return

        # Stage 1c: EC filter (CLEAN)
        if args.ec is not None:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False, prefix="ec_filtered_"
            ) as tmp:
                ec_filtered = Path(tmp.name)
            cleanups.append(lambda p=ec_filtered: p.unlink(missing_ok=True))
            n_after_ec = filter_by_ec(
                working_fasta, ec_filtered, ec_report,
                args.ec, args.clean_dir,
            )
            working_fasta = ec_filtered
        else:
            n_after_ec = n_after_pfam
            ec_report.write_text("seq_id\tmatched_ec\tall_predicted_ecs\n")

        if n_after_ec == 0 and args.ec is not None:
            nr_fasta.write_text("")
            cluster_report.write_text("representative\tmember\n")
            print_summary(n_in, n_after_motif, n_after_pfam, 0, 0,
                          args.identity, args.residues, args.pfam, args.ec,
                          nr_fasta, cluster_report, motif_report,
                          pfam_report, ec_report)
            return

        # Stage 1d: clustering
        if args.identity is None:
            print("[INFO] No --identity given — skipping clustering.",
                  file=sys.stderr)
            shutil.copyfile(working_fasta, nr_fasta)
            cluster_report.write_text("representative\tmember\n")
            n_out = n_after_ec
        else:
            if not (0 < args.identity <= 100):
                sys.exit(f"[ERROR] --identity must be in (0, 100] "
                         f"(got {args.identity})")
            identity_fraction = args.identity / 100.0
            print(f"[INFO] Clustering at {args.identity}% identity "
                  f"(min-seq-id={identity_fraction})", file=sys.stderr)
            with tempfile.TemporaryDirectory(prefix="mmseqs2_") as tmp:
                tmp_dir = Path(tmp)
                mmseqs_prefix = tmp_dir / "cluster"
                mmseqs_tmp = tmp_dir / "mmseqs_tmp"
                mmseqs_tmp.mkdir()
                run_mmseqs2(working_fasta, mmseqs_prefix,
                            identity_fraction, mmseqs_tmp)
                rep_fasta = Path(f"{mmseqs_prefix}_rep_seq.fasta")
                cluster_tsv = Path(f"{mmseqs_prefix}_cluster.tsv")
                if not rep_fasta.exists():
                    sys.exit(f"[ERROR] MMseqs2 output not found: {rep_fasta}")
                shutil.copyfile(rep_fasta, nr_fasta)
                write_clustering_report(cluster_tsv, cluster_report)
            n_out = count_sequences(nr_fasta)

        print_summary(n_in, n_after_motif, n_after_pfam, n_after_ec, n_out,
                      args.identity, args.residues, args.pfam, args.ec,
                      nr_fasta, cluster_report, motif_report, pfam_report, ec_report)

    finally:
        for fn in reversed(cleanups):
            try:
                fn()
            except Exception:
                pass


if __name__ == "__main__":
    main()
