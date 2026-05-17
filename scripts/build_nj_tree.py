#!/usr/bin/env python3

import sys
import argparse
from pathlib import Path

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import write as write_tree


def main():
    parser = argparse.ArgumentParser(
        description="Build a Neighbor-Joining tree from an aligned FASTA."
    )
    parser.add_argument("alignment", help="Input aligned FASTA (.afa)")
    parser.add_argument("output", help="Output Newick tree (.nwk)")
    parser.add_argument(
        "--model",
        default="blosum62",
        help="Distance model for DistanceCalculator (default: blosum62)",
    )
    args = parser.parse_args()

    aln = AlignIO.read(args.alignment, "fasta")
    print(f"[INFO] Loaded alignment: {len(aln)} sequences, {aln.get_alignment_length()} columns", file=sys.stderr)

    if len(aln) < 3:
        print("[WARN] Fewer than 3 sequences.", file=sys.stderr)

    calculator = DistanceCalculator(args.model)
    dm = calculator.get_distance(aln)

    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)
    n_neg = 0
    for clade in nj_tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
            n_neg += 1
    if n_neg:
        print(f"[INFO] Clamped {n_neg} negative branch length(s) to 0", file=sys.stderr)
    nj_tree.root_at_midpoint()

    nj_tree.ladderize()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    for clade in nj_tree.get_terminals():
        if clade.name and "|" in clade.name:
            clade.name = clade.name.split("|")[0]
    write_tree(nj_tree, args.output, "newick")

    print(f"[INFO] NJ tree written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
