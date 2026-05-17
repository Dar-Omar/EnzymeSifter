#!/usr/bin/env python3
"""
Assign tips of a Newick tree to k clades by cutting the k-1 longest
internal branches (equivalently, binary-searching the smallest distance
threshold that produces <= k clades).
"""

import argparse
import sys
from pathlib import Path
from Bio import Phylo


def subtree_depth(clade):
    """Maximum distance from this clade to any descendant tip."""
    if clade.is_terminal():
        return 0.0
    return max(
        (child.branch_length or 0.0) + subtree_depth(child)       #imp
        for child in clade.clades
    )


def clades_at_threshold(tree, threshold):
    """
    Return list of lists — each inner list holds the tip names of one clade.
    A clade = a maximal subtree whose depth is <= threshold.
    """
    clades = []

    def walk(clade):
        if subtree_depth(clade) <= threshold:            #imp
            tips = [t.name for t in clade.get_terminals() if t.name]
            if tips:
                clades.append(tips)
        else:
            for child in clade.clades:
                walk(child)

    walk(tree.root)
    return clades


def tree_max_depth(tree):
    return subtree_depth(tree.root)


def find_k_clades(tree, k, n_iter=60):
    """
    Binary-search the smallest threshold that produces <= k clades.
    """
    n_tips = tree.count_terminals()
    if k >= n_tips:
        clades = [[t.name] for t in tree.get_terminals() if t.name]
        return clades, 0.0, len(clades)

    if k <= 1:
        clades = [[t.name for t in tree.get_terminals() if t.name]]
        return clades, tree_max_depth(tree), 1

    lo, hi = 0.0, tree_max_depth(tree)
    best = None
    for _ in range(n_iter):
        mid = (lo + hi) / 2
        clades = clades_at_threshold(tree, mid)
        n = len(clades)
        if n > k:
            lo = mid
        else:
            hi = mid
            if n == k:
                best = (clades, mid, n)
        if hi - lo < 1e-9:
            break

    if best is not None:
        return best

    final = clades_at_threshold(tree, hi)
    return final, hi, len(final)


def main():
    p = argparse.ArgumentParser(
        description="Assign tree tips to k clades."
    )
    p.add_argument("tree", help="Input Newick tree (.nwk)")
    p.add_argument("output", help="Output TSV (tip_name, clade_id)")
    p.add_argument("-k", "--clades", type=int, required=True,
                   help="Requested number of clades")
    args = p.parse_args()

    if args.clades < 1:
        sys.exit(f"[ERROR] --clades must be >= 1 (got {args.clades})")

    tree = Phylo.read(args.tree, "newick")

    for t in tree.get_terminals():
        if t.name and "|" in t.name:
            t.name = t.name.split("|")[0]

    n_tips = tree.count_terminals()
    print(f"[INFO] Tree has {n_tips} tip(s)", file=sys.stderr)

    if args.clades > n_tips:
        print(f"[WARN] Requested {args.clades} clades but tree has only "
              f"{n_tips} tip(s); producing {n_tips} singleton clades.",
              file=sys.stderr)

    clades, threshold, actual = find_k_clades(tree, args.clades)

    if actual != args.clades:
        print(f"[WARN] Could not produce exactly {args.clades} clades — "
              f"returned {actual} (threshold={threshold:.4f}).",
              file=sys.stderr)
    else:
        print(f"[INFO] Produced {actual} clade(s) at threshold {threshold:.4f}",
              file=sys.stderr)

    clades.sort(key=lambda c: (-len(c), sorted(c)[0] if c else ""))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("tip_name\tclade_id\n")
        for i, tips in enumerate(clades, 1):
            for tip in sorted(tips):
                fh.write(f"{tip}\t{i}\n")

    print(f"[INFO] Wrote clade assignments to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
