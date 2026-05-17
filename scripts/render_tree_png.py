#!/usr/bin/env python3

import sys
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo


def load_clade_assignments(path):
    """Return dict: tip_name -> clade_id (int)."""
    mapping = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row.get("tip_name", "").strip()
            cid  = row.get("clade_id", "").strip()
            if name and cid:
                mapping[name] = int(cid)
    return mapping


def load_representatives(path):
    """Return set of id strings from clade_representatives.tsv."""
    reps = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = row.get("ID", "").strip()
            if pid and pid != "NA":
                reps.add(pid)
    return reps


def tip_to_pdb(tip_name, known_pdbs=None):
    """
    Map a tree tip back to its pdb_id.

    If `known_pdbs` is given and the tip is itself a known pdb, return
    it unchanged (single-chain case — pdb filename may end in _<digits>).
    Otherwise fall back to stripping a trailing '_<digits>' suffix
    (multi-chain tip of the form pdb_id_<N>).
    """
    if known_pdbs is not None and tip_name in known_pdbs:
        return tip_name
    parts = tip_name.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return tip_name


def build_clade_palette(n_clades):
    """Pick a categorical palette sized for the number of clades."""
    if n_clades <= 10:
        return list(plt.get_cmap("tab10").colors)
    if n_clades <= 20:
        return list(plt.get_cmap("tab20").colors)
    # >20 clades: cycle tab20 (rare; warn upstream)
    base = list(plt.get_cmap("tab20").colors)
    return [base[i % len(base)] for i in range(n_clades)]


def main():
    parser = argparse.ArgumentParser(
        description="Render a Newick tree as PNG."
    )
    parser.add_argument("tree", help="Input Newick tree (.nwk)")
    parser.add_argument("output", help="Output PNG file")
    parser.add_argument("--width", type=float, default=12)
    parser.add_argument("--height", type=float, default=8)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--highlight", default=None,
                        help="Text file with one PDB ID per line (green)")
    parser.add_argument("--clades", default=None,
                        help="TSV with columns tip_name, clade_id "
                             "— colours tips by clade")
    parser.add_argument("--representatives", default=None,
                        help="clade_representatives.tsv — marks the "
                             "representative PDB of each clade with a "
                             "gold star at the tip")
    parser.add_argument("--hit_list", default=None,
                        help="data/enzymm/hit_pdbs.txt — disambiguates "
                             "single-chain PDBs whose filename ends in "
                             "_<digits>")
    args = parser.parse_args()

    known_pdbs = None
    if args.hit_list:
        with open(args.hit_list) as fh:
            known_pdbs = {ln.strip() for ln in fh if ln.strip()}
        print(f"[INFO] Loaded {len(known_pdbs)} known pdb_id(s) from "
              f"{args.hit_list}", file=sys.stderr)

    tree = Phylo.read(args.tree, "newick")

    for clade in tree.get_terminals():
        if clade.name and "|" in clade.name:
            clade.name = clade.name.split("|")[0]
    for clade in tree.get_nonterminals():
        clade.name = None

    label_colors = {}

    # --- Clade colouring (takes precedence on colour) ---
    n_clades = 0
    if args.clades:
        tip_to_clade = load_clade_assignments(args.clades)
        unique = sorted(set(tip_to_clade.values()))
        n_clades = len(unique)
        palette = build_clade_palette(n_clades)
        clade_to_color = {c: palette[i] for i, c in enumerate(unique)}
        for clade in tree.get_terminals():
            if clade.name in tip_to_clade:
                label_colors[clade.name] = clade_to_color[tip_to_clade[clade.name]]
        print(f"[INFO] Coloured {len(tip_to_clade)} tip(s) across {n_clades} clade(s)",
              file=sys.stderr)

    # --- Highlight mode (only if --clades not used) ---
    elif args.highlight:
        highlight_ids = set()
        with open(args.highlight) as fh:
            for line in fh:
                if line.strip():
                    highlight_ids.add(line.strip())
        print(f"[INFO] Highlighting {len(highlight_ids)} enzyme(s) in green",
              file=sys.stderr)
        for clade in tree.get_terminals():
            if clade.name:
                matched = any(hid in clade.name for hid in highlight_ids)
                label_colors[clade.name] = "green" if matched else "black"

    # --- Representatives (optional overlay) ---
    representative_ids = set()
    if args.representatives:
        representative_ids = load_representatives(args.representatives)
        print(f"[INFO] Loaded {len(representative_ids)} representative PDB(s) "
              f"to mark with stars", file=sys.stderr)

    n_tips = len(tree.get_terminals())
    height = max(args.height, n_tips * 0.3)

    fig, ax = plt.subplots(figsize=(args.width, height))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_colors=label_colors if label_colors else None)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # --- Star markers for representatives ---
    # Done after Phylo.draw so we can read the final tip-label positions.
    # Stars are placed just *after* each representative's label (to its
    # right), using the measured label bounding box — this guarantees no
    # overlap regardless of figure size, dpi, or label length.
    n_marked = 0
    if representative_ids:
        # Force a draw so text extents are available from the renderer.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        x_lim = ax.get_xlim()
        gap = (x_lim[1] - x_lim[0]) * 0.01  # small data-coord gap after label

        star_xs, star_ys, star_cs = [], [], []
        for text in ax.texts:
            label = text.get_text().strip()
            if not label:
                continue
            if tip_to_pdb(label, known_pdbs) in representative_ids:
                # Measure the label's bbox in display coords, convert to data
                # coords so we can place the star at (xmax_of_label + gap).
                bbox_disp = text.get_window_extent(renderer=renderer)
                bbox_data = bbox_disp.transformed(ax.transData.inverted())
                _, y_text = text.get_position()
                star_xs.append(bbox_data.xmax + gap)
                star_ys.append(y_text)
                # Inherit the tip's clade colour; fall back to gold if the
                # representative somehow isn't in a coloured clade.
                star_cs.append(label_colors.get(label, "gold"))

        if star_xs:
            ax.scatter(
                star_xs, star_ys,
                marker="*",
                s=220,
                c=star_cs,
                edgecolors="black",
                linewidths=0.9,
                zorder=5,
                clip_on=False,   # allow stars to sit past the axes xlim
            )
            n_marked = len(star_xs)
            print(f"[INFO] Marked {n_marked} tip(s) as clade representatives",
                  file=sys.stderr)
        else:
            print("[WARN] No tips matched the representative PDB list — "
                  "no stars drawn.", file=sys.stderr)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    print(f"[INFO] Tree PNG written to {args.output} ({args.dpi} dpi)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
