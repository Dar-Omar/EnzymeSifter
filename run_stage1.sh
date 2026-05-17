#!/bin/bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: ./run_stage1.sh <input> [options]

Stage 1: filter and cluster input amino-acid sequences, producing a
non-redundant FASTA for structure prediction.

<input> can be either:
  - a single FASTA file (.fasta/.fa/.faa), or
  - a directory containing one or more FASTA files (concatenated).

Options:
  -identity <pct>     Cluster at <pct>% identity (0-100), keeping one
                      representative per cluster.
  -residues <motif>   Keep only sequences matching a motif.
                      Regex syntax: '.' = any single residue.
                      Examples: G.S.G, SHD
  -pfam <IDs>         Keep only sequences that hit any of the given
                      Pfam accessions (comma-separated, OR logic).
                      Uses Pfam gathering thresholds via hmmsearch.
                      Examples: PF07519
                                PF07519,PF00657
                      First run auto-downloads Pfam (~1.5 GB).

  -ec <IDs>           Keep only sequences whose CLEAN-predicted EC
                      number matches any of the given EC specs
                      (comma-separated, OR logic). Partial specs
                      allowed: a missing or '-' digit is a wildcard.
                      Examples: 3.13.1.8
                                3.13.-.-           (same as 3.13)
                                3.13.1.8,2.5.1.94
                      First run clones CLEAN, builds the package, and
                      downloads pretrained weights from Google Drive
                      (~few hundred MB) plus ESM-1b on first inference
                      (~7 GB to ~/.cache/torch).

When combined, filters apply in this order:
  residues -> pfam -> ec -> identity clustering

Examples:
  ./run_stage1.sh seqs.fasta -pfam PF07519
  ./run_stage1.sh seqs.fasta -residues G.S.G -ec 3.13.1.8,2.5.1.94
  ./run_stage1.sh seqs.fasta -pfam PF07519 -ec 2.5.1.94 -identity 90


Outputs:
  data/stage1/nonredundant.fasta       Non-redundant sequences
  data/stage1/clustering_report.tsv    Representative-to-member mapping
  data/stage1/motif_report.tsv         motif match results
  data/stage1/pfam_report.tsv          Pfam match results
  data/stage1/ec_report.tsv            EC prediction results
EOF
    exit 1
}

if [ -z "${1:-}" ]; then
    usage
fi

INPUT_FASTA="$1"
shift

IDENTITY=""
RESIDUES=""
PFAM=""
EC=""

while [ $# -gt 0 ]; do
    case "$1" in
        -identity)
            [ -z "${2:-}" ] && echo "Error: -identity requires a value" && exit 1
            IDENTITY="$2"
            shift 2 ;;
        -residues)
            [ -z "${2:-}" ] && echo "Error: -residues requires a value" && exit 1
            RESIDUES="$2"
            shift 2 ;;
        -pfam)
            [ -z "${2:-}" ] && echo "Error: -pfam requires a value" && exit 1
            PFAM="$2"
            shift 2 ;;
        -ec)
            [ -z "${2:-}" ] && echo "Error: -ec requires a value" && exit 1
            EC="$2"
            shift 2 ;;
        *)
            echo "Unknown option: $1"
            usage ;;
    esac
done

CONFIG="input_fasta=${INPUT_FASTA}"
[ -n "$IDENTITY" ] && CONFIG+=" identity=${IDENTITY}"
[ -n "$RESIDUES" ] && CONFIG+=" residues=${RESIDUES}"
[ -n "$PFAM" ]     && CONFIG+=" pfam=${PFAM}"
[ -n "$EC" ]       && CONFIG+=" ec=${EC}"

snakemake --snakefile Snakefile_stage1 --config ${CONFIG} --use-conda -j "$(nproc)" \
    --quiet rules progress
