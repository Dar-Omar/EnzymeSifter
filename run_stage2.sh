#!/bin/bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]

Filter options (all optional):
  -solubility <min>        Minimum predicted solubility   (e.g. 0.4)
  -usability  <min>        Minimum predicted usability    (e.g. 0.35)
  -tm  <min>  | <lo:hi>    Melting temp: cutoff OR interval
                             e.g. -tm 55         → Tm >= 55, higher better
                             e.g. -tm 50:80      → 50 <= Tm <= 80, mid best
  -topt  <lo:hi>           Optimal temp interval (required form)
                             e.g. -topt 35:42
  -phopt <lo:hi>           Optimal pH interval (required form)
                             e.g. -phopt 7.0:8.5

Clade options (optional):
  -clades <n>              Partition the NJ tree into <n> clades.
                           If combined with filter options, writes
                           predictions_output/clade_representatives.tsv
                           with the best-scoring enzyme per clade.

Scoring (when filters are given):
  Each passing enzyme receives a score in [0, 1]:
    cutoff properties:   normalised across the passing pool (higher=better)
    interval properties: closeness to the interval midpoint
  Combined score = mean of per-property scores over user-specified
  properties with non-NA values for that enzyme.

Examples:
  ./run_stage2.sh /data/pdbs
  ./run_stage2.sh /data/pdbs -usability 0.35 -solubility 0.4
  ./run_stage2.sh /data/pdbs -tm 55
  ./run_stage2.sh /data/pdbs -tm 50:80 -topt 35:42 -phopt 7:8
  ./run_stage2.sh /data/pdbs -clades 7
  ./run_stage2.sh /data/pdbs -usability 0.35 -tm 50:80 -topt 35:42 -clades 10
EOF
    exit 1
}

if [ -z "${1:-}" ]; then
    usage
fi

PDB_DIR="$1"
shift

# Parse optional filter flags
FILTER_SOLUBILITY=""
FILTER_USABILITY=""
FILTER_PHOPT=""
FILTER_TOPT=""
FILTER_TM=""
CLADES=""

while [ $# -gt 0 ]; do
    case "$1" in
        -solubility)
            [ -z "${2:-}" ] && echo "Error: -solubility requires a value" && exit 1
            FILTER_SOLUBILITY="$2"
            shift 2 ;;
        -usability)
            [ -z "${2:-}" ] && echo "Error: -usability requires a value" && exit 1
            FILTER_USABILITY="$2"
            shift 2 ;;
        -phopt)
            [ -z "${2:-}" ] && echo "Error: -phopt requires a value" && exit 1
            FILTER_PHOPT="$2"
            shift 2 ;;
        -topt)
            [ -z "${2:-}" ] && echo "Error: -topt requires a value" && exit 1
            FILTER_TOPT="$2"
            shift 2 ;;
        -tm)
            [ -z "${2:-}" ] && echo "Error: -tm requires a value" && exit 1
            FILTER_TM="$2"
            shift 2 ;;
        -clades)
            [ -z "${2:-}" ] && echo "Error: -clades requires a value" && exit 1
            CLADES="$2"
            shift 2 ;;
        *)
            echo "Unknown option: $1"
            usage ;;
    esac
done

CONFIG="pdb_dir=${PDB_DIR}"
[ -n "$FILTER_SOLUBILITY" ] && CONFIG+=" filter_solubility=${FILTER_SOLUBILITY}"
[ -n "$FILTER_USABILITY" ]  && CONFIG+=" filter_usability=${FILTER_USABILITY}"
[ -n "$FILTER_PHOPT" ]      && CONFIG+=" filter_phopt=${FILTER_PHOPT}"
[ -n "$FILTER_TOPT" ]       && CONFIG+=" filter_topt=${FILTER_TOPT}"
[ -n "$FILTER_TM" ]         && CONFIG+=" filter_tm=${FILTER_TM}"
[ -n "$CLADES" ]            && CONFIG+=" clades=${CLADES}"

snakemake --snakefile Snakefile_stage2 --config ${CONFIG} --use-conda -j "$(nproc)" \
    --quiet rules progress
