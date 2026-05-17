# EnzymeSifter

A two-stage Snakemake pipeline for sifting through large enzyme sequence collections and identifying promising candidates for downstream characterisation — from raw FASTA all the way to ranked, clade-representative enzymes annotated with predicted solubility, optimal pH, optimal temperature, and melting temperature.

EnzymeSifter wraps a number of state-of-the-art tools (CLEAN, EnzyMM, NetSolP, pHoptNN, Seq2Topt, MMseqs2, HMMER, MUSCLE) behind a single command-line interface and handles all of their setup, environment management, and inter-tool data plumbing for you.

---

## Table of contents

- [Overview](#overview)
- [Pipeline architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Stage 1 — sequence filtering & clustering](#stage-1--sequence-filtering--clustering)
- [Between the stages — structure prediction](#between-the-stages--structure-prediction)
- [Stage 2 — structural & biophysical screening](#stage-2--structural--biophysical-screening)
- [Output structure](#output-structure)
- [External tools](#external-tools)
- [Disk-space & runtime expectations](#disk-space--runtime-expectations)
- [Troubleshooting](#troubleshooting)
- [Citing EnzymeSifter](#citing-enzymesifter)
- [License](#license)

---

## Overview

EnzymeSifter is built for the common bioinformatics scenario where you have:

- a large collection of candidate enzyme sequences (anywhere from hundreds to hundreds of thousands), and
- a need to triage them down to a tractable, diverse, biophysically promising shortlist before committing to wet-lab work or expensive structural characterisation.

The pipeline runs in two stages with a structure-prediction step (carried out by the user, externally) in between:

1. **Stage 1 — sequence-level triage.** Filter by catalytic-residue motif, Pfam domain, and/or CLEAN-predicted EC number, then cluster at a user-defined identity threshold using MMseqs2.
2. **(User step) Structure prediction.** Generate PDB structures for the Stage 1 representatives using a tool of your choice (AlphaFold, ESMFold, ColabFold, etc.).
3. **Stage 2 — structural & biophysical screening.** Confirm enzymatic activity with EnzyMM, predict solubility/usability (NetSolP), optimal pH (pHoptNN), and optimal/melting temperatures (Seq2Topt/Seq2Tm), build an NJ phylogenetic tree, optionally partition it into clades, and select the best-scoring representative per clade according to your filter criteria.

---

## Pipeline architecture

```
                          ┌────────────────────────────────┐
                          │   Input FASTA (or directory)   │
                          └───────────────┬────────────────┘
                                          │
                  ╔═══════════════════════▼═══════════════════════╗
                  ║                    STAGE 1                    ║
                  ║  motif → Pfam (HMMER) → EC (CLEAN) → MMseqs2  ║
                  ╚═══════════════════════╤═══════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │   data/stage1/nonredundant.fa  │
                          └───────────────┬────────────────┘
                                          │
                  ╔═══════════════════════▼═══════════════════════╗
                  ║    (user) predict 3D structures externally    ║
                  ║       AlphaFold / ESMFold / ColabFold …       ║
                  ╚═══════════════════════╤═══════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │     directory of .pdb files    │
                          └───────────────┬────────────────┘
                                          │
   ╔══════════════════════════════════════▼══════════════════════════════════════╗
   ║                                  STAGE 2                                    ║
   ║   EnzyMM ── filter to enzymatic hits ──┬── NetSolP   (solubility, usable)  ║
   ║                                        ├── pHoptNN   (pH optimum)          ║
   ║                                        ├── Seq2Topt  (T optimum)           ║
   ║                                        ├── Seq2Tm    (T melting)           ║
   ║                                        └── MUSCLE ─► NJ tree ─► clades ──► ║
   ║                                                       representatives      ║
   ╚══════════════════════════════════════╤══════════════════════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │  predictions_output/*.tsv      │
                          │  data/trees/nj_tree_*.png      │
                          └────────────────────────────────┘
```

---

## Requirements

- Linux
- [Conda](https://docs.conda.io/en/latest/miniconda.html)
- [Snakemake](https://snakemake.readthedocs.io/)
- CPU is sufficient for all steps.

EnzymeSifter creates and manages its own conda environments automatically via `snakemake --use-conda`. You do not need to install any of the underlying tools yourself.

---

## Installation

```bash
git clone https://github.com/<your-user>/EnzymeSifter.git
cd EnzymeSifter

# Set up the controller environment (snakemake + conda frontend)
conda create -n enzymesifter -c conda-forge -c bioconda \
    snakemake mamba -y
conda activate enzymesifter

chmod +x run_stage1.sh run_stage2.sh scripts/*.sh
```

That is all you need to do. On the first invocation of either stage, Snakemake will build all of the required tool-specific conda environments under `.snakemake/conda/`, and the relevant setup scripts (`scripts/setup_*.sh`) will fetch external databases and model weights into `external/`.

---

## Quickstart

Sift a FASTA of candidate sequences, predict structures yourself, then run structural/biophysical screening:

```bash
# Stage 1: pull out the Pfam PF07519 hits, dereplicate at 90% identity
./run_stage1.sh test_input/my_sequences.fasta -pfam PF07519 -identity 90

# (predict structures for data/stage1/nonredundant.fasta using
#  your structure-prediction tool of choice — place the resulting
#  .pdb files into ./predicted_pdbs/ )

# Stage 2: confirm enzymatic activity, predict properties,
#          partition the resulting NJ tree into 8 clades, and
#          pick the best representative per clade
./run_stage2.sh ./predicted_pdbs \
    -tm 50:80 -topt 35:42 -phopt 7:8.5 \
    -usability 0.35 -solubility 0.4 \
    -clades 8
```

---

## Stage 1 — sequence filtering & clustering

```
Usage: ./run_stage1.sh <input> [options]
```

`<input>` can be a single FASTA file (`.fasta`, `.fa`, `.faa`) or a directory containing one or more such files (which will be concatenated transparently).

### Options

| Flag | Description |
|---|---|
| `-residues <motif>` | Keep sequences matching a regex motif. `.` matches any single residue. E.g. `G.S.G`, `SHD`. Case-insensitive. |
| `-pfam <IDs>` | Comma-separated Pfam accessions (`PFXXXXX`). Keeps sequences with at least one hit at Pfam gathering thresholds via `hmmsearch`. OR-logic across multiple IDs. |
| `-ec <IDs>` | Comma-separated EC numbers; partial specifications allowed (`3.13.-.-` ≡ `3.13`). OR-logic. Uses CLEAN's max-separation inference. |
| `-identity <pct>` | Cluster at `<pct>`% identity using MMseqs2 (`-c 0.8`, cluster-mode 2) and keep one representative per cluster. |

Filters apply in order: **residues → Pfam → EC → identity clustering**. Any subset can be omitted; if no filters are given, sequences pass through unchanged.

### Examples

```bash
# Motif only
./run_stage1.sh seqs.fasta -residues G.S.G

# Pfam + identity clustering
./run_stage1.sh seqs.fasta -pfam PF07519 -identity 90

# Two Pfam IDs (OR), then EC filter
./run_stage1.sh seqs.fasta -pfam PF07519,PF00657 -ec 3.13.1.8

# Everything
./run_stage1.sh seqs.fasta \
    -residues G.S.G -pfam PF07519 -ec 3.13.-.- -identity 95
```

### First-run downloads

- `-pfam` triggers download of the current Pfam-A.hmm (~1.5 GB) into `external/pfam/`.
- `-ec` triggers cloning of [CLEAN](https://github.com/tttianhao/CLEAN) and download of pretrained weights from Google Drive (~few hundred MB) into `external/CLEAN/`. The first actual inference call additionally downloads ESM-1b (~7 GB) into `~/.cache/torch/`.

---

## Between the stages — structure prediction

Stage 2 needs a directory of PDB files corresponding to the Stage 1 non-redundant FASTA. EnzymeSifter is deliberately agnostic to **how** you obtain those structures — fold them with AlphaFold2/3, ESMFold, ColabFold, OmegaFold, or download experimental structures from the PDB, whichever fits your workflow.

Conventions Stage 2 expects:

- One `.pdb` file per protein (multi-chain PDBs are supported).
- PDB filenames are used as the canonical `ID` throughout the rest of the pipeline (no extension). Avoid spaces or special characters; underscores are fine.
- Each PDB should contain SEQRES records — these are what `scripts/extract_sequences.py` reads to recover the amino-acid sequence for the per-chain predictors.

---

## Stage 2 — structural & biophysical screening

```
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]
```

### Filter options

All optional. When supplied, they (a) filter the merged predictions table into a `*_filtered.tsv` and (b) define the scoring rubric used by clade-representative selection.

| Flag | Form | Meaning |
|---|---|---|
| `-solubility <min>` | cutoff | Keep enzymes with `predicted_solubility ≥ min`. |
| `-usability <min>` | cutoff | Keep enzymes with `predicted_usability ≥ min`. |
| `-tm <min>` *or* `<lo:hi>` | cutoff or interval | Either a lower bound on Tm, or an inclusive interval (mid-point preferred when scoring). |
| `-topt <lo:hi>` | interval (required) | Keep enzymes with `Topt` inside `[lo, hi]`. |
| `-phopt <lo:hi>` | interval (required) | Keep enzymes with `pH optimum` inside `[lo, hi]`. |

### Clade options

| Flag | Description |
|---|---|
| `-clades <n>` | Cut the NJ tree into `n` clades (by binary-searching the smallest distance threshold that yields ≤ `n` clades). When combined with filter options, also writes `predictions_output/clade_representatives.tsv` with the best-scoring enzyme per clade. |

### Scoring rubric (when filter flags are present)

Every enzyme that passes all filters receives a score in `[0, 1]` per property:

- **Cutoff properties (`solubility`, `usability`, `tm` in cutoff mode):** normalised across the passing pool — `(value − pool_min) / (pool_max − pool_min)`, so higher is better.
- **Interval properties (`phopt`, `topt`, `tm` in interval mode):** `1 − |value − midpoint| / half_width`, so being close to the midpoint of the interval is best.

The combined score is the arithmetic mean of the per-property scores over the user-specified properties for which that enzyme has a non-NA prediction. Multi-sequence PDBs (heteromers) are excluded from representative selection because per-chain scores cannot be meaningfully averaged into a single structure-level score.

### Examples

```bash
# Bare minimum: just run EnzyMM, predictors, MUSCLE, and the NJ tree
./run_stage2.sh /data/pdbs

# Single-cutoff filter
./run_stage2.sh /data/pdbs -usability 0.35 -solubility 0.4

# Goldilocks-zone interval filtering
./run_stage2.sh /data/pdbs -tm 50:80 -topt 35:42 -phopt 7:8

# Full mode: filters + clade-based representative selection
./run_stage2.sh /data/pdbs \
    -usability 0.35 -tm 50:80 -topt 35:42 -clades 10
```

### First-run downloads

- NetSolP-1.0 from DTU (~5.6 GB) into `external/NetSolP-1.0/`
- pHoptNN cloned from GitHub into `external/pHoptNN/`
- Seq2Topt cloned from GitHub + weight files into `external/Seq2Topt/`

EnzyMM is installed as a pip package via `envs/enzymm.yaml`.

---

## Output structure

```
EnzymeSifter/
├── data/
│   ├── stage1/
│   │   ├── nonredundant.fasta           # ← input to your structure predictor
│   │   ├── clustering_report.tsv        # representative ↔ member mapping
│   │   ├── motif_report.tsv             # per-seq motif hit positions
│   │   ├── pfam_report.tsv              # per-seq Pfam hits
│   │   └── ec_report.tsv                # CLEAN EC predictions
│   ├── enzymm/
│   │   ├── <pdb>.tsv                    # per-PDB EnzyMM output
│   │   └── hit_pdbs.txt                 # PDBs with ≥1 enzymatic hit
│   ├── sequences/
│   │   ├── <pdb>.fasta                  # SEQRES-derived chain sequences
│   │   └── all_hits.fasta               # merged FASTA for downstream tools
│   ├── alignments/muscle.afa            # MUSCLE multiple-sequence alignment
│   ├── trees/
│   │   ├── nj_tree.nwk                  # Newick NJ tree (blosum62 distances)
│   │   ├── nj_tree.png                  # plain tree render
│   │   └── nj_tree_clades.png           # tree coloured by clade + ★ representatives
│   ├── clades/clade_assignments.tsv     # tip_name ↔ clade_id
│   └── predictions/
│       ├── netsolp.tsv                  # raw NetSolP output
│       ├── phoptnn.tsv                  # raw pHoptNN output
│       ├── seq2topt_topt.tsv            # raw Seq2Topt output
│       └── seq2topt_tm.tsv              # raw Seq2Tm output
├── predictions_output/
│   ├── all_predictions.tsv              # one row per PDB (single-chain case)
│   ├── all_predictions_structure.tsv    # structure-level table (multi-chain case)
│   ├── all_predictions_chains.tsv       # per-chain table (multi-chain case)
│   ├── *_filtered.tsv                   # threshold-passing subsets
│   └── clade_representatives.tsv        # one best enzyme per clade
└── logs/                                # per-rule stderr captures
```

### Column reference for the merged tables

| Column | Source | Notes |
|---|---|---|
| `ID` | PDB filename stem | Canonical identifier. |
| `chain_id` | derived from SEQRES | Only present in the multi-chain table. |
| `predicted_solubility` | NetSolP | 0–1, higher = more soluble. |
| `predicted_usability` | NetSolP | 0–1, higher = more usable. |
| `predicted_ph_opt` | pHoptNN | Predicted optimal pH (structure-level). |
| `predicted_topt_C` | Seq2Topt | Predicted optimum temperature, °C. |
| `predicted_tm_C` | Seq2Tm | Predicted melting temperature, °C. |
| `score` | EnzymeSifter | Combined score over user-specified filters, `[0, 1]`. |
| `n_eligible_in_clade` | EnzymeSifter | How many single-sequence PDBs in this clade passed filters. |
| `n_members_in_clade` | EnzymeSifter | Total tips in this clade. |

`NA` indicates that a particular predictor could not score that enzyme (e.g. an excessively long sequence skipped by CLEAN's ESM-1b 1022-aa limit).

---

## External tools

EnzymeSifter orchestrates the following tools. Please cite their authors if you use the corresponding output in published work.

| Tool | Purpose | Reference |
|---|---|---|
| [CLEAN](https://github.com/tttianhao/CLEAN) | EC-number prediction | Yu, Lu, Tianhao *et al.* *Science* 2023 |
| [Pfam / HMMER](http://hmmer.org/) | Domain assignment | Mistry *et al.* *NAR* 2021; Eddy *PLoS Comput. Biol.* 2011 |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Fast clustering | Steinegger & Söding *Nat. Biotechnol.* 2017 |
| [EnzyMM](https://pypi.org/project/enzymm/) | Catalytic-site detection | — |
| [NetSolP-1.0](https://services.healthtech.dtu.dk/services/NetSolP-1.0/) | Solubility/usability | Thumuluri *et al.* *Bioinformatics* 2022 |
| [pHoptNN](https://github.com/kuenzelab/pHoptNN) | Optimal pH | Künzel lab |
| [Seq2Topt](https://github.com/SizheQiu/Seq2Topt) | Topt / Tm | Qiu *et al.* |
| [MUSCLE 5](https://drive5.com/muscle/) | Multiple-sequence alignment | Edgar *Nat. Commun.* 2022 |
| [Biopython](https://biopython.org/) | NJ tree, parsers | Cock *et al.* *Bioinformatics* 2009 |

---

## Disk-space & runtime expectations

| Component | Disk | Notes |
|---|---|---|
| Pfam-A.hmm | ~1.5 GB | One-time download. |
| CLEAN + weights | < 1 GB | Excludes ESM-1b cache. |
| ESM-1b cache | ~7 GB | Stored under `~/.cache/torch/`. Downloaded on first CLEAN inference. |
| NetSolP-1.0 | ~5.6 GB | One-time download. |
| Seq2Topt / pHoptNN | < 1 GB | Both small. |
| `.snakemake/conda/` | ~10 GB | Tool-specific environments. |

Runtime is dominated by (a) HMMER + CLEAN at Stage 1 (scales with the number of input sequences) and (b) NetSolP + Seq2Topt at Stage 2 (scales with the number of unique chains). Both stages parallelise transparently with `snakemake -j`.

---

## Troubleshooting

**Stage 1 with `-ec` fails immediately with an LD/CUDA error.**
The CLEAN environment ships `torch==1.11.0` (CPU). The wrapper exports `LD_LIBRARY_PATH=${CONDA_PREFIX}/lib`. If you see `libstdc++.so.6: GLIBCXX_3.4.X not found`, ensure your conda environment was activated and that no system-level `LD_LIBRARY_PATH` is taking precedence.

**`hmmfetch: index missing` during Pfam filter.**
Re-run `bash scripts/setup_pfam.sh external/pfam` — this rebuilds the `.ssi` index after stripping ACC version suffixes.

**`No SEQRES sequences found in the PDB`.**
Several structure-prediction tools omit SEQRES records by default. Use a tool that emits them (or post-process the PDB to add them) — Stage 2 relies on SEQRES for sequence recovery.

**CLEAN skips long sequences.**
ESM-1b's input limit is 1022 residues. `run_clean_filter.py` skips longer sequences with a `WARN` and records `SKIPPED_TOO_LONG` in the EC report.

**Snakemake complains about a corrupted partial download.**
The `setup_*.sh` scripts retry up to twice on corrupt archives. If they still fail, delete the offending file under `external/` and re-run.

---

## Citing EnzymeSifter

If you use EnzymeSifter in your work, please cite this repository and the underlying tools listed in [External tools](#external-tools).

---

## License

MIT — see `LICENSE`. Please note that the external tools EnzymeSifter wraps have their own licences; consult each project for details.

## Licensing

- **Source code**: MIT License (see `LICENSE`)
- **Predicted structures** (`/pdbs/`): subject to the
  [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms).
  See `/pdbs/TERMS.txt` for details. Non-commercial use only.
  

## Installation

Clone the repository:

```bash
git clone https://github.com/Dar-Omar/EnzymeSifter.git
cd EnzymeSifter
```

Make the run scripts executable:

```bash
chmod +x run_stage1.sh run_stage2.sh
```


