# EnzymeSifter

A two-stage Snakemake pipeline for sifting through sequences and identifying promising enzyme candidates for downstream characterisation — from raw FASTA all the way to ranked, clade-representative enzymes annotated with predicted biochemical values.

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
- [License](#license)

---

## Overview

The pipeline runs in two stages with a structure-prediction step (carried out externally by the user) in between. Users are free to choose any of the optional filters:

1. **Stage 1 — filtering of sequences.** Filter by catalytic-residue motif, one or more Pfam families, and/or CLEAN-predicted EC number(s), then cluster at a user-defined identity threshold using MMseqs2.
2. **(User step) Structure prediction.** Generate PDB structures for the Stage 1 filtered sequences.
3. **Stage 2 — structural & biophysical screening.** Confirm enzymatic activity with EnzyMM, predict solubility/usability (NetSolP), optimal pH (pHoptNN), and optimal/melting temperatures (Seq2Topt/Seq2Tm), build an NJ tree, optionally partition it into clades, and select the best-scoring representative per clade according to your filtering criteria.

---

## Pipeline architecture

```
                          ┌────────────────────────────────┐
                          │           Input FASTA          │
                          └───────────────┬────────────────┘
                                          │
                  ╔═══════════════════════▼═══════════════════════╗
                  ║                    STAGE 1                    ║
                  ║  motif → Pfam (HMMER) → EC (CLEAN) → MMseqs2  ║
                  ╚═══════════════════════╤═══════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │ data/stage1/nonredundant.fasta │
                          └───────────────┬────────────────┘
                                          │
                  ╔═══════════════════════▼═══════════════════════╗
                  ║              predict 3D structures            ║
                  ║                done by the user               ║
                  ╚═══════════════════════╤═══════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │     directory of .pdb files    │
                          └───────────────┬────────────────┘
                                          │
   ╔══════════════════════════════════════▼══════════════════════════════════════╗
   ║                                  STAGE 2                                    ║
   ║   EnzyMM ── filter to enzymatic hits ──┬── NetSolP   (solubility, usable)   ║
   ║                                        ├── pHoptNN   (pH optimum)           ║
   ║                                        ├── Seq2Topt  (T optimum)            ║
   ║                                        ├── Seq2Tm    (T melting)            ║
   ║                                        └── MUSCLE ─► NJ tree ─► clades ──►  ║
   ║                                                       representatives       ║
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

Clone the repository:

```bash
git clone https://github.com/Dar-Omar/EnzymeSifter.git
cd EnzymeSifter
```

Make the run scripts executable:

```bash
chmod +x run_stage1.sh run_stage2.sh
```

That is all you need to do. On the first invocation of either stage, Snakemake will build all of the required tool-specific conda environments under `.snakemake/conda/`, and the relevant setup scripts (`scripts/setup_*.sh`) will fetch external databases and model weights into `external/`.

---

## Stage 1 — sequence filtering & clustering

```
Usage: ./run_stage1.sh <input> [options]
```

`<input>` can be a single FASTA file (`.fasta`, `.fa`, `.faa`) or a directory containing one or more such files (which will be concatenated transparently).

---

## Between the stages — structure prediction

Stage 2 needs a directory of PDB files corresponding to the Stage 1 non-redundant FASTA.

---

## Stage 2 — structural & biophysical screening

```
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]
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

## Troubleshooting

**`No SEQRES sequences found in the PDB`.**
Several structure-prediction tools omit SEQRES records by default. Use a tool that emits them (or post-process the PDB to add them) — Stage 2 relies on SEQRES for sequence recovery.

**CLEAN skips long sequences.**
ESM-1b's input limit is 1022 residues. `run_clean_filter.py` skips longer sequences with a `WARN` and records `SKIPPED_TOO_LONG` in the EC report.

**Snakemake complains about a corrupted partial download.**
The `setup_*.sh` scripts retry up to twice on corrupt archives. If they still fail, delete the offending file under `external/` and re-run.

---


## License

- **Source code**: MIT License (see `LICENSE`)
- **Predicted structures** (`/pdbs/`): subject to the
  [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms).
  See `/pdbs/TERMS.txt` for details. Non-commercial use only.
  


