# EnzymeSifter

A two-stage Snakemake pipeline for sifting through protein sequences and identifying promising enzyme candidates for downstream characterisation — from raw FASTA all the way to ranked, clade-representative enzymes annotated with predicted biochemical values.

---

## Table of contents

- [Overview](#overview)
- [Pipeline architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#Usage)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Overview

The pipeline runs in two stages with a structure-prediction step (carried out externally by the user). Users are free to choose any of the optional filters:

1. **Stage 1 — filtering of sequences.** Filter by catalytic-residue motif, one or more Pfam families, and/or CLEAN-predicted EC number(s), then cluster at a user-defined identity threshold using MMseqs2.
2. **(User step) Structure prediction.** Generate PDB structures for Stage 1 filtered sequences.
3. **Stage 2 — structural screening.** Confirm enzymatic activity with EnzyMM, predict solubility/usability (NetSolP), optimal pH (pHoptNN), and optimal/melting temperatures (Seq2Topt), build an NJ tree, optionally partition it into clades, and select the best-scoring representative per clade according to your filtering criteria.

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

## Usage


### Stage 1

```
./run_stage1.sh <input> [options]
```
#### Example
```bash
./run_stage1.sh ~/soil_proteins_renamed.fasta -residues GDSGGP -pfam PF00089 -identity 50
```

### Between the stages — structure prediction

Stage 2 needs a directory of PDB files corresponding to the Stage 1 non-redundant FASTA.

### Stage 2

```
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]
```
#### Example
```bash
./run_stage2.sh ~/trypsin_pdbs/ -solubility 0.6 -tm 55 -phopt 7:9 -topt 30:45 -clades 11
```

See tutorial for complete features of the pipeline

---

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
  


