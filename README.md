# EnzymeSifter

Snakemake pipeline integrating EnzyMM, NetSolP, pHoptNN, and Seq2Topt to identify enzymes from PDB structures and select top candidates per clade.

## Licensing

- **Source code**: MIT License (see `LICENSE`)
- **Predicted structures** (`/pdbs/`): subject to the
  [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms).
  See `/pdbs/TERMS.txt` for details. Non-commercial use only.
  
## Requirements
TBA

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

## Usage

### Stage 1 — sequence filtering & clustering

```bash
./run_stage1.sh <input.fasta|input_dir> [-identity <pct>] [-residues <motif>] [-pfam <PFXXXXX>]
```

Examples:
```bash
./run_stage1.sh seqs.fasta -pfam PF00657 -residues G.S.G -identity 90
```

### Stage 2 — enzyme detection, property prediction & clade selection

```bash
./run_stage2.sh /path/to/pdbs [filter options] [-clades <n>]
```

Examples:
```bash
./run_stage2.sh /pdbs -usability 0.35 -tm 50:80
```
