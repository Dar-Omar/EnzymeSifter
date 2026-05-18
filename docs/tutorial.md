# EnzymeSifter tutorial
 
This tutorial walks you through a complete, real EnzymeSifter run from start to finish and then documents every option available in both stages so you can adapt the workflow to your own project. The test run used in the published EnzymeSifter paper is demonstrated here as an example.

---
 
## Table of contents
 
- [Worked example: trypsins from soil samples](#worked-example-trypsins-from-soil-samples)
  - [The starting material](#the-starting-material)
  - [Stage 1 — filtering sequences](#stage-1--filtering-sequences)
  - [Between the stages — structure prediction](#between-the-stages--structure-prediction)
  - [Stage 2 — predicting properties and picking clade champions](#stage-2--predicting-properties-and-picking-clade-champions)
  - [Reading the outputs](#reading-the-outputs)
- [Stage 1 options reference](#stage-1-options-reference)
- [Stage 2 options reference](#stage-2-options-reference)
- [Tips for choosing parameter values](#tips-for-choosing-parameter-values)
---
 
## Worked example: trypsins from soil samples
 
### The starting material
 
The input was a single multi-FASTA of **2,330,712 protein sequences** available at. Sequence headers were cut at the first (–) to avoid very long headers, the gene number was kept so no duplicates appear in the initial dataset. The command used: 
 
```
sed -E 's/^(>[^-]+)-NODE-[^_]+_([0-9]+).*/\1_\2/' soil_proteins_combined.fasta > soil_proteins.fasta
```
 
### Stage 1 — filtering sequences
 
```bash
./run_stage1.sh ~/soil_proteins_renamed.fasta -residues GDSGGP -pfam PF00089 -identity 50
```
 
- **`-residues GDSGGP`** keeps only sequences containing the literal six-residue motif `G-D-S-G-G-P`, A highly conserved motif in trypsins. This is a fast, high-specificity first sieve. The `-residues` flag accepts regular expressions, so you could write `G.SGGP` to allow any residue in position 2.
- **`-pfam PF00089`** runs `hmmsearch` against Pfam HMM. The motif alone catches false positives that happen to contain the hexapeptide; requiring a Pfam Trypsin hit confirms the domain architecture.
- **`-identity 50`** clusters the remaining sequences at 50% identity with MMseqs2 and keeps one representative per cluster.

Stage 1 acted as a funnel and reduced the number of sequences from > 2.3 million to 122 trypsins using the above identified filters.

### Between the stages — structure prediction
 
Only the file extension (`.pdb`) and the presence of SEQRES records matter to Stage 2. After submitting our 122 sequences to AlphaFold server, the jobs were named as the headers but with a - instead of the dot, gene numbers were also removed from the job name as no duplicates remained in the filtered dataset. PDBs are available at 

### Stage 2 — structural features
 
```bash
./run_stage2.sh ~/trypsin_pdbs/ -solubility 0.6 -tm 55 -phopt 7:9 -topt 30:45 -clades 11
```
 
- **`-solubility 0.6`** NetSolP outputs solubility on a 0–1 scale.
- **`-tm 55`** asks for predicted melting temperatures of at least 55 °C. Because we passed it as a single number (cutoff mode), the score for this property is normalised across the *passing* pool, so higher Tm is rewarded.
- **`-phopt 7:9`** keeps enzymes whose predicted pH optimum sits in `[7, 9]`, with the scoring rubric most rewarding values near the midpoint (pH 8).
- **`-topt 30:45`** Same logic as in pH optimum but for optimum temperature.
- **`-clades 11`** partitions the resulting NJ tree into 11 clades.

 
### Reading the outputs
 
After the run, `predictions_output/` contains three files:
 
```
predictions_output/
├── all_predictions.tsv             # one row per PDB, all five predicted properties
├── all_predictions_filtered.tsv    # subset that passed your -solubility/-tm/-phopt/-topt filters
└── clade_representatives.tsv       # one winner per clade among the filtered enzymes
```
 
`clade_representatives.tsv` is the answer to "give me a structurally diverse shortlist of trypsins that satisfy my biophysical criteria". Each row reports the winning PDB for one clade, its combined score in `[0, 1]`, how many filter-passing enzymes were available in that clade (`n_eligible_in_clade`), how many total members the clade has (`n_members_in_clade`), and the predicted values for that enzyme. Clades with no passing members appear with `ID = NA`, so it tells you that no member of this clade passed the filters.
 
The coloured tree PNG with stars marking the clade representatives lives in `data/trees/`. A typical follow-up at this point is to open `data/trees/nj_tree_clades.png`, scan which clades produced winners and which didn't, and use that to decide whether your filters were too strict, whether to repeat Stage 2 with different filters, or whether the winners are ready for synthesis.
 
---
 
## Stage 1 options reference
 
```
Usage: ./run_stage1.sh <input> [options]
```
 
`<input>` is **required** and can be either a single FASTA (`.fasta` / `.fa` / `.faa`) or a directory containing one or more such files. When given a directory, EnzymeSifter concatenates all FASTAs transparently before processing.
 
All filtering options are independent; supplying none of them runs Stage 1 as a passthrough that simply writes the input back out.
 
### `-residues <motif>`
 
Keep only sequences whose amino-acid string matches a regex motif (case-insensitive).
 
- `.` matches any single residue.
- Anything else is interpreted as a literal residue.
- The motif is searched anywhere in the sequence.
Examples:
 
```bash
-residues GDSGGP        # literal hexapeptide
-residues G.S.G         # GxSxG
```
 
The matching position is recorded in `data/stage1/motif_report.tsv`, which is handy when you want to confirm the motif sits where you expected (e.g. nucleophile-elbow loops are usually around residue 100–200 in α/β hydrolases).
 
### `-pfam <IDs>`
 
Comma-separated list of Pfam accessions in `PFXXXXX` form. Sequences are kept if they hit **any one** of the listed accessions (OR logic) at Pfam gathering thresholds.
 
```bash
-pfam PF00089                   # Trypsin
-pfam PF00089,PF07519           # Trypsin OR Tannase
```
 
On the first run, `setup_pfam.sh` downloads `Pfam-A.hmm` into `external/pfam/` and builds the `hmmfetch` index. Subsequent runs reuse it.
 
### `-ec <IDs>`
 
Comma-separated list of EC numbers. Partial specs are wildcards:
 
| Spec | Means |
|---|---|
| `3.13.1.8` | Exactly EC 3.13.1.8 |
| `3.13.-.-` | Any EC starting with `3.13.` |
| `3.13` | Same as `3.13.-.-` |
| `3` | Any EC in class 3 (Hydrolases) |
 
Multiple specs use OR logic:
 
```bash
-ec 3.13.-.-,2.5.1.94
```
 
- CLEAN's underlying ESM-1b embedder caps sequences at 1022 residues. Longer sequences are *skipped* with a warning and marked `SKIPPED_TOO_LONG` in the EC report.
- First-run downloads: CLEAN repo + pretrained weights, and ESM-1b on first inference.

### `-identity <pct>`
 
Cluster sequences with MMseqs2 at `<pct>`% identity and keep one representative per cluster. The full cluster membership ends up in `data/stage1/clustering_report.tsv` with two columns: `representative` and `member`.
 
### Filter order
 
When multiple flags are combined, they always apply in this fixed order:
 
```
-residues  →  -pfam  →  -ec  →  -identity
```
 
The cheap filters run first so the expensive predictors only see sequences that already passed the simpler tests.
 
---
 
## Stage 2 options reference
 
```
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]
```
 
`/path/to/pdbs` is **required** and must be a directory containing the `.pdb` files. The filename becomes the canonical `ID` downstream, so make sure filenames are descriptive and stable. Multi-chain PDBs are supported — Stage 2 will produce per-chain plus per-structure output tables automatically.
 
If you pass no other flags, Stage 2 still runs all four property predictors plus EnzyMM, the MUSCLE alignment, and the NJ tree — it just doesn't filter or pick representatives.
 
### Filter options
 
All filters are optional. Supplying any of them does **two** things:
 
1. Filters `predictions_output/all_predictions.tsv` (or the multi-chain pair) into a `*_filtered.tsv`.
2. Defines the scoring rubric used for clade-representative selection (best when combined with `-clades`).
#### `-solubility <min>` *(cutoff only)*
 
Keep enzymes whose NetSolP-predicted solubility is ≥ `<min>`.
 
#### `-usability <min>` *(cutoff only)*
 
Keep enzymes whose NetSolP-predicted *usability* is ≥ `<min>`.
 
#### `-tm <value>` *(cutoff **or** interval)*
 
The only flag that accepts both forms:
 
- `-tm 40` — cutoff: keep enzymes with `Tm ≥ 55 °C`. Higher = better in the score.
- `-tm 50:80` — interval: keep enzymes with `50 °C ≤ Tm ≤ 80 °C`. Closer to the midpoint (65 °C) = better in the score.
Use the cutoff when "more thermostable is always better". Use the interval when you have an upper bound.
 
#### `-topt <lo:hi>` *(interval only)*
 
Predicted optimal operating temperature in °C. Must be an interval.
 
#### `-phopt <lo:hi>` *(interval only)*
 
Predicted optimal pH. Must be an interval.
 
### Clade options
 
#### `-clades <n>`
 
Partition the NJ tree into `n` clades by binary-searching the smallest tip-distance threshold that yields ≤ `n` clades, and:
 
- write `data/clades/clade_assignments.tsv` (one row per tip)
- render `data/trees/nj_tree_clades.png` with tips coloured by clade
- if filter options are also present, write `predictions_output/clade_representatives.tsv` containing the best-scoring filter-passing enzyme per clade, and overlay ★ markers on the tree at each representative
The pipeline emits a warning if `n` clades can't be produced exactly — usually because the tree topology forces a slightly different cut.
 
### How scoring works
 
Whenever filter flags are given, every filter-passing enzyme receives a score in `[0, 1]` per specified property:
 
- **Cutoff properties** (`solubility`, `usability`, `tm` in cutoff mode): score is the property's value normalised across the passing pool — `(value − pool_min) / (pool_max − pool_min)`. So `1.0` goes to the best enzyme in the pool.
- **Interval properties** (`phopt`, `topt`, `tm` in interval mode): score is `1 − |value − midpoint| / half_width`, so being at the centre of the requested range = `1.0`, being at the edge = `0.0`.
The combined score is the arithmetic mean of the per-property scores across the user-specified properties, ignoring any property the predictor returned as `NA` for that enzyme. Multi-sequence PDBs (heteromeric assemblies) are excluded from representative selection because their per-chain scores can't be meaningfully averaged into a single structure-level value.
 
---
