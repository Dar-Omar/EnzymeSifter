# EnzymeSifter tutorial
 
This tutorial walks you through a complete, real EnzymeSifter run from start to finish and then documents every option available in both stages so you can adapt the workflow to your own project. The test run used in the published EnzymeSifter paper is demonstrated here as an example.

---
 
## Table of contents
 
- [Worked example: trypsin-like proteases from a soil metaproteome](#worked-example-trypsin-like-proteases-from-a-soil-metaproteome)
  - [The starting material](#the-starting-material)
  - [Stage 1 — filtering 2.3 M sequences down to 122 representatives](#stage-1--filtering-23-m-sequences-down-to-122-representatives)
  - [Between the stages — structure prediction](#between-the-stages--structure-prediction)
  - [Stage 2 — predicting properties and picking clade champions](#stage-2--predicting-properties-and-picking-clade-champions)
  - [Reading the outputs](#reading-the-outputs)
- [Stage 1 options reference](#stage-1-options-reference)
- [Stage 2 options reference](#stage-2-options-reference)
- [Tips for choosing parameter values](#tips-for-choosing-parameter-values)
---
 
## Worked example: trypsin-like proteases from a soil metaproteome
 
### The starting material
 
The input was a single multi-FASTA of **2,330,712 protein sequences** translated from a soil metaproteomics dataset, with sequence headers normalised to unique IDs:
 
```
~/soil_proteins/soil_proteins_renamed.fasta
```
 
Our goal: produce a short, structurally diverse list of trypsin-like serine proteases that are predicted to be (a) soluble, (b) thermostable, and (c) active under mildly alkaline, ambient-temperature conditions — i.e. plausibly useful enzymes for downstream characterisation.
 
### Stage 1 — filtering 2.3 M sequences down to 122 representatives
 
```bash
./run_stage1.sh ~/soil_proteins/soil_proteins_renamed.fasta \
    -residues GDSGGP \
    -pfam PF00089 \
    -identity 50
```
 
What each flag does, and why these specific values:
 
- **`-residues GDSGGP`** keeps only sequences containing the literal six-residue stretch `G-D-S-G-G-P`. This is the conserved motif around the catalytic serine of chymotrypsin-family serine proteases, so it's a cheap, high-specificity first sieve. The `-residues` flag accepts regular expressions, so you could equally write `G.SGGP` to allow any residue in position 2.
- **`-pfam PF00089`** runs `hmmsearch` (with Pfam gathering thresholds) against the *Trypsin* Pfam HMM. The motif alone catches false positives that happen to contain the hexapeptide; requiring a Pfam Trypsin hit confirms the domain architecture.
- **`-identity 50`** clusters the remaining sequences at 50% identity with MMseqs2 and keeps one representative per cluster. 50% is aggressive — appropriate when you expect substantial sequence diversity in a metagenome and want a tractable number of structures to fold.
The pipeline reports each filtering step in turn:
 
```
Input sequences:         2330712
Motif filter ('GDSGGP'): 367 matched
Pfam filter (PF00089):   218 matched
Non-redundant sequences: 122
Identity threshold:      50.0%
```
 
A few things worth noticing from this funnel:
 
- The motif-to-Pfam ratio (367 → 218) is the expected pattern: the motif alone catches some non-trypsin hits, and the Pfam HMM removes them.
- Reducing 218 sequences to 122 representatives at 50% identity means roughly half of the Pfam hits were near-duplicates of one another — typical for metagenomic data where the same taxon contributes many similar sequences.
- 122 structures is a sensible number to fold next: enough to populate a phylogeny meaningfully, but not so many that ESMFold/AlphaFold will take days.
Reports for inspection: `data/stage1/clustering_report.tsv` (which sequence collapsed into which cluster), `data/stage1/motif_report.tsv` (where the motif was found in each sequence), and `data/stage1/pfam_report.tsv` (which Pfam HMM scored each hit).
 
### Between the stages — structure prediction
 
EnzymeSifter is deliberately agnostic about how you obtain structures. For this example we folded `data/stage1/nonredundant.fasta` with ESMFold on a GPU node and put the resulting 122 PDB files into `~/trypsin_pdbs/`. AlphaFold, ColabFold, OmegaFold, or experimental PDB downloads all work equally well — only the file extension (`.pdb`) and the presence of SEQRES records matter to Stage 2.
 
> The fact that Stage 2 launched **100 jobs** of EnzyMM rather than 122 is just Snakemake batching — each EnzyMM call is one job per PDB; the remaining 22 ran in the next wave shown in the log.
 
### Stage 2 — predicting properties and picking clade champions
 
```bash
./run_stage2.sh ~/trypsin_pdbs/ \
    -solubility 0.6 \
    -tm 55 \
    -phopt 7:9 \
    -topt 30:45 \
    -clades 11
```
 
Why these values:
 
- **`-solubility 0.6`** is a fairly stringent NetSolP cutoff. NetSolP outputs solubility on a 0–1 scale and ≥ 0.6 generally separates "likely to express solubly in *E. coli*" from "probably needs refolding or a fusion partner". Loosen this to 0.4 if your panel is small.
- **`-tm 55`** asks for predicted melting temperatures of at least 55 °C — a soft "is this a mesophile or better?" cut. Because we passed it as a single number (cutoff mode), the score for this property is normalised across the *passing* pool, so higher Tm is rewarded.
- **`-phopt 7:9`** keeps enzymes whose predicted pH optimum sits in `[7, 9]`, with the scoring rubric most rewarding values near the midpoint (pH 8). Trypsin family enzymes are classically active at slightly alkaline pH so this is biologically aligned.
- **`-topt 30:45`** asks for an ambient-to-warm temperature optimum (mid 37 °C). If you wanted a thermophilic variant you'd push this to `55:75`.
- **`-clades 11`** partitions the resulting NJ tree into 11 clades by binary-searching the smallest tip-distance threshold that produces ≤ 11 clades. A rough rule of thumb is `clades ≈ √N` for `N` input structures — with 122 PDBs, 11 sits right at √122 ≈ 11.05 and gives you a structurally diverse shortlist of around a dozen candidates with one representative per clade.
While the pipeline runs you'll see a long burst of `Activating conda environment` lines as Snakemake spins up 122 parallel EnzyMM jobs, then the checkpoint resolves and per-chain predictors (NetSolP, Seq2Topt/Seq2Tm) run on the merged FASTA while pHoptNN runs on the EnzyMM-passing PDBs in parallel. That's the expected pattern — EnzyMM is the bottleneck on first pass; the per-property predictors are fast once their environments are built.
 
### Reading the outputs
 
After the run, `predictions_output/` contains the three files that matter:
 
```
predictions_output/
├── all_predictions.tsv             # one row per PDB, all five predicted properties
├── all_predictions_filtered.tsv    # subset that passed your -solubility/-tm/-phopt/-topt filters
└── clade_representatives.tsv       # one winner per clade among the filtered enzymes
```
 
`clade_representatives.tsv` is the answer to "give me a structurally diverse shortlist of trypsins that satisfy my biophysical criteria". Each row reports the winning PDB for one clade, its combined score in `[0, 1]`, how many filter-passing enzymes were available in that clade (`n_eligible_in_clade`), and how many total members the clade has (`n_members_in_clade`). Clades with no passing members appear with `ID = NA` — useful information, because it tells you that an entire branch of your phylogeny failed the filters.
 
The full table sources are in `data/predictions/` (raw per-tool outputs), and the phylogenetic context — tree, clade assignments, coloured tree PNG with stars marking the clade representatives — lives in `data/trees/` and `data/clades/`.
 
A typical follow-up at this point is to open `data/trees/nj_tree_clades.png`, scan which clades produced winners and which didn't, and use that to decide whether your filters were too strict, whether to repeat Stage 2 with different intervals, or whether the winners are ready for synthesis.
 
---
 
## Stage 1 options reference
 
```
Usage: ./run_stage1.sh <input> [options]
```
 
`<input>` is **required** and can be either a single FASTA (`.fasta` / `.fa` / `.faa`) or a directory containing one or more such files. When given a directory, EnzymeSifter concatenates all FASTAs transparently before processing.
 
All filtering options are independent; supplying none of them runs Stage 1 as a passthrough that simply writes the input back out (useful if you only want to dereplicate, in which case pass `-identity` alone).
 
### `-residues <motif>`
 
Keep only sequences whose amino-acid string matches a regex motif (case-insensitive).
 
- `.` matches any single residue.
- Anything else is interpreted as a literal residue.
- The motif is searched anywhere in the sequence — there's no implicit anchoring.
Examples:
 
```bash
-residues GDSGGP        # literal hexapeptide (classic trypsin signature)
-residues G.S.G         # GxSxG (esterase/lipase nucleophile elbow)
-residues SHD           # catalytic triad in linear order (rare but possible)
-residues H.{1,3}E      # His followed by 1-3 of any residue then Glu
```
 
The matching position is recorded in `data/stage1/motif_report.tsv`, which is handy when you want to confirm the motif sits where you expected (e.g. nucleophile-elbow loops are usually around residue 100–200 in α/β hydrolases).
 
### `-pfam <IDs>`
 
Comma-separated list of Pfam accessions in `PFXXXXX` form. Sequences are kept if they hit **any one** of the listed accessions (OR logic) at Pfam gathering thresholds.
 
```bash
-pfam PF00089                   # Trypsin
-pfam PF07519                   # Tannase / feruloyl esterase
-pfam PF00089,PF00112,PF00326   # Trypsin OR Peptidase_C1 OR Peptidase_S9 — any serine protease family
```
 
How to find the right accession: search [Pfam at InterPro](https://www.ebi.ac.uk/interpro/entry/pfam/) by family name or by a known reference protein. The accession is the `PFXXXXX` code on the family page (drop any version suffix — EnzymeSifter strips `.N` versions itself).
 
On the first run, `setup_pfam.sh` downloads `Pfam-A.hmm` (~1.5 GB) into `external/pfam/` and builds the `hmmfetch` index. Subsequent runs reuse it.
 
### `-ec <IDs>`
 
Comma-separated list of EC numbers. Sequences are kept if [CLEAN's](https://github.com/tttianhao/CLEAN) max-separation predictor assigns them an EC that matches any of the listed specs. Partial specs are wildcards:
 
| Spec | Means |
|---|---|
| `3.13.1.8` | Exactly EC 3.13.1.8 |
| `3.13.-.-` | Any EC starting with `3.13.` |
| `3.13` | Same as `3.13.-.-` |
| `3` | Any EC in class 3 (hydrolases) |
 
Multiple specs use OR logic:
 
```bash
-ec 3.13.1.8                # one specific EC
-ec 3.13.-.-                # all of subclass 3.13
-ec 3.13.1.8,2.5.1.94       # this EC or this one
```
 
Caveats worth knowing:
 
- CLEAN's underlying ESM-1b embedder caps sequences at 1022 residues. Longer sequences are *skipped* with a warning and marked `SKIPPED_TOO_LONG` in the EC report — they do not pass through, they're dropped. Important if you're filtering very long proteins.
- First-run downloads: CLEAN repo + pretrained weights (~few hundred MB from Google Drive), plus ESM-1b on first inference (~7 GB into `~/.cache/torch/`).
### `-identity <pct>`
 
Cluster sequences with MMseqs2 at `<pct>`% identity and keep one representative per cluster. The clustering uses:
 
- `--min-seq-id <pct/100>`
- `-c 0.8` (80% bidirectional coverage)
- `--cov-mode 0`
- `--cluster-mode 2` (greedy in-cluster set cover — fast and good for redundancy reduction)
Typical values:
 
| Value | Use case |
|---|---|
| 95 | Remove near-identical duplicates only (e.g. trim a redundant database) |
| 90 | Conservative dereplication |
| 70 | Standard redundancy reduction for downstream structural analysis |
| 50 | Aggressive — diverse representatives for phylogenetic spread |
| 30 | Family-level representatives only |
 
The full cluster membership ends up in `data/stage1/clustering_report.tsv` with two columns: `representative` and `member`.
 
### Filter order
 
When multiple flags are combined, they always apply in this fixed order:
 
```
-residues  →  -pfam  →  -ec  →  -identity
```
 
The cheap filters run first so the expensive predictors (Pfam HMMER, CLEAN) only see sequences that already passed the simpler tests. This is why the example above goes from 2.3 M → 367 → 218 → 122: motif first (cheapest), Pfam second, identity clustering last.
 
---
 
## Stage 2 options reference
 
```
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]
```
 
`/path/to/pdbs` is **required** and must be a directory containing one or more `.pdb` files. The filename (sans extension) becomes the canonical `ID` everywhere downstream, so make sure filenames are descriptive and stable. Multi-chain PDBs are supported — Stage 2 will produce per-chain plus per-structure output tables automatically.
 
If you pass no other flags, Stage 2 still runs all four property predictors plus EnzyMM, the MUSCLE alignment, and the NJ tree — it just doesn't filter or pick representatives.
 
### Filter options
 
All filter options are optional. Supplying any of them does **two** things:
 
1. Filters `predictions_output/all_predictions.tsv` (or the multi-chain pair) into a `*_filtered.tsv`.
2. Defines the scoring rubric used for clade-representative selection (only meaningful when combined with `-clades`).
#### `-solubility <min>` *(cutoff only)*
 
Keep enzymes whose NetSolP-predicted solubility is ≥ `<min>`. NetSolP outputs a score on `[0, 1]`. Rough interpretive guide:
 
| Value | Interpretation |
|---|---|
| ≥ 0.7 | Strong soluble expression in *E. coli* expected |
| 0.5–0.7 | Probably soluble; minor optimisation may help |
| 0.3–0.5 | Likely needs fusion tag, refolding, or alternate host |
| < 0.3 | Probably insoluble |
 
#### `-usability <min>` *(cutoff only)*
 
Keep enzymes whose NetSolP-predicted *usability* is ≥ `<min>`. Usability combines solubility with secretory-pathway / pro-peptide considerations — useful when picking proteins for higher-throughput biophysical assays.
 
#### `-tm <value>` *(cutoff **or** interval)*
 
The only flag that accepts both forms:
 
- `-tm 55` — cutoff: keep enzymes with `Tm ≥ 55 °C`. Higher = better in the score.
- `-tm 50:80` — interval: keep enzymes with `50 °C ≤ Tm ≤ 80 °C`. Closer to the midpoint (65 °C) = better in the score.
Use the cutoff when "more thermostable is always better" (e.g. industrial biocatalysis). Use the interval when you have an upper bound for some reason (e.g. matching a host's growth temperature, or avoiding cofactor degradation).
 
#### `-topt <lo:hi>` *(interval only)*
 
Predicted optimal operating temperature in °C. Must be an interval — there's no defensible "more is always better" for `Topt` (an enzyme with `Topt = 90 °C` is *bad* if you want activity at room temperature).
 
```bash
-topt 30:45     # mesophilic / ambient
-topt 55:75     # thermophilic
```
 
#### `-phopt <lo:hi>` *(interval only)*
 
Predicted optimal pH. Must be an interval, for the same reason as `-topt`.
 
```bash
-phopt 5:7      # mildly acidic to neutral (e.g. food-grade applications)
-phopt 7:9      # neutral to mildly alkaline (e.g. trypsin-family)
-phopt 9:11     # alkaline (e.g. detergent enzymes)
```
 
### Clade options
 
#### `-clades <n>`
 
Partition the NJ tree into `n` clades by binary-searching the smallest tip-distance threshold that yields ≤ `n` clades, and:
 
- write `data/clades/clade_assignments.tsv` (one row per tip)
- render `data/trees/nj_tree_clades.png` with tips coloured by clade
- if filter options are also present, write `predictions_output/clade_representatives.tsv` containing the best-scoring filter-passing enzyme per clade, and overlay ★ markers on the tree at each representative
The pipeline emits a warning if `n` clades can't be produced exactly — usually because the tree topology forces a slightly different cut.
 
### How scoring works
 
Whenever filter flags are given, every filter-passing enzyme receives a score in `[0, 1]` per specified property:
 
- **Cutoff properties** (`solubility`, `usability`, `tm` in cutoff mode): score is the property's value normalised across the passing pool — `(value − pool_min) / (pool_max − pool_min)`. So `1.0` goes to the best enzyme in the pool, not in some abstract scale.
- **Interval properties** (`phopt`, `topt`, `tm` in interval mode): score is `1 − |value − midpoint| / half_width`, so being at the centre of the requested range = `1.0`, being at the edge = `0.0`.
The combined score is the arithmetic mean of the per-property scores across the user-specified properties, ignoring any property the predictor returned as `NA` for that enzyme. Multi-sequence PDBs (heteromeric assemblies) are excluded from representative selection because their per-chain scores can't be meaningfully averaged into a single structure-level value.
 
---
 
## Tips for choosing parameter values
 
A few rules of thumb that have helped on real datasets:
 
**For `-identity` at Stage 1.** Start aggressive (50–70%) if your input is large and metagenomic; you can always re-run more permissively. Going below 30% risks merging genuinely different sub-families into single clusters.
 
**For `-clades` at Stage 2.** `clades ≈ √N` (where `N` is the number of folded structures) is a sensible default — it scales the number of representatives with the diversity of the input, and tends to produce clades large enough that each contains at least one filter-passing enzyme. If many clades end up with `ID = NA` in the representatives table, your filters are too strict for the tree's diversity; lower `clades` or loosen filters.
 
**For Tm / Topt / pH intervals.** Make them as biologically meaningful as your application allows. A wide interval (`-tm 30:90`) is essentially a passthrough that contributes nothing to scoring; a narrow one (`-tm 65:67`) will rank-order enzymes very tightly around the midpoint. The "sweet spot" is usually a range of about 15–30 units for temperature and 1.5–2.5 units for pH.
 
**On combining filters.** Each filter you add reduces the number of enzymes available for clade-representative selection. With aggressive filters, you may want to *increase* the number of clades — counter-intuitively — to ensure each clade has enough filter-passing members to choose a winner from. Conversely, if you find that almost every enzyme passes, your filters aren't doing much and you'll likely benefit from tightening them.
 
**On reproducibility.** Stage 2's `-clades` cut is deterministic given the same tree and the same `n`, but the tree itself depends on the MUSCLE alignment, which depends on the input order. If you need bit-identical reproducibility, keep your input PDB filenames stable between runs.
