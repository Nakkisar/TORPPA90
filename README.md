# TORPPA90

**T**opology-**O**riented **R**esidue **P**reference **P**robability **A**rray — an **asymmetric, structure-aware amino acid substitution matrix** for BLASTP, derived from neural-network structure-to-sequence predictions rather than from evolutionary sequence alignments.

Built to test whether substitution scores grounded in *local structural environment* can recover remote homologs that BLOSUM misses in the twilight zone below 30% sequence identity.

> This repository is the **deliverable** — the matrix itself, ready to install into BLAST+.
> For the full research record, benchmarking data and the alternative TM-Vec/GASS pipeline, see **[TORPPA-matrix](https://github.com/Nakkisar/TORPPA-matrix)**.

---

## Why an asymmetric matrix

Standard substitution matrices are symmetric: BLOSUM62 scores A→B and B→A identically, because the pair's score is derived by averaging observations in both directions. That is a convenience of construction, not a biological fact. Whether a substitution preserves the local three-dimensional environment around a residue depends on the *direction* of the swap — replacing a small residue with a bulky one in a packed core is a different event from the reverse.

TORPPA scores are computed from a directional conditional frequency, so forward and reverse differ. For example:

| Substitution | TORPPA90 score |
|---|---|
| Isoleucine → Glycine | −7.77 |
| Glycine → Isoleucine | −9.15 |

The "90" refers to the **90% maximum mutual sequence identity** of the training set, deliberately matched to BLOSUM90 so that comparisons between them minimise the confound of evolutionary distance. A second variant, **TORPPA60**, was derived at a 60% threshold on the same reasoning that produced the BLOSUM series.

## How it was built

A new [SeqPredNN](https://github.com/falategan/SeqPredNN) model was trained for this work — a neural network that predicts a residue identity from its local structural environment.

| Step | Detail |
|---|---|
| Source | RCSB PDB, culled with [PISCES](https://dunbrack.fccc.edu/pisces/) |
| Filters | X-ray only, resolution ≤ 2.2 Å, R-factor ≤ 0.25, chain length 40–10,000, max mutual identity 90% |
| Culled set | 32,558 entries (19 March 2024); 30,388 structure files retrieved |
| Training | 200 epochs, 20% test split, oversampling to correct amino-acid class imbalance |
| Matrix | 20×20 log-odds over the model's **unnormalised confusion matrix** (true residue → predicted residue), base-2, scaled by 2, following Henikoff & Henikoff (1992) |

Because the observed frequency term is conditioned on the true residue while the background term is conditioned on the predicted residue, the resulting matrix is asymmetric by construction. The derivation generalises: the same equation applied to any SeqPredNN confusion matrix yields a corresponding TORPPA matrix.

## How it compares to BLOSUM

Correlation against the evolutionary matrices — moderately strong, but not redundant:

| Comparison | Linear correlation | Relative entropy |
|---|---|---|
| TORPPA90 vs BLOSUM62 | 0.72 | 0.15 |
| TORPPA90 vs BLOSUM90 | 0.71 | 0.16 |

Column-wise Spearman analysis shows markedly lower agreement above the diagonal than below it, and negative correlation in specific residue columns — directional structural dependencies that are mathematically absent from any symmetric sequence-derived matrix.

### Benchmark: all-against-all structural alignment

562 diverse protein chains sampled across the seven true SCOPe classes, giving **~157,641 unique pairs**. Each pair was aligned twice, since TM-score is itself asymmetric, producing over 150,000 scored alignments. Structural ground truth from TM-align within [US-align](https://zhanggroup.org/US-align/). 54 hours on 12 CPUs.

Scores were length-normalised and min-max scaled to 0–1 for direct comparison:

| Matrix | Mean | Std dev | Min | Max | Range |
|---|---|---|---|---|---|
| **TORPPA90** | **0.630** | 0.201 | 0.185 | 0.992 | 0.807 |
| **TORPPA60** | **0.617** | 0.212 | 0.124 | 0.983 | 0.858 |
| BLOSUM90 | 0.473 | 0.189 | 0.047 | 0.994 | 0.947 |
| BLOSUM62 | 0.465 | 0.182 | 0.031 | 0.949 | 0.918 |

**The result, stated plainly:** TORPPA maintains a substantially higher scoring floor than BLOSUM in the low-similarity range — roughly 0.2 relative units higher across TM-scores of 0.0–0.5, and 0.2–0.3 units for TORPPA60 against BLOSUM62 between TM 0.2 and 0.5. Above TM 0.8 the two families converge and TORPPA offers no advantage.

### Benchmark: MALISAM

[MALISAM](http://prodata.swmed.edu/malisam/) is a hand-curated set of **non-homologous** protein pairs that share near-identical folds at essentially zero sequence identity — the hardest available test for structural recognition where evolutionary signal is deliberately useless.

Here the pattern inverts. In the all-against-all test TORPPA and BLOSUM converge as TM-score rises; on MALISAM pairs the gap between them **widens** with increasing TM-score. BLOSUM's position-independent substitution probabilities bottom out, while TORPPA continues to track the structural signal that remains.

This is the clearest evidence in the project that substitution matrices can be built without evolutionary sequence history and still recognise genuine structural relationships.

### Benchmark: BLASTP retrieval

Here the matrix did not perform. Using haemoglobin beta (6KYE chain D) as query against PDBAA:

| Matrix | Top hit | Bits | Hits returned |
|---|---|---|---|
| TORPPA90 | 1DXT_B | 278 | **238** |
| BLOSUM62 | 1DXT_B | 273 | 427 |
| BLOSUM90 | 1DXT_B | 312 | 330 |

TORPPA90 returned roughly 55% of BLOSUM62's hits — **and missed the target structural homolog, myoglobin 3RGK_A, entirely**. Forced into a direct one-to-one alignment, TORPPA90 did produce the highest proportion of "positives" (40%, against 39% and 34%) at a lower bit score, consistent with prioritising structural over evolutionary similarity. But BLASTP's heuristics — word matching, E-value estimation — are built around symmetric matrices, and an asymmetric matrix is hamstrung inside that framework.

## Honest verdict

**What was demonstrated.** A substitution matrix can be built entirely from machine-learned structural predictions, with no evolutionary sequence history, and still produce a scoring pattern closely aligned with the gold-standard BLOSUM series. In the twilight zone, where sequence identity is actively unhelpful, it tracks structural similarity better than BLOSUM does.

**What stopped it being deployed.** The elevated scoring floor that makes TORPPA sensitive is the same property that compresses its dynamic range: 0.807 against BLOSUM90's 0.947. Less contrast between a random match and a true structural hit means more baseline noise, and in a large automated search that means a high false-positive rate. The matrices were judged **not yet viable for production-level homolog identification**, and the project pivoted to an alternative retrieval strategy for the actual discovery task.

A domain-specific variant — a SeqPredNN model trained only on plastic-degrading enzymes — was scoped and abandoned: PlasticDB held roughly 300 characterised proteins at the time, too sparse to overcome SeqPredNN's recall and accuracy limits, and the likely outcome was a high-noise matrix that could not reliably separate true homologs.

**What would fix it.** Calibration aimed specifically at lowering the scoring baseline. The sensitivity is real; the problem is contrast, not signal.

---

## Installation

BLAST+ compiles substitution matrices into the binary, so TORPPA90 must be built in from source rather than supplied at runtime.

**1. Get the BLAST+ source.**

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-src.tar.gz
tar -xzf ncbi-blast-2.15.0+-src.tar.gz
```

This project used 2.14.1, archived [here](https://stellenbosch-my.sharepoint.com/:f:/g/personal/22042636_sun_ac_za/ElZKYJQcacBGib6HmJdSe4UBgASEkgZaO8vmaFO869cJoA?e=VOptg4). Later versions should work but are untested.

**2. Replace the BLOSUM90 table** with the pre-transformed TORPPA90 table in `BLAST+/` of this repository:

```bash
cp BLAST+/sm_blosum90.c  ncbi-blast-2.15.0+-src/c++/src/util/tables/
```

TORPPA90 occupies BLOSUM90's slot because BLAST+ offers no mechanism for registering a new named matrix. **After this step `-matrix BLOSUM90` invokes TORPPA90**, so keep a separate unmodified build if you need genuine BLOSUM90 for comparison.

**3. Build.**

```bash
cd ncbi-blast-2.15.0+-src/c++
./configure && make
```

**4. Search.**

```bash
makeblastdb -in pdbaa -dbtype prot -out pdbaa
blastp -query my_protein.fasta -db pdbaa -matrix BLOSUM90 -out hits.txt
```

Searches in this project used gap opening 10 and gap extension 1, against the preformatted [pdbaa](https://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz) database.

To build the table yourself from the raw matrix rather than using the supplied one, the values in `TORPPA90_matrix/` must be **transposed** first — BLAST+ reads its tables in the opposite orientation to the matrix as published here.

> Given the retrieval limitation above, TORPPA90 is better suited to **rescoring a fixed alignment** than to driving the search itself. In the benchmarks it was strongest where the alignment was already determined structurally and the matrix only had to score it.

---

## Repository layout

| Path | Contents |
|---|---|
| `TORPPA90_matrix/` | The 20×20 matrix, as published |
| `BLAST+/` | Pre-transformed `sm_blosum90.c` table, ready to drop into the BLAST+ source tree |
| `Matrices/` | Intermediate and alternative variants generated during development, including TORPPA60 |
| `SeqPredNN/` | SeqPredNN model outputs and the confusion matrix the log-odds are derived from |
| `Correlations/` | Spearman, linear correlation and relative entropy analyses against BLOSUM62 and BLOSUM90 |
| `6KYE_D native vs predicted/` | Positive-control comparison of the native 6KYE chain D sequence against SeqPredNN's prediction for the same structure |

## Built with

[SeqPredNN](https://github.com/falategan/SeqPredNN) (Lategan et al., 2023) — structure-conditioned residue prediction, the source of the substitution statistics. [PISCES](https://dunbrack.fccc.edu/pisces/) for dataset culling, [US-align](https://zhanggroup.org/US-align/) for structural ground truth.

## Citation

> Torppa, J. (2026). *Development of a Pipeline to Identify Putative Bioplastic Hydrolases*. MSc thesis, Centre for Bioinformatics and Computational Biology, Stellenbosch University. Chapter 4: A Novel Structure-Aware Substitution Matrix.

## Contact

Johannes Torppa — jetorppa@gmail.com
MSc Computational Biology, Centre for Bioinformatics and Computational Biology, Stellenbosch University
