# Detailed methods and decisions (GSE337533 clustering)

The big-picture summary is in README.md. This file keeps every step, threshold and judgment call.

Clustering of all cells in the 22 gene-expression 10x libraries of GEO GSE337533. Published
cell-type annotations and the associated publication/preprint were **not** used; annotations
below are my own, from marker genes plus an orthogonal BCR check.

## Main outputs

| file | content |
|---|---|
| `results/cell_clusters.csv` | **one row per retained cell (367,389)**: `cell_id` (`<library>_<barcode>`), `cluster` (final label, e.g. `B_lineage_5`), `annotation`, `compartment`, `global_cluster`, library/barcode/enrichment/donor, `label_source`, `transfer_confidence`, `qc_flag` |
| `results/removed_cells.csv` | 64,642 called cells removed before clustering, with reason |
| `results/cluster_annotations.csv` | 65 clusters → putative identity + key markers |
| `results/vdj_validation_by_cluster.csv` | BCR chain status per cluster (orthogonal validation) |
| `results/figures/`, `results/global/`, `results/sub/<compartment>/` | UMAPs, marker dotplots, per-cluster QC/composition/DE tables |

`cell_clusters.csv` has 432,031 − 64,642 = 367,389 rows. Every cell called in step 1 is in exactly
one of the two CSVs.

## Pipeline (`scripts/`, run with `pixi run python -I scripts/<step>.py …`)

1. **01_load.py: per-library loading + cell calling.** Of the 22 GEO matrices, 14 are Cell Ranger
   *filtered* and 8 are *raw* (~730k barcodes). For raw ones I called cells at the inflection of the
   log barcode-rank curve (GEX UMIs, floor 500). Curves had clean cliffs except BM_2_220524
   (degraded profile; kept the knee, QC handles the rest). Barcode overlap between libraries is at
   chance level, so the 22 libraries are independent 10x reactions, not re-sequencing.
2. **02_qc.py: QC, hashtag demux, doublets.** The "Antibody Capture" features include ~137 ADTs,
   but they are **all zero** (no ADT library deposited), so clustering is RNA-only. Hashtags: run
   "220524" (6 libraries) uses Cellhash1-3, the other 16 libraries use Cellhash4-6.
   * QC per library, MAD-based (depth differs 5×): log UMIs/genes within −3/+4 MAD (floors 500
     UMIs / 250 genes), mito % ≤ median + 4 MAD (clamped to 10-25%).
   * Hash demux: per hashtag, threshold = median + 5 MAD of log counts among cells dominated by
     another hashtag. My first try was a 2-component GMM, which failed when one donor was ~80% of a
     library. Hash doublets were removed; hash-negatives were kept.
   * Scrublet per library with a fixed score threshold of 0.25. The automatic threshold was unstable
     (0-1,000 calls per library). Calibration against hash doublets: 28% of them have score > 0.25
     vs 3.5% of hash singlets (most cross-donor doublets are same-cell-type and invisible to
     Scrublet) (`results/qc/doublet_calibration.txt`).
3. **02c_lineage_doublets.py: heterotypic doublet filter.** The first global clustering
   (`results/global_v1`) showed a B+granulocyte co-expressing population coming almost entirely
   from library Bcell_2_2. I flag cells where two mutually exclusive lineage programs
   (B / granulocytic-myeloid / erythroid / T-NK, each the mean of 6-8 genes) both exceed 1.0
   log-norm. This removed 5,317 cells: 9.7% of Bcell_2_2 and 0.1-2% of the other libraries.
   IGLL1/CD24 are excluded from the B score because myeloid progenitors and neutrophils express them.
4. **03_global.py: global embedding.** CP10k + log1p. 3,000 HVGs (seurat_v3 on counts,
   batch-aware), **excluding Ig/TCR V(D)J and constant genes** (clonotype, light-chain isotype and
   ambient plasma-cell Ig would otherwise drive clusters; IGLL1 is kept), mito, ribosomal, sex-linked
   genes (donor identity) and MALAT1/NEAT1. Cell-cycle genes are **kept**, because proliferation is
   biology here (pre-B bursts, progenitors). PCA 50, then Harmony on `library` using 30 PCs. Library
   is the finest technical unit; donor and run are nested in it. Enrichment (BM/CD34/B) changes
   composition but is not used as a batch. Then kNN (k=20) and Leiden at several resolutions, each
   with 5 seeds. Res 0.2 gave 12 clusters with seed ARI 0.96 and became the compartment layer.
   * kNN uses a small pynndescent wrapper (`knn.py`) because scanpy's built-in pynndescent path hits
     a numba typing error with the installed versions. It was validated at 98% recall against exact
     kNN.
5. **Compromised library.** Even after step 3, Bcell_2_2 cells had pervasive cross-lineage signal:
   30% of its B-cluster cells had myeloid score > 0.5 (other 20 libraries ≤ 2%), and 50% of its
   myeloid-cluster cells had B score > 0.5 (others ≤ 5%). It also had low UMIs and failed hashtag
   calling, and only its filtered matrix is on GEO, so ambient RNA can't be corrected.
   **Decision:** Bcell_2_2 is not used to *define* sub-clusters. Its 31k cells are labelled
   afterwards by kNN majority vote (k=30) in the global Harmony space (`label_source=knn_transfer`,
   `transfer_confidence`; median 0.87) and flagged in `qc_flag`. The same donors (4-6) are covered
   by the five other main-run B libraries.
6. **05_subcluster.py + 05b_refine.py: compartment sub-clustering.** Compartments are global
   clusters grouped by lineage (`results/compartments.json`). Each is re-processed (own HVGs, PCA 30,
   Harmony, kNN). Resolution and merging rules:
   * Resolution = the **most seed-stable** resolution (max mean ARI over 5 Leiden seeds, res ≥ 0.2,
     requiring ARI ≥ 0.75; else res 0.1). My first rule, "highest resolution with ARI ≥ 0.75",
     over-split small compartments (Plasma 19, Erythroid 24).
   * Clusters smaller than max(30 cells, 1% of compartment) merge into the nearest cluster.
   * Iterative merge of nearest pairs with < 10 DE genes (|log2FC| > 1, padj < 0.01, ≥ 20 pp
     detection difference), ignoring Ig/TCR, mito/ribo, sex genes and immediate-early/dissociation
     genes (FOS/JUN/ATF3/NR4A/HSPA1…). This way clusters that differ only by clonotype, isotype or
     handling stress get merged.
   * Result: 65 clusters (B 10, T/NK 7, HSPC/myeloid progenitors 9, monocyte 4, neutrophil 6,
     erythroid 8, plasma 11, pDC 10).
7. **06_finalize.py**: assembles the CSVs and figures. **07_vdj_validation.py**: BCR check.

## Validation

* **Batch:** every major cluster draws cells from 13-21 libraries, and no single library supplies
  more than 18% of any non-doublet cluster with >1,000 cells. The one exception is
  Neutrophil_0 (33%): neutrophils are almost only present in the 6 whole-BM libraries, because the
  CD34 and B-cell enrichments remove them. (`max_lib_frac` in `results/sub/*/sub_summary.csv`).
* **BCR (orthogonal):** Ig genes were excluded from the clustering features, so V(D)J status is
  independent evidence. It reproduces the developmental order:

  | cluster | none | IGH only | IGH+L | class-switched |
  |---|---|---|---|---|
  | B_lineage_7 / 6 (pro-B) | 94% | 1-2% | 3% | ~1% |
  | B_lineage_5 (cycling large pre-B) | 28% | **59%** | 13% | 0% |
  | B_lineage_8 (late large pre-B) | 40% | **47%** | 12% | 0% |
  | B_lineage_2 (small pre-B / immature) | 36% | 33% | 28% | 0% |
  | B_lineage_0 (naive) | 32% | 5% | **57%** | 1% |
  | B_lineage_1 / 9 (memory-like / IEG-high) | ~22% | 2% | ~60% | **14% / 32%** |
  | Plasma_0/1/5 | 2% | 1% | **95%** | **82-87%** |
  | T, monocytes, HSPC | 98-99% | 0% | ~1% | n/a |

  The 30-40% "none" in B cells reflects V(D)J capture sensitivity; non-B cells give the
  false-positive floor (~1%).

## Caveats

* 6,712 retained cells sit in 17 small clusters annotated *Doublet-like* or *ambient*. They are kept
  and labelled rather than dropped, so users can filter on `annotation`.
* Small doublet groups inside the large B compartment (≈700 B+myeloid, ≈600 B+T) and a 266-cell
  G2/M pro-B group were absorbed by the 1% size floor, mostly into B_lineage_1. Residual homotypic
  doublets (especially in Bcell_2_2) can't be detected with the available data.
* Activated/IEG-high states (B_lineage_3/4/9, T_NK_2) may partly reflect ex vivo handling. They
  survive the merge test on non-IEG genes (e.g. CD83, HLA-II, class-switch rate) and are labelled
  accordingly.
* Annotations are marker-based interpretation; cluster IDs are the primary result.

## Reproduce
```
pixi install
pixi run python -I scripts/01_load.py data/raw data/libs
pixi run python -I scripts/01b_qc_loading.py data/raw data/libs results/qc
pixi run python -I scripts/02_qc.py data/libs data/qc
pixi run python -I scripts/02c_lineage_doublets.py data/qc/merged_qc.h5ad data/qc/merged_qc_ld.h5ad results/qc
pixi run python -I scripts/03_global.py data/qc/merged_qc_ld.h5ad data/global 30
pixi run python -I scripts/04_evaluate.py data/global/global.h5ad leiden_0.2 results/global
pixi run python -I scripts/05_subcluster.py data/global/global.h5ad data/qc/merged_qc_ld.h5ad results/compartments.json data/sub
pixi run python -I scripts/05b_refine.py data/sub data/sub/subcluster_resolution_sweep.csv
pixi run python -I scripts/06_finalize.py data/global/global.h5ad data/sub/subcluster_labels.csv results/cluster_annotations.csv data/qc results/qc/lineage_doublet_calls.csv.gz results
pixi run python -I scripts/07_vdj_validation.py data/vdj results/cell_clusters.csv results/vdj_validation_by_cluster.csv
```
Raw data: GEO GSE337533 supplementary `*_{barcodes,features,matrix}` files → `data/raw/`; BCR
`*_filtered_contig_annotations.csv.gz` (GSM9857521-542) → `data/vdj/`. Runtime is about 2-3 h on
an M4 Pro with 24 GB, with peak RAM around 10 GB.
