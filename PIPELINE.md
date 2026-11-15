# Pipeline runtime

Wall-clock time per step for a clean run on an Apple M4 Pro (12 cores, 24 GB RAM). Times were
reconstructed from the start and finish timestamps of each step's successful run during
development. Steps that ran in the background are upper bounds.

| Step | Time | Share |
|---|---|---|
| 00 download (3.2 GB) | ~1.5 min | 3% |
| 01 load + cell calling (22 libraries) | 1.7 min | 3% |
| 01b barcode-curve diagnostics | ~0.3 min | <1% |
| 02 QC + hashtag demux + Scrublet | ~10 min | 17% |
| 02c lineage-doublet filter | 1.1 min | 2% |
| 03 global embedding + clustering | ~10 min | 17% |
| 04 evaluate global clusters | 0.8 min | 1% |
| **05 subclustering of the 8 cell families** | **~20 min** | **34%** |
| **05b refine (choose resolution, merge clusters)** | **~9 min** | **15%** |
| 04 evaluate the 8 cell families | ~1.3 min | 2% |
| 06 finalize + 07 receptor check | ~1 min | 2% |
| **Total** | **~57 min** | |

## Bottlenecks

- **Subclustering (05 + 05b), about half the runtime.** Each of the 8 cell families runs Leiden
  45 times (9 resolutions × 5 seeds) to measure stability. The merge step then re-tests
  neighbouring clusters over ~24k genes until no more merges happen.
- **Global clustering (03).** Harmony takes about 1 s and PCA under 30 s. The time goes to the
  Leiden stability sweep (30 runs on 367k cells), the neighbour graph, UMAP and writing the output.
- **QC (02).** Scrublet run separately on each of 22 libraries, plus writing a 1.4 GB merged file.
- **Mostly single-threaded.** CPU sat ~92% idle (one of 12 cores busy) and peak memory was ~10 GB.
  The Leiden seeds and resolutions, the 8 families and the 22 QC libraries are independent, so
  running them in parallel is the main available speed-up.
