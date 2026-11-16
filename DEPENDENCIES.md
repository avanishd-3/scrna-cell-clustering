# Dependencies

Versions are the ones installed in the pixi environment that produced the results. The pinned set is in `pixi.toml` and `pixi.lock`.

## Declared in `pixi.toml`

Every package the scripts import or call directly.

| Dependency | Version | Used for |
|---|---|---|
| python | 3.12.16 | Runtime for all pipeline scripts |
| scanpy | 1.12.5 | Core single-cell toolkit: QC metrics, normalisation, HVG selection, Scrublet doublet scores, kNN graph, Leiden, UMAP, Wilcoxon marker tests, dotplots and UMAP plots |
| anndata | 0.13.4 | Cells × genes data container; reading and writing `.h5ad` files |
| harmonypy (PyPI) | 2.1.0 | Harmony batch correction of PCA embeddings across the 22 libraries |
| python-igraph | 1.0.0 | Leiden clustering (scanpy `flavor="igraph"`) |
| scikit-misc | 0.5.3 | LOESS fit required by scanpy's `seurat_v3` HVG method |
| scikit-learn | 1.9.1 | PCA, adjusted Rand index (seed stability), nearest-neighbour label transfer, base class for the custom kNN transformer |
| pynndescent | 0.7.0 | Approximate kNN search, wrapped in `scripts/knn.py` (scanpy's built-in pynndescent path hit a numba typing error) |
| numpy | 2.5.4 | Array maths throughout (knee detection, scaling, scores, kNN voting) |
| scipy | 1.18.1 | Reading Matrix Market files (`scipy.io.mmread`) and sparse matrices |
| pandas | 3.1.1 | Tables: metadata, QC summaries, crosstabs, all CSV outputs |
| pyarrow | 26.0.1 | Arrow-backed string columns used by pandas 3 |
| matplotlib | 3.11.2 | Barcode-rank diagnostic plot; backend for scanpy figures |

## Transitive dependencies the pipeline relies on

Not declared, because they come in through scanpy, anndata or python-igraph, but they do work in
the pipeline.

| Dependency | Version | Used for |
|---|---|---|
| umap-learn | 0.6.0 | Fuzzy kNN connectivities and UMAP embeddings (via scanpy) |
| numba | 0.68.1 | JIT compilation underlying umap-learn and pynndescent |
| igraph (C library) | 1.0.1 | Compiled graph library behind python-igraph; runs the Leiden algorithm |
| h5py | 3.16.0 | HDF5 I/O behind anndata's `.h5ad` files |
| seaborn | 0.13.2 |   Scanpy depends on it