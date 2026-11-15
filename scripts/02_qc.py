"""Step 2: per-library QC, hashtag demultiplexing, doublet detection; write merged counts.

Decisions:
- Only GEX features are used for clustering. ADT features are all-zero (no ADT library on GEO).
  Cellhash features are used only for donor demux + cross-donor doublet calls.
- QC thresholds are per-library and data-driven (MAD-based), because libraries differ a lot
  in depth (median UMI 1.2k-5.6k), with absolute floors.
- Hash demux: per library, per hashtag, threshold = median + 5 MAD of log1p(counts) among cells
  dominated by another hashtag. 1 positive -> singlet (donor), >=2 -> hash doublet, 0 -> negative.
  Hash-negatives are kept (they are often just low-depth cells), hash doublets removed.
- Scrublet (scanpy implementation) run per library; doublet_score > 0.25 removed.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp

LIBS, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sc.settings.verbosity = 1
SCRUB_THR = 0.25

def mad(x):
    return np.median(np.abs(x - np.median(x))) * 1.4826

def hash_demux(H, names, k=5.0):
    """H: cells x hashtags counts (only hashtags used in this library).
    Background for hashtag j is estimated from cells whose dominant hashtag is NOT j
    (robust to strong donor imbalance, unlike a per-hashtag 2-component mixture).
    Threshold_j = median + k*MAD of log1p(count_j) in that background (k=5: conservative,
    since the background itself contains ~5-10% doublets)."""
    L = np.log1p(H)
    top = L.argmax(1)
    thr = np.zeros(L.shape[1])
    for _ in range(2):  # refine once: recompute dominant call among positives
        for j in range(L.shape[1]):
            bg = L[top != j, j]
            thr[j] = np.median(bg) + k * max(mad(bg), 0.25)
        pos = L > thr
    npos = pos.sum(1)
    call = np.where(npos == 1, np.array(names)[pos.argmax(1)], np.where(npos == 0, "Negative", "Doublet"))
    return call, thr

adatas = []; summ = []
for f in sorted(LIBS.glob("*.h5ad")):
    a = ad.read_h5ad(f)
    lib = a.obs.library.iloc[0]
    hmask = a.var_names.str.startswith("Cellhash")
    H = a[:, hmask].X.toarray()
    used = H.mean(0) > 20                      # hashtags actually used in this library
    hnames = [n.replace("Cellhash", "donor") for n in a.var_names[hmask][used]]
    a.obs["hash_call"], hthr = hash_demux(H[:, used], hnames)
    a.obs["hash_total"] = H[:, used].sum(1)
    a = a[:, a.var.feature_type == "Gene Expression"].copy()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    a.var["ribo"] = a.var_names.str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True)
    o = a.obs
    lc, lg = o.log1p_total_counts, o.log1p_n_genes_by_counts
    lo_c = max(np.median(lc) - 3 * mad(lc), np.log1p(500))
    hi_c = np.median(lc) + 4 * mad(lc)
    lo_g = max(np.median(lg) - 3 * mad(lg), np.log1p(250))
    mt_thr = min(max(np.median(o.pct_counts_mt) + 4 * mad(o.pct_counts_mt), 10), 25)
    qc_pass = (lc >= lo_c) & (lc <= hi_c) & (lg >= lo_g) & (o.pct_counts_mt <= mt_thr)
    a.obs["qc_pass"] = qc_pass.values
    # Scrublet on QC-passing cells
    b = a[qc_pass.values].copy()
    # fixed score threshold (scrublet's automatic bimodality threshold was unstable across
    # libraries: 0-1000 calls); 0.25 is calibrated against hash doublets in 02b_doublet_check
    sc.pp.scrublet(b, expected_doublet_rate=0.08, threshold=SCRUB_THR, random_state=0)
    a.obs["doublet_score"] = np.nan; a.obs["scrublet_doublet"] = False
    a.obs.loc[b.obs_names, "doublet_score"] = b.obs.doublet_score
    a.obs.loc[b.obs_names, "scrublet_doublet"] = b.obs.predicted_doublet.astype(bool)
    keep = a.obs.qc_pass & ~a.obs.scrublet_doublet & (a.obs.hash_call != "Doublet")
    s = dict(library=lib, n_in=a.n_obs, qc_fail=int((~a.obs.qc_pass).sum()),
             min_umi=float(np.expm1(lo_c)), max_umi=float(np.expm1(hi_c)), min_genes=float(np.expm1(lo_g)),
             mt_thr=float(mt_thr), hash_thr=np.round(np.expm1(hthr)).tolist(),
             hash_doublet=int(((a.obs.hash_call == "Doublet") & a.obs.qc_pass).sum()),
             scrublet_doublet=int(a.obs.scrublet_doublet.sum()),
             hash_negative=int(((a.obs.hash_call == "Negative") & keep).sum()),
             n_out=int(keep.sum()))
    s.update({f"frac_{k}": round(v, 3) for k, v in a.obs.hash_call[keep].value_counts(normalize=True).items()})
    print(s, flush=True); summ.append(s)
    a.obs.to_csv(OUT / f"obs_{lib}.csv.gz")
    a = a[keep.values].copy()
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    adatas.append(a)

pd.DataFrame(summ).to_csv(OUT / "qc_summary.csv", index=False)
adata = ad.concat(adatas, join="outer", merge="first")
adata.obs["run"] = np.where(adata.obs.library.str.contains("220524"), "run_220524", "run_main")
sc.pp.filter_genes(adata, min_cells=20)
print(adata)
adata.write_h5ad(OUT / "merged_qc.h5ad", compression="gzip")
