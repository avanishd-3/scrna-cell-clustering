"""Step 5b: refine compartment sub-clusterings (re-uses the graphs computed in 05_subcluster.py).

Why: the first rule ("highest resolution with seed-ARI >= 0.75") over-split small compartments
(Plasma 3.6k cells -> 19 clusters, pDC 2.2k -> 12, Erythroid -> 24), and the DE merge test was
fooled by clonal / isotype Ig genes and immediate-early (dissociation) genes.

Refined rule, applied uniformly to all compartments:
  1. resolution = argmax seed-ARI over res >= 0.2, provided that ARI >= 0.75; else res = 0.1.
     (prefer the most reproducible partition, not the finest acceptable one)
  2. Leiden (seed 0) at that resolution on the stored Harmony kNN graph.
  3. clusters smaller than max(30, 1% of compartment) are merged into the nearest cluster
     (Harmony-centroid distance).
  4. iterative merge of nearest pairs with < 10 DE genes (|log2FC|>1, padj<0.01, detection
     difference >= 20 pp), where Ig/TCR, mito, ribo, sex-linked and immediate-early stress genes
     are ignored.
usage: 05b_refine.py sub_dir sweep.csv
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

SUB, SWEEP = Path(sys.argv[1]), sys.argv[2]
sweep = pd.read_csv(SWEEP).drop_duplicates(["compartment", "resolution"])
IGNORE = re.compile(r"^(IG[HKL][VDJC]|IGHM$|IGHD$|IGHG|IGHA|IGHE|IGLL5$|JCHAIN$|TR[ABDG][VDJC]|MT-|RP[SL]\d|"
                    r"XIST$|TSIX$|RPS4Y|DDX3Y$|UTY$|KDM5D$|EIF1AY$|USP9Y$|ZFY$|MALAT1$|NEAT1$|"
                    r"FOS$|FOSB$|JUN$|JUNB$|JUND$|ATF3$|EGR[1-3]$|IER[235]$|DUSP1$|ZFP36$|ZFP36L1$|KLF[246]$|"
                    r"NR4A[1-3]$|HSPA1[AB]$|HSPA8$|HSPH1$|DNAJB1$|PPP1R15A$|CD69$|RGS[12]$|TRIB1$|BTG2$|SOCS3$)")
MIN_DE = 10

def n_de(a, c1, c2, key, genes):
    s = a[a.obs[key].isin([c1, c2])].raw.to_adata()[:, genes]
    s.obs[key] = a.obs.loc[s.obs_names, key].values
    sc.tl.rank_genes_groups(s, key, groups=[c1], reference=c2, method="wilcoxon")
    d = sc.get.rank_genes_groups_df(s, c1).set_index("names")
    in1 = (s.obs[key] == c1).values
    p1 = pd.Series(np.asarray((s.X[in1] > 0).mean(0)).ravel(), index=s.var_names)
    p2 = pd.Series(np.asarray((s.X[~in1] > 0).mean(0)).ravel(), index=s.var_names)
    dp = (p1 - p2).reindex(d.index)
    sig = (d.pvals_adj < 0.01) & (((d.logfoldchanges > 1) & (dp >= 0.2)) | ((d.logfoldchanges < -1) & (dp <= -0.2)))
    return int(sig.sum())

def nearest(a, cl):
    cen = np.stack([a.obsm["X_harmony"][a.obs["sub"].values == c].mean(0) for c in cl])
    D = np.linalg.norm(cen[:, None] - cen[None], axis=2); np.fill_diagonal(D, np.inf)
    return D

labels, log = [], []
for f in sorted(SUB.glob("*.h5ad")):
    name = f.stem
    a = ad.read_h5ad(f)
    sw = sweep[sweep.compartment == name]
    cand = sw[(sw.resolution >= 0.2) & (sw.mean_ARI >= 0.75)]
    res = float(cand.loc[cand.mean_ARI.idxmax(), "resolution"]) if len(cand) else 0.1
    sc.tl.leiden(a, resolution=res, key_added="sub", random_state=0, flavor="igraph", n_iterations=2)
    a.obs["sub"] = a.obs["sub"].astype(str)
    k0 = a.obs["sub"].nunique()
    genes = a.raw.var_names[~a.raw.var_names.str.match(IGNORE)]
    min_size = max(30, int(0.01 * a.n_obs))
    # 3. size merge
    while True:
        vc = a.obs["sub"].value_counts(); cl = list(vc.index)
        small = [c for c in cl if vc[c] < min_size]
        if not small or len(cl) < 2: break
        D = nearest(a, cl); i = cl.index(small[-1]); j = int(np.argmin(D[i]))
        log.append(dict(compartment=name, action="size_merge", src=cl[i], dst=cl[j], n=int(vc[cl[i]])))
        a.obs.loc[a.obs["sub"] == cl[i], "sub"] = cl[j]
    # 4. DE merge
    while True:
        cl = list(a.obs["sub"].unique())
        if len(cl) < 2: break
        D = nearest(a, cl); merged = False
        for i in np.argsort(D.min(1)):
            j = int(np.argmin(D[i])); nd = n_de(a, cl[i], cl[j], "sub", genes)
            if nd < MIN_DE:
                log.append(dict(compartment=name, action="de_merge", src=cl[i], dst=cl[j], n_de=nd))
                a.obs.loc[a.obs["sub"] == cl[i], "sub"] = cl[j]; merged = True; break
        if not merged: break
    order = a.obs["sub"].value_counts().index
    a.obs["sub"] = a.obs["sub"].map({c: str(i) for i, c in enumerate(order)})
    a.obs["subcluster"] = name + "_" + a.obs["sub"]
    print(f"{name}: n={a.n_obs} res={res} leiden_k={k0} final_k={a.obs['sub'].nunique()}", flush=True)
    log.append(dict(compartment=name, action="final", resolution=res, leiden_k=k0, final_k=a.obs["sub"].nunique()))
    a.write_h5ad(f, compression="gzip")
    labels.append(a.obs[["subcluster"]])

pd.concat(labels).to_csv(SUB / "subcluster_labels.csv")
pd.DataFrame(log).to_csv(SUB / "refine_log.csv", index=False)
