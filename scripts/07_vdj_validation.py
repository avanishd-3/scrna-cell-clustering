"""Step 7: orthogonal validation of B-lineage clusters with BCR (V(D)J) contigs.
Ig genes were excluded from clustering features, so receptor status is independent evidence.
Expected: pro-B ~no productive IGH; pre-B productive IGH without light chain; immature/mature
paired IGH + IGK/IGL; non-B cells ~no productive contigs.
usage: 07_vdj_validation.py vdj_dir cell_clusters.csv out_csv
"""
import sys, glob, re
from pathlib import Path
import pandas as pd, numpy as np
VDJ, CC, OUT = sys.argv[1:4]
rows = []
for f in glob.glob(f"{VDJ}/*_filtered_contig_annotations.csv.gz"):
    lib = re.match(r"GSM\d+_(.+)_filtered_contig_annotations", Path(f).name).group(1)
    d = pd.read_csv(f, usecols=["barcode", "is_cell", "chain", "productive", "high_confidence", "c_gene"])
    d = d[d.high_confidence.astype(str).str.lower() == "true"]
    d["is_prod"] = d.productive.astype(str).str.lower() == "true"
    d["h"] = (d.chain == "IGH") & d.is_prod
    d["l"] = d.chain.isin(["IGK", "IGL"]) & d.is_prod
    d["hc"] = d.c_gene.where(d.h)
    g = d.groupby("barcode")
    s = pd.DataFrame({"any_contig": g.size() > 0, "IGH_prod": g.h.any(), "L_prod": g.l.any(),
                      "IGH_c": g.hc.agg(lambda x: ",".join(sorted(set(x.dropna()))))})
    s.index = lib + "_" + s.index
    rows.append(s)
v = pd.concat(rows)
c = pd.read_csv(CC, index_col=0)
c = c.join(v)
for k in ["any_contig", "IGH_prod", "L_prod"]:
    c[k] = c[k].fillna(False).astype(bool)
c["state"] = np.select([c.IGH_prod & c.L_prod, c.IGH_prod & ~c.L_prod, ~c.IGH_prod & c.L_prod],
                       ["H+L", "H only", "L only"], "none")
c["switched"] = c.IGH_c.fillna("").str.contains("IGHG|IGHA|IGHE")
ref = c[c.label_source == "clustered"]
t = pd.crosstab(ref.cluster, ref.state, normalize="index").round(3)
t["n"] = ref.cluster.value_counts()
t["class_switched"] = ref.groupby("cluster").switched.mean().round(3)
t["annotation"] = ref.groupby("cluster").annotation.first()
t = t[["n", "none", "H only", "L only", "H+L", "class_switched", "annotation"]]
t.to_csv(OUT)
pd.set_option("display.width", 200)
sel = t.index.str.startswith(("B_lineage", "Plasma", "T_NK_0", "T_NK_1", "Monocyte_0", "HSPC_MyeloidProg_0", "HSPC_MyeloidProg_6", "HSPC_MyeloidProg_7"))
print(t[sel].to_string())
