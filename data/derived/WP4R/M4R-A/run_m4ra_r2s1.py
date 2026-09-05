#!/usr/bin/env python
"""M4R-A ROUND2 step1: NCBI GFF -> gene InterPro map + domain IPR resolution from data."""
import hashlib, os, re
import pandas as pd

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4R/M4R-A"
GFF = f"{ROOT}/data/raw/ncbi-datasets/GCF_000002765.6/ncbi_dataset/data/GCF_000002765.6/genomic.gff"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

gene_ipr = {}
gene_product = {}
n_cds = 0
with open(GFF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        n_cds += 1
        a = c[8]
        m = re.search(r"locus_tag=(PF3D7_\d+)", a)
        if not m:
            continue
        g = m.group(1)
        iprs = set(re.findall(r"InterPro:(IPR\d+)", a))
        gene_ipr.setdefault(g, set()).update(iprs)
        pm = re.search(r"product=([^;]+)", a)
        if pm and g not in gene_product:
            gene_product[g] = pm.group(1).replace("%2C", ",")
print("cds_lines:", n_cds, "genes_with_ipr:", len(gene_ipr), flush=True)

rows = [{"gene": g, "interpro": ";".join(sorted(v)), "n_ipr": len(v),
         "prod": gene_product.get(g, "")} for g, v in gene_ipr.items()]
M = pd.DataFrame(rows)
M.to_csv(f"{OUT}/M4RA_interpro_map.tsv", sep="\t", index=False)
print("map rows:", len(M), "genes with >=1 IPR:",
      int((M.n_ipr > 0).sum()), flush=True)

# resolve domain IPRs from data: which IPRs co-occur with AP2 / pumilio product text
low = M["prod"].str.lower()
is_ap2 = low.str.contains("ap2|apetala", na=False)
is_pum = low.str.contains("pumilio|puf", na=False)
print("product-AP2 genes:", int(is_ap2.sum()), "product-PUM genes:", int(is_pum.sum()), flush=True)
from collections import Counter
ca, cp = Counter(), Counter()
for s in M[is_ap2].interpro: ca.update(s.split(";") if s else [])
for s in M[is_pum].interpro: cp.update(s.split(";") if s else [])
print("AP2-product IPRs:", ca.most_common(12), flush=True)
print("PUM-product IPRs:", cp.most_common(12), flush=True)
# background frequency of those IPRs
allc = Counter()
for s in M.interpro:
    if s: allc.update(s.split(";"))
cand = set([k for k, _ in ca.most_common(12)] + [k for k, _ in cp.most_common(12)])
for k in sorted(cand):
    print(f"{k}: bg_n={allc.get(k,0)} ap2_n={ca.get(k,0)} pum_n={cp.get(k,0)}", flush=True)
print("STEP1 DONE", flush=True)
