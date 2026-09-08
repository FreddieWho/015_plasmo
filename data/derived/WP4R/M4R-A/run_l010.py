#!/usr/bin/env python3
"""L-010: 真 GO 人工注释下 AP2 富集复核（CLM09 方法上限）。
QuickGO taxon 36329 全量注释下载 -> UniProt->PF3D7 映射（NCBI GFF Dbxref）->
预指定 GO 集 Fisher（与 round-2 同口径）-> IEA/非IEA 分层。
seed 20260905; sc env.
"""
import gzip, json, math, os, re, sys, time
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats

SEED = 20260905
ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4R/M4R-A"
RAW = f"{OUT}/L010_quickgo_pages.jsonl"
os.makedirs(OUT, exist_ok=True)

# 预指定 GO 集（标签运行时经 QuickGO term API 核对，不背号码即采信原文）
GOSETS = {
    "TF_DNAbinding": ["GO:0003700", "GO:0140110"],
    "RNA_binding": ["GO:0003723"],
    "chromatin": ["GO:0003682", "GO:0006338"],
    "reg_transcription": ["GO:0006355"],
}

def get(url, accept="application/json", timeout=60, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": accept})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:
            print(f"  retry {i+1}: {type(e).__name__}", flush=True)
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"fetch failed: {url}")

# 0. 核对 GO 标签
print("== GO term label check ==", flush=True)
for sname, goids in GOSETS.items():
    for g in goids:
        d = json.loads(get(f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{g}"))
        name = d["results"][0].get("name")
        print(f"  {g} = {name}  [set {sname}]", flush=True)

# 1. 下载全量注释
print("== download QuickGO taxon 36329 ==", flush=True)
d0 = json.loads(get("https://www.ebi.ac.uk/QuickGO/services/annotation/search?taxonId=36329&page=1&limit=100"))
total = d0["numberOfHits"]
npages = (total + 99) // 100
print(f"  hits={total} pages={npages}", flush=True)
rows = []
t0 = time.time()
import sys as _sys
if "--reuse-raw" in _sys.argv and os.path.exists(RAW):
    rows = [json.loads(l) for l in open(RAW)]
    print(f"  reused {RAW} n={len(rows)}", flush=True)
else:
    for p in range(1, npages + 1):
        d = json.loads(get(f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?taxonId=36329&page={p}&limit=100"))
        rows.extend(d["results"])
        if p % 25 == 0 or p == npages:
            print(f"  page {p}/{npages} rows={len(rows)} t={time.time()-t0:.0f}s", flush=True)
    with open(RAW, "w") as h:
        for r in rows:
            h.write(json.dumps(r) + "\n")
    print(f"  saved {RAW} n={len(rows)}", flush=True)

# 2. UniProt -> PF3D7 映射（NCBI GFF Dbxref）
print("== map UniProt->PF3D7 ==", flush=True)
gff = f"{ROOT}/data/raw/ncbi-datasets/GCF_000002765.6/ncbi_dataset/data/GCF_000002765.6/genomic.gff"
up2pf, pf_check = {}, set()
with open(gff) as h:
    for line in h:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        attr = c[8]
        m1 = re.search(r"UniProtKB/(?:TrEMBL|Swiss-Prot):([A-Z0-9]+)", attr)
        m2 = re.search(r"locus_tag=([A-Z0-9_]+)", attr)
        if m1 and m2 and m2.group(1).startswith("PF3D7_"):
            up2pf[m1.group(1)] = m2.group(1)
            pf_check.add(m2.group(1))
print(f"  map entries={len(up2pf)} pf genes={len(pf_check)}", flush=True)
assert len(pf_check) > 5000, "mapping too small, GFF wrong?"

# 3. 解析注释 -> 每基因 GO 集
recs = []
for r in rows:
    up = (r.get("geneProductId") or "").split(":")[-1]
    pf = up2pf.get(up)
    if not pf:
        continue
    quals = r.get("qualifier") or ""
    if "NOT" in quals.split("|"):
        continue
    recs.append(dict(gene=pf, go=r.get("goId"), ev=r.get("goEvidence") or "?",
                      aspect=r.get("goAspect") or "?"))
A = pd.DataFrame(recs)
A.to_csv(f"{OUT}/M4RA_l010_annotations.tsv", sep="\t", index=False)
print(f"  mapped annotations={len(A)} genes={A.gene.nunique()} ev={A.ev.value_counts().to_dict()}", flush=True)

# 4. 富集（与 round-2 同口径：universe 5285，Asn q90）
gff2 = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff"
gids = set()
with open(gff2) as _h:
    for _line in _h:
        if _line.startswith("#"):
            continue
        _c = _line.rstrip("\n").split("\t")
        if len(_c) > 8 and _c[2] == "protein_coding_gene":
            _m = re.search(r"ID=(PF3D7_\d+)", _c[8])
            if _m:
                gids.add(_m.group(1))
print(f"  gff genes={len(gids)}", flush=True)
FAA = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta"
def _gkey(h):
    import re as _re
    m = _re.search(r"PF3D7_\d+", str(h))
    return m.group(0) if m else None
prots = set()
hid = None
with open(FAA) as fh:
    for line in fh:
        if line.startswith(">"):
            g = _gkey(line.rstrip())
            if g:
                prots.add(g)
F = pd.read_csv(f"{OUT}/../M4R-D/M4RD_D5_pf_features.tsv", sep="\t", usecols=["gene", "asn_frac"])
F = F.drop_duplicates(subset="gene", keep="first")
F = F[F.gene.isin(prots) & F.gene.isin(gids)].copy()
assert len(F) == 5285, len(F)
thr = F.asn_frac.quantile(0.90)
print(f"  universe check: n={len(F)} q90={thr:.4f} (round-2: 0.1837)", flush=True)
assert abs(thr - 0.1837) < 1e-3, thr
top = set(F.loc[F.asn_frac >= thr, "gene"])
print(f"  universe={len(F)} q90={thr:.4f} top10 n={len(top)}", flush=True)

def fisher(geneset):
    s = set(geneset) & set(F.gene)
    a = len(s & top)
    b = len(s - top)
    c = len(top - s)
    d = len(F) - len(s | top)
    or_, p = stats.fisher_exact([[a, b], [c, d]])
    # Woolf CI
    if min(a, b, c, d) > 0:
        se = math.sqrt(1/a + 1/b + 1/c + 1/d)
        lo, hi = math.exp(math.log(or_) - 1.96*se), math.exp(math.log(or_) + 1.96*se)
    else:
        lo, hi = float("nan"), float("nan")
    return dict(n=len(s), in_top=a, OR=or_, lo=lo, hi=hi, p=p)

out = []
for sname, goids in GOSETS.items():
    g = A[A.go.isin(goids)]
    r = fisher(g.gene.unique()); r.update(set=sname, stratum="all", n_ann=len(g))
    out.append(r)
    for strat, evs in [("nonIEA", None), ("IEA", None)]:
        gg = g[g.ev != "IEA"] if strat == "nonIEA" else g[g.ev == "IEA"]
        r2 = fisher(gg.gene.unique()); r2.update(set=sname, stratum=strat, n_ann=len(gg))
        out.append(r2)
R = pd.DataFrame(out)
ps = R.p.fillna(1.0).values
order = np.argsort(ps)
q = np.empty_like(ps); q[order] = np.minimum.accumulate((ps[order] * len(ps) / (np.arange(len(ps)) + 1))[::-1])[::-1]
R["q_bh"] = np.clip(q, 0, 1)
R.to_csv(f"{OUT}/M4RA_l010_enrichment.tsv", sep="\t", index=False)
print(R[["set", "stratum", "n", "in_top", "OR", "p", "q_bh"]].to_string(), flush=True)
print("DONE", flush=True)
