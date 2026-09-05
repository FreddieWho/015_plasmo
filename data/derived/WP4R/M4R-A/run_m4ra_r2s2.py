#!/usr/bin/env python
"""M4R-A ROUND2 step2: InterPro-anchored formal sets + enrichment + logit, side-by-side vs regex."""
import math, os, re
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4R/M4R-A"
FAA = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta"
GFF = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff"

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s))
    return m.group(0) if m else None

# ---- protein features (same recipe as round-1) ----
prots = {}
with open(FAA) as fh:
    hid, seq = None, []
    def flush():
        if hid is None: return
        g = gkey(hid)
        if not g: return
        s = "".join(seq).upper().replace("*", "")
        if len(s) > len(prots.get(g, "")): prots[g] = s
    for line in fh:
        if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
        else: seq.append(line.strip())
    flush()

def feats(s):
    n = len(s); asn = s.count("N") / n if n else 0.0
    mt, cur = 0, 0
    for ch in s:
        cur = cur + 1 if ch == "N" else 0
        mt = max(mt, cur)
    w, th, low = 20, 1.8, 0
    if n >= w:
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        cnt = np.zeros(256, int)
        for c in arr[:w]: cnt[c] += 1
        def ent():
            p = cnt[cnt > 0] / w
            return float(-(p * np.log2(p)).sum())
        low = int(ent() < th)
        for i in range(w, n):
            cnt[arr[i]] += 1; cnt[arr[i - w]] -= 1
            p = cnt[cnt > 0] / w
            if float(-(p * np.log2(p)).sum()) < th: low += 1
        lcr = low / (n - w + 1)
    else:
        lcr = 0.0
    return n, asn, mt, lcr

P = pd.DataFrame({g: feats(s) for g, s in prots.items()},
                 index=["length", "asn_frac", "polyn_max", "lcr_frac"]).T
P.index.name = "gene"; P = P.reset_index()

genes, prod = {}, {}
for line in open(GFF):
    if line.startswith("#") or "\tprotein_coding_gene\t" not in line: continue
    f = line.rstrip("\n").split("\t")
    m = re.search(r"ID=(PF3D7_\d+)", f[8]); d = re.search(r"description=([^;]+)", f[8])
    if m:
        g = m.group(1)
        genes[g] = True
        prod[g] = (d.group(1) if d else "").replace("%2C", ", ").replace("%27", "'").lower()
D = P.merge(pd.DataFrame({"gene": list(genes), "prod": [prod.get(g, "") for g in genes]}),
            on="gene", how="inner")
q90 = D.asn_frac.quantile(0.9)
D["asn_top10"] = D.asn_frac >= q90
print("universe:", len(D), "q90:", round(float(q90), 4), flush=True)

M = pd.read_csv(f"{OUT}/M4RA_interpro_map.tsv", sep="\t")
M["iprset"] = M.interpro.fillna("").apply(lambda s: set(s.split(";")) if s else set())
g2ipr = dict(zip(M.gene, M.iprset))

# ---- regex classes (round-1 recipe, for anchor discovery only) ----
SETS = {
    "ApiAP2": r"\bap2|apiap2",
    "PUF": r"puf1|puf2|pumilio",
    "RNA_binding_broad": r"rna[- ]binding|rrm|helicase|zinc finger|alba|dozi|cith|caf\d|not\d|ccr4",
    "chromatin_reg": r"chromatin|histone|bromodomain|set domain|gcn5|sir2|hdac|acetyltransferase|methyltransferase|morne?c|hp1",
    "CCR4_NOT": r"ccr4|caf1|not1|not2|not5|pop2|caf40",
    "transcription_reg": r"transcription factor|transcriptional|rna polymerase|general transcription|mediator",
    "sexual_gametocyte": r"gametocyte|gamete|oocyst|ookinete|zygote|fertili|p47|p48/45|ap2-g|ap2g|gdv1|mdv1",
    "proteostasis": r"heat shock|chaperone|dnaj|hsp|proteasome|ubiquitin|unfolded",
}
for k, pat in SETS.items():
    D[k + "_rx"] = D["prod"].str.contains(pat, regex=True)

# ---- anchor IPRs per class: per-IPR Fisher(class vs rest), BH ----
all_iprs = sorted({i for v in g2ipr.values() for i in v})
print("distinct IPRs:", len(all_iprs), flush=True)
D["iprset"] = D.gene.map(lambda g: g2ipr.get(g, set()))
anchors = {}
for k in SETS:
    hit = D[k + "_rx"].values
    if hit.sum() < 3:
        anchors[k] = set(); continue
    cand = sorted({i for v in D[hit].iprset for i in v})
    ps, ors, ns = [], [], []
    for ipr in cand:
        car = D.iprset.map(lambda s, i=ipr: i in s).values
        a = int((hit & car).sum()); b = int(((~hit) & car).sum())
        c = int((hit & ~car).sum()); d = int(((~hit) & ~car).sum())
        if min(a, b, c, d) == 0:
            ps.append(1.0); ors.append(float("nan")); ns.append(a); continue
        o, p = stats.fisher_exact([[a, b], [c, d]])
        ps.append(p); ors.append(o); ns.append(a)
    q = bh(np.array(ps))
    keep = {ipr for ipr, qq, o, nn in zip(cand, q, ors, ns)
            if qq < 0.05 and (o > 2) and nn >= 3}
    anchors[k] = keep
    print(f"{k}: regex_n={int(hit.sum())} anchor_IPRs={len(keep)}", flush=True)

# overrides from data-resolved domains
anchors["ApiAP2"] = {"IPR001471"}          # 24/24 product-AP2, bg 24: clean
anchors["PUF_strict"] = {"IPR001313"}      # 2/3 product-PUM, bg 2
# Pumilio-repeat broader (explicit caveat: generic repeat IPRs)
anchors["PUF_repeatbroad"] = {"IPR001313", "IPR011989", "IPR016024", "IPR033133"}

aset_rows = []
for k, aiprs in anchors.items():
    D[k + "_ip"] = D.iprset.map(lambda s: bool(s & aiprs))
    aset_rows.append({"formal_set": k, "anchor_iprs": ";".join(sorted(aiprs)),
                      "n": int(D[k + "_ip"].sum())})
pd.DataFrame(aset_rows).to_csv(f"{OUT}/M4RA_interpro_sets.tsv", sep="\t", index=False)
print(pd.DataFrame(aset_rows).to_string(), flush=True)

# ---- enrichment: formal vs regex side-by-side ----
top = D.asn_top10.values
erows = []
def enrich(col):
    h = D[col].values
    a = int((h & top).sum()); b = int(((~h) & top).sum())
    c = int((h & ~top).sum()); d = int(((~h) & ~top).sum())
    o, p = stats.fisher_exact([[a, b], [c, d]])
    if min(a, b, c, d) > 0:
        se = math.sqrt(1/a + 1/b + 1/c + 1/d)
        lo, hi = math.exp(math.log(o) - 1.96*se), math.exp(math.log(o) + 1.96*se)
    else:
        lo, hi = (0.0, float("inf")) if o in (0, float("inf")) else (float("nan"), float("nan"))
    return o, lo, hi, p, int(h.sum()), a

tests = []
for k in list(anchors):
    o, lo, hi, p, n, a = enrich(k + "_ip")
    erows.append({"set": k + " [InterPro-formal]", "n_set": n, "n_in_top10": a,
                  "OR": o, "lo": lo, "hi": hi, "p": p})
    tests.append(p)
    if k in SETS:
        o2, lo2, hi2, p2, n2, a2 = enrich(k + "_rx")
        erows.append({"set": k + " [regex-baseline]", "n_set": n2, "n_in_top10": a2,
                      "OR": o2, "lo": lo2, "hi": hi2, "p": p2})
        tests.append(p2)
E = pd.DataFrame(erows)
E["q_bh"] = bh(E.p.values)
E.to_csv(f"{OUT}/M4RA_interpro_enrichment.tsv", sep="\t", index=False)
print(E.to_string(), flush=True)

# ---- strict PUF table ----
puf = D[D["PUF_strict_ip"]][["gene", "prod", "length", "asn_frac", "polyn_max", "lcr_frac", "asn_top10"]].copy()
puf.to_csv(f"{OUT}/M4RA_strictPUF.tsv", sep="\t", index=False)
print("strictPUF n:", len(puf), flush=True)
print(puf.to_string(), flush=True)

# ---- logit counterfactual with InterPro regulator flag ----
D["loglen"] = np.log(D.length)
D["is_reg_ip"] = (D["ApiAP2_ip"] | D["PUF_strict_ip"] | D["chromatin_reg_ip"] | D["CCR4_NOT_ip"]).astype(int)
X = sm.add_constant(D[["is_reg_ip", "loglen", "lcr_frac"]])
mdl = sm.Logit(D.asn_top10.astype(int), X).fit(disp=0)
adj = pd.DataFrame({"coef": mdl.params, "se": mdl.bse, "OR": np.exp(mdl.params),
                    "OR_lo": np.exp(mdl.params - 1.96*mdl.bse),
                    "OR_hi": np.exp(mdl.params + 1.96*mdl.bse), "p": mdl.pvalues})
adj.to_csv(f"{OUT}/M4RA_interpro_logit.tsv", sep="\t")
print(adj.to_string(), flush=True)
print("STEP2 DONE", flush=True)
