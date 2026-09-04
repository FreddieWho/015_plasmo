#!/usr/bin/env python
# L1b: interpret the 5 double-positive codons (D-022: interpret AFTER significance)
# 1) position-1 composition confound: re-test with gc1 covariate added
# 2) driver genes with product descriptions (PlasmoDB GFF)
# 3) tRNA correspondence from literature annotation (quantitative supply quarantined)
import gzip, hashlib, json, os, re
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
RAW = f"{ROOT}/data/raw"
L1 = f"{ROOT}/data/derived/WP4/L1_codon_x_stress"
OUT = f"{ROOT}/data/derived/WP4/L1b_interpret"
os.makedirs(OUT, exist_ok=True)
SEED = 20260904
rng = np.random.default_rng(SEED)
B = 1000
HITS = ["AAT_freq", "GTT_freq", "GTC_freq", "GCT_freq", "GCC_freq"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# --- rebuild per-gene codon counts (same as L1) ---
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STD = dict(zip(CODONS, AA))
SENSE = [c for c in CODONS if STD[c] != "*"]

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s))
    return m.group(0) if m else None

acc = "GCF_000002765.6"
path = f"{RAW}/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna"
best = {}
with open(path) as fh:
    hid, seq = None, []
    def flush():
        if hid is None: return
        m = re.search(r"\[locus_tag=([^\]]+)\]", hid)
        if not m: return
        g = gkey(m.group(1))
        if not g: return
        s = "".join(seq).upper().replace("U", "T")
        if "[partial" in hid or len(s) < 150 or len(s) % 3 or not set(s) <= set("ATGC"): return
        cod = [s[i:i+3] for i in range(0, len(s), 3)]
        if "*" in [STD.get(c, "*") for c in cod][:-1]: return
        if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
    for line in fh:
        if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
        else: seq.append(line.strip())
    flush()
rows = {}
for g, (_, s) in best.items():
    cod = [s[i:i+3] for i in range(0, len(s), 3) if STD.get(s[i:i+3], "*") != "*"]
    if not cod: continue
    n = len(cod)
    from collections import Counter
    cc = Counter(cod)
    gc3 = sum(1 for c in cod if c[2] in "GC") / n
    gc1 = sum(1 for c in cod if c[0] in "GC") / n
    rows[g] = {"n_cod": n, "gc3": gc3, "gc1": gc1, **{c: cc.get(c, 0) / n for c in SENSE}}
CDF = pd.DataFrame(rows).T; CDF.index.name = "gene"
print("genes:", len(CDF), flush=True)

# --- bootstrap with gene clusters, configurable covariates ---
def boot_slope(x, C, y, groups, B=B):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(C).all(axis=1)
    x, C, y, groups = x[ok], C[ok], y[ok], np.asarray(groups)[ok]
    n = len(y)
    if n < 100 or np.std(x) == 0 or np.std(y) == 0: return np.nan, np.nan, n
    u, inv = np.unique(groups, return_inverse=True)
    def wb(w):
        sw = np.sqrt(w); X = np.column_stack([x * sw, C * sw[:, None], sw])
        b, *_ = np.linalg.lstsq(X, y * sw, rcond=None); return b[0]
    b0 = wb(np.ones(n)); bs = np.empty(B)
    for i in range(B):
        cnt = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u))
        try: bs[i] = wb(cnt[inv].astype(float))
        except Exception: bs[i] = np.nan
    bs = bs[np.isfinite(bs)]
    p = 2 * min((np.sum(bs <= 0) + 1) / (len(bs) + 1), (np.sum(bs >= 0) + 1) / (len(bs) + 1))
    return b0, p, n

def bh(p):
    p = np.asarray(p, float); ok = np.isfinite(p); q = np.full(len(p), np.nan)
    po = p[ok]; n = len(po)
    if n == 0: return q
    o = np.argsort(po); qq = np.empty(n)
    qq[o] = np.minimum.accumulate((po[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q[ok] = np.minimum(qq, 1); return q

RESP = pd.read_csv(f"{L1}/L1_response_by_gene.tsv", sep="\t")
DTE = pd.read_csv(f"{L1}/L1_GSE226632_dTE.tsv", sep="\t")

# --- 1) pos1-adjusted retest ---
out = []
for f in HITS:
    cod = f.replace("_freq", "")
    fv = CDF[cod]
    for bg in ["Dd2_WT", "Dd2_R539T", "Dd2_C580Y"]:
        sub = RESP[(RESP.bg == bg) & (RESP.time <= 6) & RESP.usable]
        j = sub.merge(CDF[[cod, "gc3", "gc1"]].reset_index().rename(columns={cod: "f"}), on="gene")
        x = j.f.values; y = j.resp.values
        b3, p3, n = boot_slope(x, np.column_stack([j.gc3.values, j.time.values]), y, j.gene.values)
        b13, p13, _ = boot_slope(x, np.column_stack([j.gc1.values, j.gc3.values, j.time.values]), y, j.gene.values)
        out.append({"feature": f, "layer": "discovery_DHA", "bg": bg, "slope_gc3only": b3,
                    "p_gc3only": p3, "slope_gc1gc3": b13, "p_gc1gc3": p13, "n": n})
    j = DTE.merge(CDF[[cod, "gc3", "gc1"]].reset_index().rename(columns={cod: "f"}), on="gene")
    b3, p3, n = boot_slope(j.f.values, j.gc3.values.reshape(-1, 1), j.dTE.values, j.gene.values)
    b13, p13, _ = boot_slope(j.f.values, j[["gc1", "gc3"]].values, j.dTE.values, j.gene.values)
    out.append({"feature": f, "layer": "replication_dTE", "bg": "-", "slope_gc3only": b3,
                "p_gc3only": p3, "slope_gc1gc3": b13, "p_gc1gc3": p13, "n": n})
R1 = pd.DataFrame(out)
for layer in R1.layer.unique():
    m = R1.layer == layer
    R1.loc[m, "q_gc1gc3"] = bh(R1.loc[m, "p_gc1gc3"].values)
R1["survives_pos1"] = (R1.q_gc1gc3 < 0.05) & (np.sign(R1.slope_gc1gc3) == np.sign(R1.slope_gc3only))
R1.to_csv(f"{OUT}/L1b_pos1_adjusted.tsv", sep="\t", index=False)
print("pos1 retest done", flush=True)

# --- 2) driver genes with PlasmoDB product descriptions ---
gff = f"{RAW}/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff"
prod = {}
with open(gff) as fh:
    for line in fh:
        if line.startswith("#") or "\tgene\t" not in line: continue
        parts = line.rstrip("\n").split("\t")
        attrs = parts[8]
        m = re.search(r"ID=([^;]+)", attrs); d = re.search(r"description=([^;]+)", attrs)
        if m:
            g = gkey(m.group(1))
            if g: prod[g] = (d.group(1) if d else "").replace("%2C", ",").replace("%3B", ";")
print("GFF products:", len(prod), flush=True)
drv = []
for f in HITS:
    cod = f.replace("_freq", "")
    resp_wt = RESP[(RESP.bg == "Dd2_WT") & (RESP.time <= 6) & RESP.usable].groupby("gene").resp.mean()
    d = pd.DataFrame({"f": CDF[cod], "resp": resp_wt}).dropna()
    d["score"] = d.f.rank(pct=True) * d.resp.abs().rank(pct=True)
    top = d.nlargest(25, "score")
    for g, r in top.iterrows():
        drv.append({"feature": f, "gene": g, "feature_value": r.f, "mean_resp_Dd2WT": r.resp,
                    "product": prod.get(g, "")})
pd.DataFrame(drv).to_csv(f"{OUT}/L1b_driver_genes.tsv", sep="\t", index=False)

# --- 3) tRNA correspondence (literature annotation; quantitative supply quarantined) ---
trna = pd.DataFrame([
    {"feature": "AAT_freq", "aa": "Asn", "codon": "AAT", "decoder": "tRNA-Asn-GUU (wobble reads AAT/AAC)",
     "known_mods": "no s2U (not U34 Lys/Glu/Gln class); queuosine at U34 of Asn tRNA reported in other eukaryotes",
     "lit_anchor": "Small-Saunders 2024 NatMicrobiol (s2U class); Q-tRNA literature"},
    {"feature": "GTT_freq", "aa": "Val", "codon": "GTT", "decoder": "tRNA-Val-GAC / IAC",
     "known_mods": "inosine I34 class reads T/C-ending; no s2U", "lit_anchor": "standard wobble annotation"},
    {"feature": "GTC_freq", "aa": "Val", "codon": "GTC", "decoder": "tRNA-Val-GAC",
     "known_mods": "same Val family", "lit_anchor": "standard wobble annotation"},
    {"feature": "GCT_freq", "aa": "Ala", "codon": "GCT", "decoder": "tRNA-Ala-CGC/GGC family",
     "known_mods": "no s2U; GC-start codons, low abundance in AT-rich Pf", "lit_anchor": "standard wobble annotation"},
    {"feature": "GCC_freq", "aa": "Ala", "codon": "GCC", "decoder": "tRNA-Ala-GGC",
     "known_mods": "same Ala family", "lit_anchor": "standard wobble annotation"},
])
trna.to_csv(f"{OUT}/L1b_tRNA_annotation.tsv", sep="\t", index=False)

# --- deliverables ---
surv = R1[R1.layer != "replication_dTE"].groupby("feature").survives_pos1.sum()
surv_r = R1[R1.layer == "replication_dTE"].set_index("feature").survives_pos1
summ = {f: {"discovery_bgs_surviving_pos1": int(surv.get(f, 0)), "replication_survives_pos1": bool(surv_r.get(f, False))}
        for f in HITS}
json.dump(summ, open(f"{OUT}/L1b_qc.json", "w"), indent=2)
open(f"{OUT}/claim_impact.md", "w").write(f"""# L1b claim impact addendum (interpretation of EVID-L1-002 double-positives)

- date: 2026-09-04
- pos1-composition adjustment result: `{json.dumps(summ, indent=1)}`
- interpretation: features surviving gc1+gc3 covariates in BOTH layers are beyond-composition signals;
  those collapsing are position-1 composition echoes (still real, but compositional, not selective).
- tRNA layer: literature annotation only (quantitative supply quarantined, adversarial review B3).
- status: exploratory interpretation aid; no claim change.
""")
mani = pd.DataFrame({"file": [path, gff, f"{L1}/L1_response_by_gene.tsv", f"{L1}/L1_GSE226632_dTE.tsv"]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(os.path.join(OUT, f))]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
