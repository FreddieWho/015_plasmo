#!/usr/bin/env python
# L1: gene-level synonymous-codon x stress interaction (exploratory, D-018/D-019)
# Discovery: GSE151189 (DHA x K13 time course) ; Replication: GSE226632 (AA starvation, polysome/total TE)
# Model: resp ~ feature + gc3 + time (gene-cluster bootstrap B=1000); replication: dTE ~ feature + gc3
import gzip, hashlib, json, math, os, re, tarfile
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
RAW = f"{ROOT}/data/raw"
OUT = f"{ROOT}/data/derived/WP4/L1_codon_x_stress"
os.makedirs(OUT, exist_ok=True)
ACC_PF = "GCF_000002765.6"
SEED = 20260904
rng = np.random.default_rng(SEED)
B = 1000

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STD_CODE = dict(zip(CODONS, AA))
SENSE = [c for c in CODONS if STD_CODE[c] != "*"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s))
    return m.group(0) if m else None

def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(len(p), np.nan)
    if ok.sum() == 0: return q
    po = p[ok]; n = len(po); o = np.argsort(po)
    qq = np.empty(n); qq[o] = np.minimum.accumulate((po[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q[ok] = np.minimum(qq, 1)
    return q

# ---------- 1. Pf CDS codon counts, LCR dual-track ----------
def lcr_mask(prot):
    n = len(prot); masked = np.zeros(n, bool)
    if n < 30: return masked
    arr = np.array(list(prot))
    for i in range(n - 29):
        w = arr[i:i + 30]; cnt = Counter(w)
        ent = -sum((v / 30) * math.log(v / 30, 2) for v in cnt.values())
        if ent < 1.5: masked[i:i + 30] = True
    return masked

def parse_cds(acc):
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
            if "[partial" in hid or len(s) < 150 or len(s) % 3: return
            if not set(s) <= set("ATGC"): return
            cod = [s[i:i + 3] for i in range(0, len(s), 3)]
            if "*" in [STD_CODE.get(c, "*") for c in cod][:-1]: return
            if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
        for line in fh:
            if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
            else: seq.append(line.strip())
        flush()
    rows = {}
    for g, (_, s) in best.items():
        cod = [s[i:i + 3] for i in range(0, len(s), 3) if STD_CODE.get(s[i:i + 3], "*") != "*"]
        if not cod: continue
        prot = "".join(STD_CODE[c] for c in cod)
        mask = lcr_mask(prot)
        cc = Counter(cod); ccm = Counter(c for i, c in enumerate(cod) if not mask[i])
        gc3 = sum(1 for c in cod if c[2] in "GC") / len(cod)
        rows[g] = {"n_cod": len(cod), "n_codm": int((~mask).sum()), "gc3": gc3,
                   "lcr_frac": float(mask.mean()),
                   **{f"cod_{c}": cc.get(c, 0) for c in SENSE},
                   **{f"codm_{c}": ccm.get(c, 0) for c in SENSE}}
    df = pd.DataFrame(rows).T
    df.index.name = "gene"
    return df

CDF = parse_cds(ACC_PF)
print("CDS genes:", len(CDF), flush=True)

def features(df, pre, n):
    """focus features from count columns with prefix pre, denominator column n"""
    d = df[n].replace(0, np.nan)
    f = pd.DataFrame(index=df.index)
    f["AAA_freq"] = df[pre + "AAA"] / d
    f["AAG_freq"] = df[pre + "AAG"] / d
    f["AAA_share"] = df[pre + "AAA"] / (df[pre + "AAA"] + df[pre + "AAG"]).replace(0, np.nan)
    f["ATT_freq"] = df[pre + "ATT"] / d
    f["ATCATA_freq"] = (df[pre + "ATC"] + df[pre + "ATA"]) / d
    f["gc3"] = df["gc3"]
    return f

FEAT_U = features(CDF, "cod_", "n_cod")
FEAT_M = features(CDF, "codm_", "n_codm")
FOCUS = ["AAA_freq", "AAG_freq", "AAA_share", "ATT_freq", "ATCATA_freq"]
EXPL = [f"{c}_freq" for c in SENSE if f"{c}_freq" not in FOCUS]
EXPL_DF = pd.DataFrame({f"{c}_freq": CDF[f"cod_{c}"] / CDF["n_cod"].replace(0, np.nan) for c in SENSE})
EXPL_DF["gc3"] = CDF["gc3"]
print("focus+expl features ready", flush=True)

# ---------- 2. platform probe -> PF3D7 (parsed from RAW soft, provenance fix) ----------
soft = f"{RAW}/geo/GPL18893/GPL18893_family.soft.gz"
p2g = {}
with gzip.open(soft, "rt", errors="replace") as fh:
    in_tab = False; hdr = None; i_id = i_orf = None
    for line in fh:
        if line.startswith("!platform_table_begin"):
            in_tab = True; continue
        if line.startswith("!platform_table_end"):
            break
        if not in_tab: continue
        if hdr is None:
            hdr = line.rstrip("\n").split("\t")
            hdr = [h.strip().strip('"') for h in hdr]
            i_id = hdr.index("ID"); i_orf = hdr.index("ORF")
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(i_id, i_orf): continue
        g = gkey(p[i_orf])
        if g: p2g[p[i_id].strip('"')] = g
print("probes mapped:", len(p2g), flush=True)

# ---------- 3. GSE151189 series matrix -> per-gene log2 ----------
mx_path = f"{RAW}/geo/GSE151189/GSE151189_series_matrix.txt.gz"
titles = []
with gzip.open(mx_path, "rt") as fh:
    for line in fh:
        if line.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            break
df = pd.read_csv(mx_path, sep="\t", comment="!", header=0, index_col=0, low_memory=False)
df.columns = titles
def parse_sample(t):
    dha = "_DHA_" in t
    parts = t.split("_DHA_")[0].split("_") if dha else t.rsplit("_rep", 1)[0].split("_")
    return "_".join(parts[:2]), float(parts[2][:-1]), dha, t.rsplit("_rep", 1)[1]
meta = pd.DataFrame([parse_sample(t) for t in titles], columns=["bg", "time", "dha", "rep"], index=titles)

g2rows = {}
for pid, g in p2g.items():
    if pid in df.index: g2rows.setdefault(g, []).append(pid)
genes = sorted(set(g2rows) & set(CDF.index))
G = pd.DataFrame({g: df.loc[g2rows[g]].mean(axis=0) for g in genes}).T.reindex(columns=titles)
print("gene x sample matrix:", G.shape, flush=True)

# ---------- 4. per-gene DHA response per bg x time ----------
resp_rows = []
for (bg, tm), grp in meta.groupby(["bg", "time"]):
    d = grp[grp.dha].index.tolist(); c = grp[~grp.dha].index.tolist()
    if not d or not c: continue
    dm = G[d].mean(axis=1, skipna=True); cm = G[c].mean(axis=1, skipna=True)
    nn = G[d].notna().sum(axis=1); nc = G[c].notna().sum(axis=1)
    r = pd.DataFrame({"gene": G.index, "bg": bg, "time": tm, "resp": dm - cm,
                      "n_dha": nn.values, "n_ctrl": nc.values,
                      "usable": (nn >= 2) & (nc >= 2)})
    resp_rows.append(r)
RESP = pd.concat(resp_rows)
RESP.to_csv(f"{OUT}/L1_response_by_gene.tsv", sep="\t", index=False)
strata = RESP.groupby(["bg", "time"])["usable"].sum()
print("usable strata:\n", strata.to_string(), flush=True)

# ---------- 5. slope tests with bootstrap ----------
def ols_beta(x, C, y):
    X = np.column_stack([x, C, np.ones(len(x))])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return b[0]

def boot_slope(x, C, y, groups=None, B=B):
    """gene-cluster bootstrap: resample groups (genes), all rows of a sampled gene enter together (weighted lstsq)"""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(C).all(axis=1)
    x, C, y = x[ok], C[ok], y[ok]
    n = len(y)
    if n < 100 or np.std(x) == 0 or np.std(y) == 0: return np.nan, np.nan, np.nan, n
    groups = np.arange(n) if groups is None else np.asarray(groups)[ok]
    u, inv = np.unique(groups, return_inverse=True)
    def wbeta(w):
        sw = np.sqrt(w)
        X = np.column_stack([x * sw, C * sw[:, None], sw])
        b, *_ = np.linalg.lstsq(X, y * sw, rcond=None)
        return b[0]
    b0 = wbeta(np.ones(n))
    bs = np.empty(B)
    for i in range(B):
        counts = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u))
        try: bs[i] = wbeta(counts[inv].astype(float))
        except Exception: bs[i] = np.nan
    bs = bs[np.isfinite(bs)]
    if len(bs) < B * 0.9: return b0, np.nan, np.nan, n
    p = 2 * min((np.sum(bs <= 0) + 1) / (len(bs) + 1), (np.sum(bs >= 0) + 1) / (len(bs) + 1))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return b0, p, (hi - lo) / 2, n

def slope_table(FEATS_DF, feat_list, tag):
    rows = []
    Dd2 = ["Dd2_WT", "Dd2_R539T", "Dd2_C580Y"]
    for f in feat_list:
        fv = FEATS_DF[f]
        slopes = {}
        for bg in Dd2:
            sub = RESP[(RESP.bg == bg) & (RESP.time <= 6) & RESP.usable]
            j = sub.merge(pd.DataFrame({"gene": FEATS_DF.index.values, "f": fv.values, "gc3": FEATS_DF["gc3"].values}), on="gene")
            x = j["f"].values.astype(float)
            C = np.column_stack([j["gc3"].values.astype(float), j["time"].values.astype(float)])
            y = j["resp"].values.astype(float)
            b, p, ci, n = boot_slope(x, C, y, groups=j["gene"].values)
            slopes[bg] = b
            rho = stats.spearmanr(j["f"], j["resp"], nan_policy="omit").statistic if n >= 100 else np.nan
            rows.append({"feature": f, "track": tag, "bg": bg, "slope": b, "ci95_half": ci,
                         "p_boot": p, "n_genes": n, "spearman": rho})
        for mut in ["Dd2_R539T", "Dd2_C580Y"]:
            if np.isfinite(slopes.get(mut, np.nan)) and np.isfinite(slopes.get("Dd2_WT", np.nan)):
                rows.append({"feature": f, "track": tag, "bg": f"{mut}-WT_interaction",
                             "slope": slopes[mut] - slopes["Dd2_WT"], "ci95_half": np.nan,
                             "p_boot": np.nan, "n_genes": np.nan, "spearman": np.nan})
    R = pd.DataFrame(rows)
    m = R.bg.isin(Dd2)
    R.loc[m, "q_bh"] = bh(R.loc[m, "p_boot"].values)
    return R

FO_U = slope_table(FEAT_U, FOCUS, "unmasked")
FO_M = slope_table(FEAT_M, FOCUS, "masked")
FOCUS_ALL = pd.concat([FO_U, FO_M])
FOCUS_ALL.to_csv(f"{OUT}/L1_focus_slopes.tsv", sep="\t", index=False)
EX_U = slope_table(EXPL_DF, EXPL, "exploratory_unmasked")
EX_U.to_csv(f"{OUT}/L1_exploration_slopes.tsv", sep="\t", index=False)
print("slopes done", flush=True)

# ---------- 6. replication: GSE226632 TE under AA starvation ----------
tar = tarfile.open(f"{RAW}/geo/GSE226632/GSE226632_RAW.tar")
HT = {}
for m in tar.getnames():
    base = m.split("/")[-1]
    if "htseq" in base:
        parts = base.replace(".tabular.txt.gz", "").split("_htseq-count-")[1].split("-")
        fh = tar.extractfile(m); d = {}
        with gzip.open(fh, "rt") as z:
            for line in z:
                p = line.split()
                if len(p) >= 2 and p[0].startswith("PF3D7_"): d[p[0]] = float(p[1])
        HT[(parts[0], parts[1], parts[2])] = d
tar.close()
te = {}
for cond in ["CM", "AA_free"]:
    per_rep = []
    for r in ["replicate1", "replicate2", "replicate3"]:
        tot, pol = HT[(cond, "total", r)], HT[(cond, "polysome", r)]
        gs = set(tot) & set(pol)
        per_rep.append({g: math.log2((pol[g] + 1) / (tot[g] + 1)) for g in gs})
    gs = set.intersection(*[set(d) for d in per_rep])
    te[cond] = {g: float(np.mean([d[g] for d in per_rep])) for g in gs}
genes_te = sorted(set(te["CM"]) & set(te["AA_free"]))
DTE = pd.DataFrame({"gene": genes_te,
                    "dTE": [te["AA_free"][g] - te["CM"][g] for g in genes_te]})
DTE.to_csv(f"{OUT}/L1_GSE226632_dTE.tsv", sep="\t", index=False)

# replicate: FULL 61-codon scan + focus extras (D-022: discovery scans unrestricted; interpretation after significance)
ALLF = EXPL_DF.copy()
ALLF["AAA_share"] = FEAT_U["AAA_share"]; ALLF["ATCATA_freq"] = FEAT_U["ATCATA_freq"]
rep_feats = [f"{c}_freq" for c in SENSE] + ["AAA_share", "ATCATA_freq"]
Dd2bgs = ["Dd2_WT", "Dd2_R539T", "Dd2_C580Y"]
def find_hits(T):
    d = T[T.bg.isin(Dd2bgs)]
    out = {}
    for f, s in d.groupby("feature"):
        sig = s[s.q_bh < 0.05]
        if len(sig) >= 2 and pd.Series(np.sign(sig.slope)).nunique() == 1:
            out[f] = float(np.sign(sig.slope.iloc[0]))
    return out
hits = find_hits(pd.concat([FO_U, EX_U]))
rep_rows = []
for f in rep_feats:
    fv = ALLF[f]
    j = DTE.merge(pd.DataFrame({"gene": ALLF.index.values, "f": fv.values, "gc3": ALLF["gc3"].values}), on="gene")
    x = j["f"].values.astype(float); C = np.column_stack([j["gc3"].values.astype(float)])
    y = j["dTE"].values.astype(float)
    b, p, ci, n = boot_slope(x, C, y)
    rho = stats.spearmanr(j["f"], j["dTE"], nan_policy="omit").statistic
    rep_rows.append({"feature": f, "dTE_slope": b, "ci95_half": ci, "p_boot": p, "n_genes": n,
                     "spearman": rho, "discovery_hit": f in hits,
                     "discovery_sign": hits.get(f, np.nan)})
REP = pd.DataFrame(rep_rows)
REP["q_bh"] = bh(REP["p_boot"].values)
REP["replicated"] = REP.apply(lambda r: bool(np.isfinite(r.discovery_sign)
    and np.sign(r.dTE_slope) == r.discovery_sign and r.q_bh < 0.10), axis=1)
REP.to_csv(f"{OUT}/L1_TE_replication.tsv", sep="\t", index=False)
print("replication done; discovery hits:", sorted(hits), flush=True)

# ---------- 7. K13 control row ----------
k13 = "PF3D7_1343700"
k13_rows = [{"feature": f, "value": float(FEAT_U.loc[k13, f]) if k13 in FEAT_U.index else np.nan}
            for f in FOCUS + ["gc3"]]
kr = RESP[RESP.gene == k13][["bg", "time", "resp", "usable"]].to_dict("records")
pd.DataFrame({"k13_features": json.dumps(k13_rows, ensure_ascii=False),
              "k13_responses": json.dumps(kr, ensure_ascii=False)}, index=[0]).to_csv(
    f"{OUT}/L1_K13_row.tsv", sep="\t", index=False)

# ---------- 8. verdict + deliverables ----------
replicated = REP[REP.replicated].feature.tolist()
verdict = ("DOUBLE_POSITIVE" if replicated else
           "DISCOVERY_ONLY" if hits else "DOUBLE_NEGATIVE")
qc = {"genes_cds": int(len(CDF)), "genes_with_probe": int(len(genes)),
      "strata": {f"{a}@{b}": int(v) for (a, b), v in strata.items()},
      "cam3ii_dha_design": "all Cam3II DHA strata are rep1-only (n=1) by experimental design -> excluded from inference, descriptive only",
      "discovery_hits": sorted(hits), "replicated": replicated, "verdict": verdict,
      "bootstrap_B": B, "seed": SEED, "model": "resp ~ feature + gc3 + time; dTE ~ feature + gc3"}
json.dump(qc, open(f"{OUT}/L1_qc.json", "w"), indent=2)
open(f"{OUT}/claim_impact.md", "w").write(f"""# L1 claim impact (exploratory, D-018/D-019 — no preregistration)

- date: 2026-09-04
- verdict: **{verdict}**
- discovery hits (FULL 61-codon scan + focus extras; >=2 Dd2 backgrounds q<0.05 same sign): {sorted(hits)}
- replicated in GSE226632 dTE (same sign, q<0.10): {replicated}
- rule: double-positive -> evidence card, claimable in main text; discovery-only -> supplementary, labeled single-dataset; double-negative -> H3 gene-level no-evidence.
- CLM04 impact: none until user review; Gate C PIVOT untouched.
- caveats: microarray log2ratio vs 3D7 pool; Cam3II DHA rep1-only by design (descriptive only); gc3 covariate included (beyond-composition signal); GENE-CLUSTER bootstrap B={B} (all rows of a gene resampled together); exploratory label per D-019; ATT_freq discovery sign (negative under DHA) FLIPS in starvation dTE (positive) -> cross-perturbation sign-inconsistent, hence discovery-only.
""")
mani = pd.DataFrame({"file": [mx_path, soft,
                              f"{RAW}/ncbi-datasets/{ACC_PF}/ncbi_dataset/data/{ACC_PF}/cds_from_genomic.fna",
                              f"{RAW}/geo/GSE226632/GSE226632_RAW.tar"]})
mani["sha256"] = mani.file.map(sha)
mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(os.path.join(OUT, f))]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", verdict, flush=True)
