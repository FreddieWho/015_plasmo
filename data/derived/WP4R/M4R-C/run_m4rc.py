#!/usr/bin/env python3
"""M4R-C ROUND1: stage-adjusted DHA, acute-vs-chronic partition, Asn feature
decomposition, independent ART validation, dTE interaction upgrade attempt.

Task packet: docs/tasks/M4R-C_TASK_PACKET.md (M4R-C-ROUND1).
Seed: 20260905. Read-only inputs; writes only to data/derived/WP4R/M4R-C/.
Labels: transcript/mRNA contrasts are NEVER translation evidence (M4R-X01 stands).
"""
import gzip, json, math, os, re, zipfile
import numpy as np
import pandas as pd
from scipy import stats

SEED = 20260905
rng = np.random.default_rng(SEED)
OUT = "data/derived/WP4R/M4R-C"
os.makedirs(OUT, exist_ok=True)

def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out

# ---------------- A. protein features ----------------
print("A: protein features...")
recs, hdr, seq = [], None, []
with open("data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta") as h:
    for line in h:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(seq))); seq = []
            hdr = line
        else:
            seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))

def shannon(s):
    from collections import Counter
    c = Counter(s); n = len(s)
    return -sum((v/n)*math.log2(v/n) for v in c.values())

W = 12
feat_rows = []
for hdr, s in recs:
    m = re.search(r"gene=(PF3D7_[A-Za-z0-9]+)", hdr)
    mp = re.search(r"gene_product=([^|]+)", hdr)
    gene = m.group(1); prod = (mp.group(1).strip() if mp else "")
    L = len(s)
    nN = s.count("N"); asn = nN/L if L else np.nan
    maxrun, run, cnt5 = 0, 0, 0
    for aa in s:
        if aa == "N": run += 1
        else:
            if run >= 5: cnt5 += 1
            maxrun = max(maxrun, run); run = 0
    if run >= 5: cnt5 += 1
    maxrun = max(maxrun, run)
    # LCR windows
    lc = np.zeros(L, dtype=bool)
    if L >= W:
        for i in range(L - W + 1):
            if shannon(s[i:i+W]) < 1.5:
                lc[i:i+W] = True
    lcr_frac = lc.mean()
    n_in = sum(1 for i, aa in enumerate(s) if aa == "N" and lc[i])
    asn_in = n_in/L; asn_out = (nN-n_in)/L
    pl = prod.lower()
    feat_rows.append(dict(gene=gene, length=L, asn_frac=asn, polyN_maxrun=maxrun,
        polyN_runs_ge5=cnt5, lcr_frac=lcr_frac, asn_in_lcr=asn_in,
        asn_out_lcr=asn_out, is_AP2=("ap2" in pl),
        is_PUF=("puf" in pl or "pumilio" in pl),
        is_chromatin=("chromatin" in pl or "histone" in pl),
        is_RNAbind=("rna-binding" in pl or "rna binding" in pl),
        product=prod))
FEAT = pd.DataFrame(feat_rows).set_index("gene")
FEAT.to_csv(f"{OUT}/M4RC_protein_features.tsv", sep="\t")
asn90 = FEAT.asn_frac.quantile(0.9)
print(f"  n={len(FEAT)} AsnTop10 thr={asn90:.4f} AP2={FEAT.is_AP2.sum()} PUF={FEAT.is_PUF.sum()} CHROM={FEAT.is_chromatin.sum()} RNAB={FEAT.is_RNAbind.sum()}")

# ---------------- feature-test engine ----------------
def ftests(y, label, feat=FEAT, extra_covar=None):
    """y: pd.Series indexed by gene. Returns list of dict rows."""
    d = pd.DataFrame({"y": y}).join(feat, how="inner")
    d = d[np.isfinite(d.y)]
    # length-residualized y
    x = np.log10(d.length.values); yy = d.y.values
    X = np.column_stack([np.ones(len(d)), x])
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    d = d.assign(y_lenres=yy - X@beta)
    top = d.asn_frac >= d.asn_frac.quantile(0.9)
    rows = []
    num_feats = ["asn_frac", "asn_out_lcr", "asn_in_lcr", "polyN_maxrun", "lcr_frac", "length"]
    for f in num_feats:
        v = d[f].values
        r1 = stats.spearmanr(v, d.y.values, nan_policy="omit")
        r2 = stats.spearmanr(v, d.y_lenres.values, nan_policy="omit")
        a = d.loc[top, "y"].values; b = d.loc[~top, "y"].values
        try: mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        except Exception: mw = None
        rows.append(dict(contrast=label, feature=f, n=len(d),
            spearman_rho=r1.statistic, spearman_p=r1.pvalue,
            lenadj_rho=r2.statistic, lenadj_p=r2.pvalue,
            MW_AsnTop10_p=(mw.pvalue if mw else np.nan),
            median_Top10=np.median(a), median_rest=np.median(b)))
    # regulator flags: MW of y in set vs rest
    for f in ["is_AP2", "is_PUF", "is_chromatin", "is_RNAbind"]:
        a = d.loc[d[f], "y"].values; b = d.loc[~d[f], "y"].values
        if len(a) < 3:
            rows.append(dict(contrast=label, feature=f, n=len(d), spearman_rho=np.nan,
                spearman_p=np.nan, lenadj_rho=np.nan, lenadj_p=np.nan,
                MW_AsnTop10_p=np.nan, median_Top10=np.nan, median_rest=np.nan)); continue
        try: mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        except Exception: mw = None
        rows.append(dict(contrast=label, feature=f, n=len(d), spearman_rho=np.nan,
            spearman_p=np.nan, lenadj_rho=np.nan, lenadj_p=np.nan,
            MW_AsnTop10_p=(mw.pvalue if mw else np.nan),
            median_Top10=np.median(a), median_rest=np.median(b)))
    return rows

ALL = []

def checkpoint(tag):
    global ALL
    if not ALL:
        print(f"  [ckpt {tag}] nothing yet"); return
    RES = pd.DataFrame(ALL)
    pcols = [c for c in RES.columns if c.endswith("_p")]
    for c in pcols:
        RES[c.replace("_p", "_q")] = RES.groupby("contrast")[c].transform(lambda v: bh_fdr(v.fillna(1).values))
    RES.to_csv(f"{OUT}/M4RC_feature_tests_all.tsv", sep="\t", index=False)
    qcol = "spearman_q" if "spearman_q" in RES.columns else "spearman_p"
    piv = RES[RES.feature == "asn_frac"][["contrast", "n", "spearman_rho", qcol]]
    piv.to_csv(f"{OUT}/M4RC_validation_table.tsv", sep="\t", index=False)
    print(f"  [ckpt {tag}] contrasts={RES.contrast.nunique()} rows={len(RES)}")

# ---------------- B. MCA stage reference ----------------
print("B: MCA stage reference...")
z = zipfile.ZipFile("data/raw/mca/pf-ch10x-set1/pf-ch10x-set1.zip")
meta = pd.read_csv(z.open("pf-ch10x-set1-ch10x-data.csv"), index_col=0)
print("  stage_HR:", meta.STAGE_HR.value_counts().to_dict())
exp = pd.read_csv(z.open("pf-ch10x-set1-ch10x-exp.csv"), index_col=0)
exp = exp.loc[:, exp.columns.isin(meta.index)]
meta = meta.loc[exp.columns]
stage_prof = exp.T.groupby(meta.STAGE_HR.values).mean().T  # gene x stage
stage_prof.to_csv(f"{OUT}/M4RC_MCA_stage_pseudobulk.tsv", sep="\t")
stages = list(stage_prof.columns)
ring_c = [c for c in stages if "ring" in c.lower()]
sch_c = [c for c in stages if "schizont" in c.lower()]
stagevec = (stage_prof[sch_c].mean(axis=1) - stage_prof[ring_c].mean(axis=1)) if (ring_c and sch_c) else None
# schizont-peak genes (internal stage-confound control definition)
if len(stages) >= 3:
    peak = stage_prof.idxmax(axis=1)
    schiz_peak = set(peak[peak.str.lower().str.contains("schizont")].index)
else:
    schiz_peak = set()
print(f"  genes={stage_prof.shape[0]} stages={stages} schizont-peak n={len(schiz_peak)}")
with open(f"{OUT}/M4RC_schizont_peak_genes.txt", "w") as h:
    h.write("\n".join(sorted(schiz_peak)) + "\n")

# ---------------- C. stage-adjusted DHA (GSE151189 via L1 resp) ----------------
print("C: stage-adjusted DHA...")
RESP = pd.read_csv("data/derived/WP4/L1_codon_x_stress/L1_response_by_gene.tsv", sep="\t")
sv = stagevec.dropna() if stagevec is not None else None
dha_rows = []
for (bg, tm), j in RESP.groupby(["bg", "time"]):
    y = j.set_index("gene")["resp"].dropna()
    y = y[y.index.isin(FEAT.index)]
    if len(y) < 100: continue
    lab = f"GSE151189_DHA-CTRL_{bg}_{tm}h"
    ALL += ftests(y, lab)
    dha_rows.append(dict(contrast=lab, n=len(y), median_resp=y.median()))
    if sv is not None:
        common = y.index.intersection(sv.index)
        if len(common) > 100:
            yy = y.loc[common].values; xx = sv.loc[common].values
            X = np.column_stack([np.ones(len(common)), xx])
            beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
            resid = pd.Series(yy - X@beta, index=common)
            r2 = (1 - np.var(yy - X@beta)/np.var(yy))
            ALL += ftests(resid, lab + "_STAGEADJ")
            dha_rows.append(dict(contrast=lab + "_STAGEADJ", n=len(common),
                median_resp=resid.median(), stage_R2=r2))
DHA = pd.DataFrame(dha_rows)
DHA.to_csv(f"{OUT}/M4RC_stage_adjusted_DHA.tsv", sep="\t", index=False)
checkpoint("C_DHA")

# ---------------- D1. persistence GSE225340 ----------------
print("D1: GSE225340 persistence...")
with gzip.open("data/raw/geo/GSE225340/GSE225340_all_sample_fpkm_table.txt.gz", "rt") as h:
    fpkm = pd.read_csv(h, sep="\t", index_col=0)
fpkm.index = fpkm.index.str.replace(r"\.\d+$", "", regex=True)
ring_cols = [c for c in fpkm.columns if "Rings" in c]
late_cols = [c for c in fpkm.columns if re.search(r"D([5-9]|1[0-2])$", c)]
early_cols = [c for c in fpkm.columns if re.search(r"D[1-4]$", c)]
print(f"  ring={len(ring_cols)} early={len(early_cols)} late={len(late_cols)}")
def dorm_fc(day_cols):
    lfc = {}
    R = fpkm[ring_cols]
    for g in fpkm.index:
        r = R.loc[g].dropna().values; d = fpkm.loc[g, day_cols].dropna().values
        if len(r) >= 2 and len(d) >= 3:
            lfc[g] = np.log2(np.median(d)+1) - np.log2(np.median(r)+1)
    return pd.Series(lfc)
for nm, cols in [("DORM_LATE_vs_Rings", late_cols), ("DORM_EARLY_vs_Rings", early_cols)]:
    y = dorm_fc(cols); print(f"  {nm} n={len(y)} med={y.median():.3f}")
    ALL += ftests(y, f"GSE225340_{nm}")
checkpoint("D1_persist")

# ---------------- D2. K13 scRNA acute x genotype ----------------
print("D2: K13 scRNA...")
dge_dir = "data/raw/k13scrna/M4R_K13_SCRNA_2026/K13_mt_v_wt_sc_v1/data/raw/dge"
smap = pd.read_csv("data/raw/k13scrna/M4R_K13_SCRNA_2026/K13_mt_v_wt_sc_v1/data/raw/sample_map.csv")
pb = {}
for _, r in smap.iterrows():
    f = f"{dge_dir}/{r['sample']}_gene_exon_dge_n20000.txt.gz"
    try:
        with gzip.open(f, "rt") as h:
            m = pd.read_csv(h, sep="\t", index_col=0)
    except FileNotFoundError:
        print(f"  MISSING {f}"); continue
    m.index = m.index.str.replace(r"\.\d+$", "", regex=True)
    s = m.sum(axis=1); s = s[s >= 10]
    pb[r["sample"]] = s
    print(f"  {r['sample']}: cells={m.shape[1]} genes_ok={len(s)}")
S = pd.DataFrame(pb).fillna(0)
cpm = np.log2(S.div(S.sum(axis=0), axis=1)*1e6 + 1)
def sc_contrast(a, b, lab):
    global ALL
    y = (cpm[a] - cpm[b]).dropna()
    y = y[y.index.isin(FEAT.index)]
    print(f"  {lab} n={len(y)} med={y.median():.3f}")
    ALL += ftests(y, lab)
for t in ["2h", "4h", "6h"]:
    sc_contrast(f"MRA1250_{t}DHA", f"MRA1250_{t}DMSO", f"scRNA_WT_{t}_DHAvsDMSO")
    sc_contrast(f"MRA1251_{t}DHA", f"MRA1251_{t}DMSO", f"scRNA_580Y_{t}_DHAvsDMSO")
    # genotype x treatment interaction (DiD), transcript-only
    y = ((cpm[f"MRA1251_{t}DHA"]-cpm[f"MRA1251_{t}DMSO"]) - (cpm[f"MRA1250_{t}DHA"]-cpm[f"MRA1250_{t}DMSO"])).dropna()
    y = y[y.index.isin(FEAT.index)]
    print(f"  scRNA_DiD_{t} n={len(y)} med={y.median():.3f}")
    ALL += ftests(y, f"scRNA_DiD_580YxDHA_{t}")
checkpoint("D2_scRNA")

# ---------------- E1. Mok chronic (reuse L1h_long) ----------------
print("E1: Mok...")
L1h = pd.read_csv("data/derived/WP4/L1h_mok_genotype_context/L1h_long.tsv", sep="\t")
for ct, j in L1h.groupby("contrast"):
    y = j.set_index("gene")["log2FC"].dropna()
    y = y[y.index.isin(FEAT.index)]
    ALL += ftests(y, f"Mok2021_{ct}")
checkpoint("E1_Mok")

# ---------------- E2. GSE59099 chronic transcript ----------------
print("E2: GSE59099...")
pmap = pd.read_csv("data/derived/WP3/M3-03_dha/platform_GPL18893_probe2gene.tsv", sep="\t")
print("  pmap cols:", list(pmap.columns)[:6], "n=", len(pmap))
# parse series matrix TABLE
import subprocess
tbl = subprocess.run(["zcat", "data/raw/geo/GSE59099/GSE59099_series_matrix.txt.gz"],
                     capture_output=True, text=True).stdout.splitlines()
beg = next(i for i, l in enumerate(tbl) if l.startswith("!series_matrix_table_begin"))
hdr = tbl[beg+1].replace('"', "").split("\t")
dat = [l.replace('"', "").split("\t") for l in tbl[beg+2:] if l and not l.startswith("!series_matrix_table_end")]
G = pd.DataFrame(dat, columns=hdr).set_index("ID_REF")
print(f"  matrix {G.shape}")
# metadata: characteristics
meta_rows = [l.replace('"', "").split("\t") for l in tbl if l.startswith("!Sample_characteristics_ch1")]
chars = {}
for row in meta_rows:
    for v in row[1:]:
        if ":" in v:
            k, vv = v.split(":", 1); chars.setdefault(k.strip(), []).append(vv.strip())
CH = pd.DataFrame(chars)
CH.index = [c.strip('"') for c in hdr[1:]]
hl = pd.to_numeric(CH.get("parasite clearance halflife upon artemisinin treatment (h)"), errors="coerce")
tp = CH.get("timepoint"); geo = CH.get("geographic origin")
pre = (tp == "prior to artemisinin combination therapy (ACT)")
mek = geo.isin(["Pailin, Cambodia", "Mae Sot, Thailand", "Preah Vihear, Cambodia", "Binh Phuoc, Vietnam",
                "Rattanakiri, Cambodia", "Pursat, Cambodia", "Attapeu, Laos", "Shwe Kyin, Myanmar",
                "Sisakhet, Thailand", "Ranong, Thailand", "Khun Han, Thailand"])
R = (hl >= 5) & pre; S = (hl <= 3) & pre
print(f"  pre-ACT n={pre.sum()} R(HL>=5)={R.sum()} S(HL<=3)={S.sum()} MekongR={(R&mek).sum()} MekongS={(S&mek).sum()}")
# probe->gene
pc = pmap.columns
probe_col = pc[0]; gene_col = [c for c in pc if "gene" in c.lower() or "PF3D7" in c.lower()][0]
pg = dict(zip(pmap[probe_col].astype(str), pmap[gene_col].astype(str)))
Gnum = G.apply(pd.to_numeric, errors="coerce")
def gene_level(cols):
    sub = Gnum[cols]
    by = {}
    for pr in sub.index:
        g = pg.get(str(pr), "")
        if isinstance(g, str) and g.startswith("PF3D7"):
            by.setdefault(g.split(".")[0], []).append(pr)
    rows = {g: sub.loc[prs].median(axis=0).values for g, prs in by.items() if len(prs)}
    return pd.DataFrame(rows, index=cols)
def g99_test(colsR, colsS, lab):
    GR = gene_level(colsR); GS = gene_level(colsS)
    common = GR.columns.intersection(GS.columns).intersection(FEAT.index)
    y = pd.Series(GR[common].median(axis=0) - GS[common].median(axis=0), index=common)
    print(f"  {lab} genes={len(y)} med={y.median():.4f}")
    ALL += ftests(y, lab)
    # continuous HL association
    cols = [c for c in G.columns if pre.loc[c] and pd.notna(hl.loc[c])]
    GG = gene_level(cols)
    common2 = GG.columns.intersection(FEAT.index)
    hv = hl.loc[cols].values.astype(float)
    rhos = {}
    for g in common2:
        v = GG[g].values
        if np.isfinite(v).sum() > 50:
            rhos[g] = stats.spearmanr(hv[np.isfinite(v)], v[np.isfinite(v)], nan_policy="omit").statistic
    yc = pd.Series(rhos)
    print(f"  {lab}_HLspearman genes={len(yc)} med={yc.median():.4f}")
    ALL += ftests(yc, lab + "_HLspearman")
Rcols = [c for c in G.columns if R.loc[c]]; Scols = [c for c in G.columns if S.loc[c]]
g99_test(Rcols, Scols, "GSE59099_RvsS_preACT")
g99_test([c for c in Rcols if mek.loc[c]], [c for c in Scols if mek.loc[c]], "GSE59099_RvsS_Mekong")
checkpoint("E2_GSE59099")

# ---------------- F. dTE interaction upgrade ----------------
print("F: dTE interaction...")
os.makedirs("/tmp/m4rc/dte", exist_ok=True)
import tarfile
tar = tarfile.open("data/raw/geo/GSE226632/GSE226632_RAW.tar")
tar.extractall("/tmp/m4rc/dte")
conds = {"CM": [], "AA_free": []}
count_files = {"CM": {}, "AA_free": {}}
for rep in [1, 2, 3]:
    for cond, tag in [("CM", "CM"), ("AA_free", "AA_free")]:
        for frac in ["total", "polysome"]:
            fn = f"GSM7080{531+ (rep-1)*4 + (0 if cond=='CM' else 2) + (0 if frac=='total' else 1)}_htseq-count-{tag}-{frac}-replicate{rep}.tabular.txt.gz"
import glob
fmap = {}
for f in glob.glob("/tmp/m4rc/dte/*.tabular.txt.gz"):
    b = os.path.basename(f)
    m = re.search(r"htseq-count-(CM|AA_free)-(total|polysome)-replicate(\d)", b)
    if m: fmap[(m.group(1), m.group(2), int(m.group(3)))] = f
print(f"  files={len(fmap)}")
CNT = {}
for k, f in fmap.items():
    d = pd.read_csv(f, sep="\t", header=None, index_col=0, names=["count"])
    d.index = d.index.str.replace(r"\.\d+$", "", regex=True)
    CNT[k] = d["count"]
C = pd.DataFrame(CNT).fillna(0)
C = C.loc[:, sorted(C.columns)]
# median-ratio size factors
gm = np.exp(np.log(C.where(C > 0).replace(0, np.nan)).mean(axis=1, skipna=True))
sf = (C.div(gm, axis=0)).median(axis=0)
Y = np.log2(C.div(sf, axis=1) + 1)
# design: cond(0/1) frac(0/1) rep dummies + interaction
cols = list(Y.columns)
X = pd.DataFrame({"icept": 1}, index=cols)
X["cond"] = [1 if c[0] == "AA_free" else 0 for c in cols]
X["frac"] = [1 if c[1] == "polysome" else 0 for c in cols]
X["CxF"] = X.cond * X.frac
for r in [2, 3]:
    X[f"rep{r}"] = [1 if c[2] == r else 0 for c in cols]
Xm = X.values; XtXinv = np.linalg.inv(Xm.T @ Xm)
df_res = len(cols) - Xm.shape[1]
inter, tstat, pse = {}, {}, {}
for g in Y.index:
    yv = Y.loc[g].values
    if yv.max() < 1: continue
    beta, *_ = np.linalg.lstsq(Xm, yv, rcond=None)
    resid = yv - Xm@beta
    s2 = (resid@resid)/df_res
    se = math.sqrt(s2*XtXinv[3, 3])
    inter[g] = beta[3]; pse[g] = beta[3]/se if se > 0 else 0.0
I = pd.Series(inter); T = pd.Series(pse)
p_int = pd.Series({g: 2*stats.t.sf(abs(t), df_res) for g, t in T.items()})
q_int = pd.Series(bh_fdr(p_int.values), index=p_int.index)
DTE = pd.DataFrame({"dTE_interaction": I, "t": T, "p": p_int, "q": q_int})
# compare with ratio-based dTE
RATIO = pd.read_csv("data/derived/WP4/L1_codon_x_stress/L1_GSE226632_dTE.tsv", sep="\t").set_index("gene")["dTE"]
common = DTE.index.intersection(RATIO.index)
agree = np.sign(DTE.loc[common, "dTE_interaction"]) == np.sign(RATIO.loc[common])
sp = stats.spearmanr(DTE.loc[common, "dTE_interaction"], RATIO.loc[common], nan_policy="omit")
DTE.to_csv(f"{OUT}/M4RC_dTE_interaction.tsv", sep="\t")
print(f"  genes={len(DTE)} common={len(common)} sign_agree={agree.mean():.3f} spearman={sp.statistic:.3f} (p={sp.pvalue:.2g})")
print(f"  q<0.05 n={(q_int<0.05).sum()}")
# Asn-feature test on interaction dTE
ALL += ftests(DTE.dTE_interaction, "GSE226632_dTEinteraction")
dte_verdict = dict(sign_agreement=float(agree.mean()), spearman=float(sp.statistic),
                   spearman_p=float(sp.pvalue), n_q05=int((q_int < 0.05).sum()))
with open(f"{OUT}/M4RC_dTE_upgrade.json", "w") as h:
    json.dump(dte_verdict, h, indent=2)

# ---------------- finalize ----------------
checkpoint("FINAL")
RES = pd.DataFrame(ALL)
RES.to_pickle(f"{OUT}/M4RC_feature_tests_all.pkl")
print("DONE", RES.shape)
