#!/usr/bin/env python3
"""M3-01 expression-weighted demand D(c,s,t) + translation-output tests.
Task M3-01 | Aim3/WP3/M3/GATE_C tests H3. Seed 42 throughout.
Reads raw only; writes data/derived/WP3/M3-01_demand/ (+ data/interim/M3-01/ rebuildable).
"""
import os, re, json, gzip, tarfile, zipfile, hashlib, io
import numpy as np
import pandas as pd
from scipy import stats as sstats

SEED = 42
rng = np.random.default_rng(SEED)
ROOT = "/home/huyudi/015_plasmo"
RAW = f"{ROOT}/data/raw"
OUT = f"{ROOT}/data/derived/WP3/M3-01_demand"
INTER = f"{ROOT}/data/interim/M3-01"
os.makedirs(OUT, exist_ok=True); os.makedirs(INTER, exist_ok=True)

STD_CODE = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
"TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
"ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
"AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
"TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
CODONS = [c for c in STD_CODE if STD_CODE[c] != "*"]
AA_FAMS = sorted(set(STD_CODE[c] for c in CODONS))
SPLIT = {"L":[("CTN",["CTT","CTC","CTA","CTG"]),("TTR",["TTA","TTG"])],
         "S":[("TCN",["TCT","TCC","TCA","TCG"]),("AGY",["AGT","AGC"])],
         "R":[("CGN",["CGT","CGC","CGA","CGG"]),("AGR",["AGA","AGG"])]}
CAND_CODON_FAMS = ["F","I","P","A","N","D","E","L"]
CAND_AAS = ["A","R","N","D","Q","G","I","K","F","P","W","Y","V"]
ACC = {"SP001":"GCF_000002765.6","SP005":"GCF_000006355.2","SP011":"GCF_900002375.2"}

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# ---------- 1. CDS parse ----------
def parse_cds(acc):
    path = f"{RAW}/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna"
    best = {}
    with open(path) as fh:
        hid, seq = None, []
        def flush():
            if hid is None: return
            m = re.search(r"\[locus_tag=([^\]]+)\]", hid)
            if not m: return
            lt = m.group(1); s = "".join(seq).upper().replace("U","T")
            if "[partial" in hid: return
            if len(s) < 150 or len(s) % 3: return
            if not set(s) <= set("ATGC"): return
            cod = [s[i:i+3] for i in range(0, len(s), 3)]
            if "*" in [STD_CODE.get(c,"*") for c in cod][:-1]: return
            if len(s) > len(best.get(lt, ("", ""))[1]): best[lt] = (hid, s)
        for line in fh:
            if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
            else: seq.append(line.strip())
        flush()
    return best

def lcr_mask_positions(prot):
    # 30aa sliding window Shannon entropy <1.5 -> masked aa positions
    from collections import Counter
    import math
    n = len(prot); masked = np.zeros(n, bool)
    if n < 30: return masked
    arr = np.array(list(prot))
    for i in range(n - 29):
        w = arr[i:i+30]
        counts = Counter(w)
        ent = -sum((v/30)*math.log(v/30, 2) for v in counts.values())
        if ent < 1.5: masked[i:i+30] = True
    return masked

CDS = {}
for sp, acc in ACC.items():
    genes = parse_cds(acc)
    rows = {}
    for lt, (_, s) in genes.items():
        cod = [s[i:i+3] for i in range(0, len(s), 3) if s[i:i+3] in STD_CODE and STD_CODE[s[i:i+3]] != "*"]
        prot = "".join(STD_CODE[c] for c in cod)
        mask = lcr_mask_positions(prot) if len(prot) >= 30 else np.zeros(len(prot), bool)
        cc = {c: 0 for c in CODONS}; ccm = {c: 0 for c in CODONS}
        gc3n = gc3d = 0
        for i, c in enumerate(cod):
            cc[c] += 1
            if not mask[i]: ccm[c] += 1
            if c[2] in "GC": gc3n += 1
            gc3d += 1
        if not cod: continue
        rows[lt] = {"n_cod": len(cod), "gc3": gc3n/max(gc3d,1),
                    "lcr_frac": mask.mean(), **{f"cod_{c}":cc[c] for c in CODONS},
                    **{f"codm_{c}":ccm[c] for c in CODONS}}
    CDS[sp] = pd.DataFrame(rows).T
    CDS[sp].index.name = "gene"
print("CDS genes:", {k: len(v) for k,v in CDS.items()}, flush=True)

# ---------- 2. M2 residuals + genome GC ----------
m2cod = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-01_codon/M2-01_codon_residual_table.tsv", sep="\t")
m2aa = pd.read_csv(f"{ROOT}/data/derived/WP2/M2-02_aa/M2-02_aa_residual_table.tsv", sep="\t")
m1comp = pd.read_csv(f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv", sep="\t")
gcol = [c for c in m1comp.columns if "genome" in c.lower() and "gc" in c.lower()][0]
GC = dict(zip(m1comp.iloc[:,0], m1comp[gcol]))
# species_id column name?
spcol = m1comp.columns[0]
GC = dict(zip(m1comp[spcol], m1comp[gcol]))
print("GC:", {s: GC.get(s) for s in ACC}, flush=True)
fam_resid = m2cod[m2cod.stratum=="full"].groupby(["species_id","aa"])["log2_obs_exp_masked"].mean().reset_index()
aa_resid = m2aa[["species_id","aa","residual_masked","cohens_d_masked","genome_gc"]].copy()

# ---------- 3. GEO ----------
geodir = f"{RAW}/geo/GSE226632"
tar = tarfile.open(f"{geodir}/GSE226632_RAW.tar")
members = tar.getnames()
def read_gz_tar(name):
    fh = tar.extractfile(name)
    d = {}
    with gzip.open(fh, "rt") as z:
        for line in z:
            p = line.split()
            if len(p) >= 2 and p[0].startswith("PF3D7_"): d[p[0]] = float(p[1])
    return d
HT, MRNA, TRNA, DEC = {}, {}, {}, {}
for m in members:
    base = m.split("/")[-1]
    if base.startswith("GSM70805") and "htseq" in base:
        parts = base.replace(".tabular.txt.gz","").split("_htseq-count-")[1].split("-")
        HT[(parts[0], parts[1], parts[2].split(".")[0])] = read_gz_tar(m)  # (cond, assay, rep)
    elif re.match(r"GSM70805\d+_[RTS]\dm_count", base):
        cp = re.search(r"_([RTS])(\d)m_count", base)
        MRNA[(cp.group(1), cp.group(2))] = read_gz_tar(m)
    elif "tRNA" in base:
        m0 = re.search(r"_tRNA-(control|AA|HF)_(\d)_-count", base)
        m1 = re.search(r"_tRNA_([RTS])(\d)_count", base)
        m2 = re.search(r"_tRNA-([RTS])(\d)-periodate-count", base)
        if m2: TRNA.setdefault("CHG_%s%s" % (m2.group(1), m2.group(2)), []).append(read_gz_tar(m))
        elif m1: TRNA.setdefault("%s%s" % (m1.group(1), m1.group(2)), []).append(read_gz_tar(m))
        elif m0: TRNA.setdefault("%s_%s" % (m0.group(1), m0.group(2)), []).append(read_gz_tar(m))
    elif re.match(r"GSM8559\d+_R\d[CA]\d+", base):
        cp = re.search(r"_(R\d)([CA])(\d+)", base)
        DEC.setdefault((cp.group(1), cp.group(2)), {})[int(cp.group(3))] = read_gz_tar(m)
print(f"HT keys={len(HT)} mRNA={len(MRNA)} tRNAconds={list(TRNA)} decay={len(DEC)}", flush=True)
tar.close()

def repmean(dicts):
    genes = set.intersection(*[set(d) for d in dicts])
    return {g: float(np.mean([d[g] for d in dicts])) for g in genes}

HTm = {(c,a): repmean([HT[(c,a,r)] for r in ["replicate1","replicate2","replicate3"]])
       for c in ["CM","AA_free"] for a in ["total","polysome"]}
MRNAm = {st: repmean([MRNA[(st,r)] for r in ["1","2","3"]]) for st in ["R","T","S"]}
TRNAm = {}
for cond, dicts in TRNA.items():
    feats = set.intersection(*[set(d) for d in dicts])
    TRNAm[cond] = {f: float(np.mean([d[f] for d in dicts])) for f in feats}
# aggregate reps: condition supply (control/AA/HF), stage supply (R/T/S), charging (CHG_R/T/S vs total)
def agg_mean(keys):
    dd = [TRNAm[k] for k in keys if k in TRNAm]
    feats = set.intersection(*[set(d) for d in dd])
    return {f: float(np.mean([d[f] for d in dd])) for f in feats}
SUP = {c: agg_mean([f"{c}_{i}" for i in ["1","2","3"]]) for c in ["control","AA","HF"]}
SUP.update({s: agg_mean([f"{s}{i}" for i in ["1","2","3"]]) for s in ["R","T","S"]})

# ---------- 4. MCA pseudo-bulk ----------
MCA = {}
for tag, zf, ef, df in [("SP001","pf-ch10x-set1","pf-ch10x-set1-ch10x-exp.csv","pf-ch10x-set1-ch10x-data.csv"),
                        ("SP005","pk-ch10x-set1","pk-ch10x-set1-ch10x-exp.csv","pk-ch10x-set1-ch10x-data.csv"),
                        ("SP011","pb-ss2-set1","pb-ss2-set1-ss2-exp.csv","pb-ss2-set1-ss2-data.csv")]:
    zp = f"{RAW}/mca/{zf}/{zf}.zip"
    exdir = f"{INTER}/mca_{zf}"; os.makedirs(exdir, exist_ok=True)
    if not (os.path.exists(f"{exdir}/{ef}") and os.path.exists(f"{exdir}/{df}")):
        with zipfile.ZipFile(zp) as z: z.extractall(exdir)
    meta = pd.read_csv(f"{exdir}/{df}")
    stagecol = "STAGE_HR" if "STAGE_HR" in meta.columns else "STAGE_LR"
    cellcol = "CELL_ID" if "CELL_ID" in meta.columns else meta.columns[0]
    try:
        exp = pd.read_csv(f"{exdir}/{ef}", index_col=0, engine="c", on_bad_lines="skip", low_memory=False)
        exp = exp.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.float32)
    except Exception:
        exp = pd.read_csv(f"{exdir}/{ef}", index_col=0, engine="python")
        exp = exp.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.float32)
    exp = exp.loc[:, exp.columns.isin(set(meta[cellcol]))]
    stages = {}
    for st, grp in meta.groupby(stagecol):
        cells = [c for c in grp[cellcol] if c in exp.columns]
        if len(cells) >= 10: stages[str(st)] = exp[cells].sum(axis=1)
    MCA[tag] = {"stages": stages, "meta_n": len(meta),
                "ncells": {k: int((meta[stagecol]==k).sum()) for k in stages}}
    print(tag, "stages:", {k: v for k,v in MCA[tag]["ncells"].items()}, flush=True)

# ---------- helpers ----------
def demand_shares(counts_df, weights, masked=False):
    """weights: Series gene->E. returns codon shares + AA shares."""
    pre = "codm_" if masked else "cod_"
    cols = [pre+c for c in CODONS]
    w = weights.reindex(counts_df.index).fillna(0).clip(lower=0)
    tot = (counts_df[cols].multiply(w, axis=0)).sum()
    tot.index = [c.replace(pre, "", 1) for c in tot.index]
    dc = tot / tot.sum()
    aa = {}
    for a in AA_FAMS:
        aa[a] = sum(dc[pre+c] if False else dc[c] for c in CODONS if STD_CODE[c]==a)
    return dc, pd.Series(aa)

def gc3_expected(counts_df, weights, families):
    """per-gene GC3 null expected counts under weights; returns obs, exp per codon."""
    pre = "cod_"  # expression-weighted null uses unmasked counts w/ GC3 rule (matches M2 primary)
    obs = (counts_df[[pre+c for c in CODONS]].multiply(weights.reindex(counts_df.index).fillna(0).clip(lower=0), axis=0)).sum()
    obs.index = [c.replace("cod_","") for c in obs.index]
    exp = {}
    for fam, cods in families.items():
        gc3 = counts_df["gc3"]; n = counts_df[[pre+c for c in cods]].sum(axis=1)
        w = weights.reindex(counts_df.index).fillna(0).clip(lower=0)
        probs = {}
        for c in cods:
            probs[c] = np.where(c[2] in "GC", gc3.values, 1-gc3.values)
        psum = sum(probs.values())
        for c in cods:
            exp[c] = float(((n.values*w.values*probs[c]/np.where(psum==0,1,psum))).sum())
    exp = pd.Series(exp)
    return obs, exp

FAMS = {a: [c for c in CODONS if STD_CODE[c]==a] for a in AA_FAMS}
FAMS_SPLIT = dict(FAMS)
for a, subs in SPLIT.items():
    del FAMS_SPLIT[a]
    for name, cods in subs: FAMS_SPLIT["%s_%s"%(a,name)] = cods

def prop_explained(obs, exp):
    uni = obs.sum()/len(obs)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi_e = ((obs-exp)**2/np.where(exp==0,1,exp)).sum()
        chi_u = ((obs-uni)**2/uni).sum()
    if not np.isfinite(chi_u) or chi_u == 0:
        return float("nan"), float(chi_e), float(chi_u), float(obs.sum())
    return float(1-chi_e/chi_u), float(chi_e), float(chi_u), float(obs.sum())

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]
    return np.minimum(q,1)

def ols_t(X, y):
    X = np.column_stack([np.ones(len(X)), np.asarray(X,float)])
    y = np.asarray(y,float)
    beta, res, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y)-X.shape[1]
    resid = y - X@beta
    s2 = (resid@resid)/max(dof,1)
    cov = s2*np.linalg.pinv(X.T@X)
    se = np.sqrt(np.diag(cov)); t = beta/se
    p = 2*sstats.t.sf(np.abs(t), dof)
    return beta, se, t, p, dof

# ---------- 5. Pf GEO analyses ----------
SP = "SP001"
cdf = CDS[SP]
lcr_flag = cdf["lcr_frac"] > 0.15
res_gene = []
universe = set(cdf.index)
# TE contexts
te = {}
for cond in ["CM","AA_free"]:
    T = pd.Series(HTm[(cond,"total")]); P = pd.Series(HTm[(cond,"polysome")])
    genes = sorted(set(T.index) & set(P.index) & universe)
    te[cond] = pd.DataFrame({"total":T.reindex(genes),"poly":P.reindex(genes)}, index=genes)
    te[cond]["TE"] = np.log2((te[cond]["poly"]+1)/(te[cond]["total"]+1))
te["dTE"] = te["AA_free"]["TE"] - te["CM"]["TE"]
# OLS traits
traits = {"TE_CM": te["CM"]["TE"], "TE_AAfree": te["AA_free"]["TE"], "dTE_AAminusCM": te["dTE"]}
# decay slopes
if DEC:
    for cond in ["C","A"]:
        reps = sorted({k[0] for k in DEC if k[1]==cond})
        slopes = {}
        times = None
        for rep in reps:
            ts = sorted(DEC[(rep,cond)].keys())
            M = pd.DataFrame({t: pd.Series(DEC[(rep,cond)][t]) for t in ts})
            M = np.log2(M+1)
            genes = sorted(set(M.index) & universe)
            M = M.reindex(genes)
            X = np.array(ts, float)
            Xc = X - X.mean()
            sl = (M.values @ Xc) / (Xc @ Xc)
            slopes[rep] = pd.Series(sl, index=genes)
        traits[f"decay_{cond}"] = pd.concat(slopes, axis=1).mean(axis=1)
feat_cod, feat_aa = {}, {}
for c in CODONS:
    v = cdf[f"cod_{c}"]/cdf["n_cod"]
    feat_cod[c] = v
for a in AA_FAMS:
    v = cdf[[f"cod_{c}" for c in CODONS if STD_CODE[c]==a]].sum(axis=1)/cdf["n_cod"]
    feat_aa[a] = v
covars = pd.DataFrame({"loglen": np.log10(cdf["n_cod"]), "gc3": cdf["gc3"], "lcr": cdf["lcr_frac"]})
te_rows = []
for tname, y in traits.items():
    genes = sorted(set(y.dropna().index) & universe)
    yv = y.reindex(genes)
    Xc = covars.reindex(genes)
    for kind, feats in [("codon", feat_cod), ("AA", feat_aa)]:
        ps, bs = [], []
        for name, f in feats.items():
            X = pd.concat([f.reindex(genes), Xc], axis=1).values
            beta, se, t, p, dof = ols_t(X, yv.values)
            te_rows.append({"trait":tname,"kind":kind,"feature":name,"beta":beta[1],"se":se[1],
                            "p":p[1],"n":len(genes),"dof":dof})
te_df = pd.DataFrame(te_rows)
for tname in te_df.trait.unique():
    m = te_df.trait==tname
    te_df.loc[m,"q"] = bh(te_df.loc[m,"p"].values)
te_df.to_csv(f"{OUT}/M3-01_TE_genemodel.tsv", sep="\t", index=False)

# expression-weighted null improvement (Pf GEO contexts + MCA stages later)
imp_rows = []
m2prop = {(r.species_id,r.aa): r.prop_explained_GC3_masked for r in
          pd.read_csv(f"{ROOT}/data/derived/WP2/M2-01_codon/M2-01_background_explained.tsv", sep="\t").itertuples()}
def improvement(sp, ctxname, weights):
    cdfx = CDS[sp]
    genes = sorted(set(weights.index) & set(cdfx.index))
    w = weights.reindex(genes)
    obs, exp = gc3_expected(cdfx.reindex(genes), w, FAMS_SPLIT)
    for fam in FAMS_SPLIT:
        o = obs[FAMS_SPLIT[fam]]; e = exp[FAMS_SPLIT[fam]]
        pe, _, _, N = prop_explained(o, e)
        key = (sp, fam if len(fam)==1 else fam.split("_")[0])
        imp_rows.append({"species":sp,"context":ctxname,"family":fam,"N":N,
                         "prop_expr_weighted":pe,
                         "prop_M2_uniform": m2prop.get(key, np.nan)})
for cond, assay in [("CM","total"),("CM","polysome"),("AA_free","total"),("AA_free","polysome")]:
    improvement("SP001", f"GEO_{cond}_{assay}", pd.Series(HTm[(cond,assay)]))
for st in ["R","T","S"]:
    improvement("SP001", f"GEO_mRNA_{st}", pd.Series(MRNAm[st]))

# ---------- 6. MCA demand + improvement ----------
stage_dem = []
for sp in ["SP001","SP005","SP011"]:
    cdfx = CDS[sp]
    for st, vec in MCA[sp]["stages"].items():
        v = vec.copy(); v.index = [str(i).split(".")[0] for i in v.index]
        v = v.groupby(level=0).sum()  # MCA matrices can repeat gene rows
        genes = sorted(set(v.index) & set(cdfx.index))
        cpm = v.reindex(genes)/v.reindex(genes).sum()*1e6
        dc, da = demand_shares(cdfx.reindex(genes), cpm)
        dc.index = ["codon_"+c for c in dc.index]; da.index = ["AA_"+a for a in da.index]
        row = {"species":sp,"context":f"MCA_{st}","ncells":MCA[sp]["ncells"][st],"ngenes":len(genes)}
        row.update(dc.to_dict()); row.update(da.to_dict())
        stage_dem.append(row)
        improvement(sp, f"MCA_{st}", cpm)
stage_df = pd.DataFrame(stage_dem)
stage_df.to_csv(f"{OUT}/M3-01_stage_demand.tsv", sep="\t", index=False)
imp_df = pd.DataFrame(imp_rows)
imp_df["delta"] = imp_df.prop_expr_weighted - imp_df.prop_M2_uniform
imp_df.to_csv(f"{OUT}/M3-01_expression_weighted_null.tsv", sep="\t", index=False)
# demand shares for GEO contexts too
dem_rows = []
for sp, ctxs in [("SP001", {f"GEO_{c}_{a}": pd.Series(HTm[(c,a)]) for c in ["CM","AA_free"] for a in ["total","polysome"]} | {f"GEO_mRNA_{s}": pd.Series(MRNAm[s]) for s in ["R","T","S"]})]:
    cdfx = CDS[sp]
    for ctx, w in ctxs.items():
        genes = sorted(set(w.index) & set(cdfx.index))
        tot = w.reindex(genes).sum()
        cpm = w.reindex(genes)/tot*1e6
        dc, da = demand_shares(cdfx.reindex(genes), cpm)
        row = {"species":sp,"context":ctx,"ngenes":len(genes)}
        row.update({"codon_"+c: dc[c] for c in dc.index})
        row.update({"AA_"+a: da[a] for a in da.index})
        dem_rows.append(row)
pd.DataFrame(dem_rows).to_csv(f"{OUT}/M3-01_GEO_demand.tsv", sep="\t", index=False)

# ---------- 7. tRNA supply + charging ----------
# TRNAm keys are per-file tags (control_1..3, AA_1..3, HF_1..3, R1..S3, CHG_R1..);
# SUP (built above) holds rep-aggregated supply dicts keyed control/AA/HF/R/T/S.
trna_rows = []
anti = {}
for cond, d in SUP.items():
    s = pd.Series(d)
    m = s.index.to_series().str.extract(r"HS_([A-Za-z]+)-([ACGTU]+)-")
    df = pd.DataFrame({"aa3":m[0],"ac":m[1],"count":s.values}).dropna()
    by_ac = df.groupby("ac")["count"].sum()
    trna_rows.append({"context":f"tRNA_{cond}","n_features":len(s),"total":float(s.sum()),
                       "n_anticodons":int(len(by_ac))})
    anti[cond] = by_ac
# charging fraction per stage: periodate(charged)/total, anticodon level
chg_rows = []
for s in ["R","T","S"]:
    tot = anti.get(s)
    chg_raw = agg_mean([f"CHG_{s}{i}" for i in ["1","2","3"]])
    if tot is not None and chg_raw:
        cs = pd.Series(chg_raw)
        m = cs.index.to_series().str.extract(r"HS_([A-Za-z]+)-([ACGTU]+)-")
        cdf2 = pd.DataFrame({"aa3":m[0],"ac":m[1],"count":cs.values}).dropna()
        chg = cdf2.groupby("ac")["count"].sum()
        frac = (chg/tot.reindex(chg.index)).fillna(0)
        for ac, v in frac.items():
            chg_rows.append({"stage":s,"anticodon":ac,"charged_frac":float(v),
                             "charged":float(chg[ac]),"total":float(tot[ac])})
chg_df = pd.DataFrame(chg_rows)
chg_df.to_csv(f"{OUT}/M3-01_tRNA_charging.tsv", sep="\t", index=False)
pd.DataFrame(trna_rows).to_csv(f"{OUT}/M3-01_tRNA_qc.tsv", sep="\t", index=False)
# supply per codon via WC complement + G:U wobble at codon pos1
comp = {"A":"U","U":"A","G":"C","C":"G"}
def pairs(codon):
    codon = codon.replace("T", "U")  # DNA->RNA; anticodon keys are RNA letters
    out = []
    cac = "".join(comp[b] for b in codon[::-1])  # WC anticodon
    out.append(cac)
    wob = {"U":["G"],"G":["U"],"C":["U"],"A":["U"]}  # codon-pos1 alternatives
    for alt in wob.get(codon[0], []):
        out.append("".join(comp[b] for b in (alt+codon[1:])[::-1]))
    return out
sup_rows = []
for cond, by_ac in anti.items():
    for c in CODONS:
        sup_rows.append({"context":f"tRNA_{cond}","codon":c,"AA":STD_CODE[c],
                         "supply":float(sum(by_ac.get(a,0) for a in pairs(c)))})
sup_df = pd.DataFrame(sup_rows)
# S/D per AA: supply(AA)/demand(AA) using GEO CM-total demand
dem = pd.read_csv(f"{OUT}/M3-01_GEO_demand.tsv", sep="\t")
dcm = dem[dem.context=="GEO_CM_total"].iloc[0]
for cond in anti:
    tot_s = sup_df[sup_df.context==f"tRNA_{cond}"].groupby("AA").supply.sum()
    tot_s = tot_s/tot_s.sum()
    for a in AA_FAMS:
        d = dcm.get("AA_"+a, np.nan)
        sup_df.loc[sup_df.context==f"tRNA_{cond}", "SD_"+a] = float(tot_s.get(a,0))/d if d else np.nan
sup_df.to_csv(f"{OUT}/M3-01_tRNA_supply.tsv", sep="\t", index=False)

# ---------- 8. primary contrast Pf-vs-Pk ----------
fam_mean = fam_resid  # species x aa signed mean log2obs/exp
def fam_val(sp, fam):
    if "_" in fam:
        a, sub = fam.split("_",1)
        sub = dict((n,c) for n,c in SPLIT[a])[sub]
        rows = m2cod[(m2cod.species_id==sp)&(m2cod.aa==a)&(m2cod.codon.isin(sub))&(m2cod.stratum=="full")]
    else:
        rows = m2cod[(m2cod.species_id==sp)&(m2cod.aa==fam)&(m2cod.stratum=="full")]
    return float(rows.log2_obs_exp_masked.mean())
rows = []
dem_all = pd.concat([stage_df[["species","context"]+ [c for c in stage_df.columns if c.startswith("codon_") or c.startswith("AA_")]]], ignore_index=True)
for fam in list(FAMS_SPLIT.keys()):
    try: rv = (fam_val("SP001",fam), fam_val("SP005",fam))
    except Exception: continue
    # demand: family share mean over MCA contexts per species
    def dshare(sp):
        sub = dem_all[dem_all.species==sp]
        if fam in FAMS_SPLIT:
            cols = ["codon_"+c for c in FAMS_SPLIT[fam]]
        else: cols = []
        return float(sub[cols].sum(axis=1).mean()) if cols else np.nan
    dv = (dshare("SP001"), dshare("SP005"))
    rows.append({"layer":"codon","family":fam,"resid_Pf":rv[0],"resid_Pk":rv[1],
                 "dresid":rv[0]-rv[1],"dem_Pf":dv[0],"dem_Pk":dv[1],"ddem":dv[0]-dv[1],
                 "match": bool(np.sign(rv[0]-rv[1])==np.sign(dv[0]-dv[1]))})
for a in CAND_AAS:
    r = aa_resid[aa_resid.aa==a].set_index("species_id").residual_masked
    def dshareA(sp):
        sub = dem_all[dem_all.species==sp]
        return float(sub["AA_"+a].mean())
    dv = (dshareA("SP001"), dshareA("SP005"))
    dr = float(r.get("SP001",np.nan)-r.get("SP005",np.nan))
    rows.append({"layer":"AA","family":a,"resid_Pf":float(r.get("SP001",np.nan)),"resid_Pk":float(r.get("SP005",np.nan)),
                 "dresid":dr,"dem_Pf":dv[0],"dem_Pk":dv[1],"ddem":dv[0]-dv[1],
                 "match": bool(np.sign(dr)==np.sign(dv[0]-dv[1]))})
pc = pd.DataFrame(rows)
pc.to_csv(f"{OUT}/M3-01_primary_contrast_Pf_vs_Pk.tsv", sep="\t", index=False)
k = int(pc.match.sum()); n = len(pc)
ci = (sstats.beta.ppf(0.025,k+1,n-k+1) if n>k else 0.0, sstats.beta.ppf(0.975,k+1,n-k) if k>0 else 1.0)
pv = float(sstats.binomtest(k,n,0.5,alternative="two-sided").pvalue)

# ---------- 9. QC ----------
qc = {
 "seed": SEED,
 "cds_genes": {k:int(len(v)) for k,v in CDS.items()},
 "geo_n_HT_files": len(HT),
 "geo_overlap_Pf": {f"{c}_{a}": int(len(set(HTm[(c,a)])&set(CDS['SP001'].index))) for c in ["CM","AA_free"] for a in ["total","polysome"]},
 "mca_overlap": {sp: int(len(set(MCA[sp]['stages'][list(MCA[sp]['stages'])[0]].index.astype(str).str.split('.').str[0])&set(CDS[sp].index))) for sp in MCA},
 "te_rep_corr": {},
 "lcr_flagged_frac": {k: float((v.lcr_frac>0.15).mean()) for k,v in CDS.items()},
 "primary_match": f"{k}/{n}", "match_CI95": list(ci), "binom_p": pv,
 "caveats": ["HS_ tRNA reference is human-annotated; anticodon-aggregated supply proxy only",
             "families are not independent trials; binomial is auxiliary descriptive",
             "n=3 species with expression; no cross-species regression fitted (I3)"]}
for (c,a) in [("CM","total"),("CM","polysome"),("AA_free","total"),("AA_free","polysome")]:
    reps = [pd.Series(HT[(c,a,r)]) for r in ["replicate1","replicate2","replicate3"]]
    idx = sorted(set.intersection(*[set(r.index) for r in reps])&set(CDS["SP001"].index))
    M = np.log2(np.column_stack([r.reindex(idx).values+1 for r in reps]))
    cc = np.corrcoef(M.T)
    qc["te_rep_corr"][f"{c}_{a}"] = float(cc[np.triu_indices(3,1)].mean())
pd.Series(qc).to_json(f"{OUT}/M3-01_qc.json", indent=2)
pd.DataFrame([{"metric":k,"value":str(v)} for k,v in qc.items() if not isinstance(v,dict)]).to_csv(f"{OUT}/M3-01_qc.tsv", sep="\t", index=False)

# key finding printout
sig = te_df[(te_df.q<0.05)]
print("TE model tests:", len(te_df), "sig:", len(sig), flush=True)
print(te_df[te_df.feature.isin(["TTT","TTC","AAA","AAG","ATT","ATC","ATA"]+CAND_AAS)].sort_values("p").head(20).to_string(), flush=True)
print("improvement head:", imp_df.groupby("family").delta.mean().sort_values(ascending=False).head(10).to_string(), flush=True)
print(f"primary match {k}/{n} CI95=({ci[0]:.2f},{ci[1]:.2f}) p={pv:.3f}", flush=True)
print("DONE", flush=True)

# ---------- 10. 04-s10 deliverables ----------
import datetime
inputs = [f"{RAW}/geo/GSE226632/GSE226632_RAW.tar",
 f"{RAW}/mca/pf-ch10x-set1/pf-ch10x-set1.zip", f"{RAW}/mca/pk-ch10x-set1/pk-ch10x-set1.zip",
 f"{RAW}/mca/pb-ss2-set1/pb-ss2-set1.zip",
 f"{ROOT}/data/derived/WP2/M2-01_codon/M2-01_codon_residual_table.tsv",
 f"{ROOT}/data/derived/WP2/M2-02_aa/M2-02_aa_residual_table.tsv",
 f"{ROOT}/data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv",
 f"{ROOT}/data/derived/WP2/M2-01_codon/M2-01_background_explained.tsv"]
for sp, acc in ACC.items():
    inputs.append(f"{RAW}/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna")
with open(f"{OUT}/input_manifest.tsv", "w") as fh:
    fh.write("file\tsha256\tsize\n")
    for p in inputs:
        fh.write(f"{p}\t{sha(p) if os.path.exists(p) else 'MISSING'}\t{os.path.getsize(p) if os.path.exists(p) else 0}\n")
with open(f"{OUT}/params.yaml", "w") as fh:
    fh.write(f"seed: {SEED}\nprimary_contrast: [SP001, SP005]\naux: [SP011]\n"
     f"lcr: {30}aa_window_entropy_lt_1.5_flag_gt_0.15\nnormalization: species_internal_CPM_first\n"
     f"fdr: BH_per_trait_track\ndate: {datetime.date.today().isoformat()}\n")
with open(f"{OUT}/README.md", "w") as fh:
    fh.write("# M3-01 demand\nQ: do 7(+L split-subfamily) codon + 13 AA M2 residuals predict expression-weighted demand beyond genome_gc (H3)?\nPrimary: Pf-vs-Pk direction match. Pb aux descriptive. Transcript!=translation; tRNA layers separate (M3-02).\n")
with open(f"{OUT}/claim_impact.md", "w") as fh:
    fh.write(f"CLM04: primary_match {k}/{n} CI95 {ci} binom_p {pv}. Positive only if direction match beyond GC with stage consistency; else H3 stays HYPOTHESIS.\n")
with open(f"{OUT}/checksums.sha256", "w") as fh:
    import glob as _g
    for p in sorted(_g.glob(f"{OUT}/*")):
        if os.path.isfile(p) and not p.endswith("checksums.sha256"):
            fh.write(f"{sha(p)}  {os.path.basename(p)}\n")
print("DELIVERABLES_WRITTEN", flush=True)
