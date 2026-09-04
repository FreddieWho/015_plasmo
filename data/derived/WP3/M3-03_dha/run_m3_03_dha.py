#!/usr/bin/env python3
"""M3-03 GSE151189 DHA x K13 time-course. Gate C HOLD info-gain task, 保A per D-016.
Q: do candidate families (7+L codons, 13 AAs, LysAAA/AAG focus) respond to DHA
over time, differently by K13 background? Pre-registered primary: early times
(<=8h, minimal stage-delay) DHA response in WT + K13 interaction (DiD).
Writes data/derived/WP3/M3-03_dha/ (+ input_manifest etc per 04-s10).
"""
import os, re, gzip, hashlib, datetime
import numpy as np
import pandas as pd
from scipy import stats as sstats

SEED = 42
ROOT = "/home/huyudi/015_plasmo"
RAW = f"{ROOT}/data/raw"
OUT = f"{ROOT}/data/derived/WP3/M3-03_dha"
os.makedirs(OUT, exist_ok=True)

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
CAND_CODON = ["F","I","P","A","N","D","E","L"]
CAND_AA = ["A","R","N","D","Q","G","I","K","F","P","W","Y","V"]
ACC_PF = "GCF_000002765.6"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s))
    return m.group(0) if m else None

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o]*n/np.arange(1, n+1))[::-1])[::-1]
    return np.minimum(q, 1)

# ---------- 1. Pf CDS codon counts ----------
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
            if "[partial" in hid: return
            if len(s) < 150 or len(s) % 3: return
            if not set(s) <= set("ATGC"): return
            cod = [s[i:i+3] for i in range(0, len(s), 3)]
            if "*" in [STD_CODE.get(c, "*") for c in cod][:-1]: return
            if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
        for line in fh:
            if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
            else: seq.append(line.strip())
        flush()
    rows = {}
    for g, (_, s) in best.items():
        cod = [s[i:i+3] for i in range(0, len(s), 3) if s[i:i+3] in STD_CODE and STD_CODE[s[i:i+3]] != "*"]
        if not cod: continue
        cc = {c: 0 for c in CODONS}
        for c in cod: cc[c] += 1
        rows[g] = {"n_cod": len(cod), **{f"cod_{c}": cc[c] for c in CODONS}}
    df = pd.DataFrame(rows).T
    df.index.name = "gene"
    return df

CDF = parse_cds(ACC_PF)
print("CDS genes:", len(CDF), flush=True)

# ---------- 2. platform probe -> PF3D7 ----------
plat = pd.read_csv("/tmp/m3plat.tsv", sep="\t", usecols=["ID", "ORF"], dtype=str)
plat["gene"] = plat["ORF"].map(gkey)
p2g = plat.dropna(subset=["gene"]).set_index("ID")["gene"].to_dict()
print("probes mapped:", len(p2g), flush=True)

# ---------- 3. series matrix ----------
mx_path = f"{RAW}/geo/GSE151189/GSE151189_series_matrix.txt.gz"
titles, accessions, datastart = [], [], None
with gzip.open(mx_path, "rt") as fh:
    for i, line in enumerate(fh):
        if line.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            accessions = [t.strip().strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            datastart = True
            break
assert len(titles) == len(accessions) and len(titles) > 0
df = pd.read_csv(mx_path, sep="\t", comment="!", header=0, index_col=0,
                 low_memory=False)
df.columns = titles  # align by order
print("matrix:", df.shape, flush=True)

def parse_sample(t):
    dha = "_DHA_" in t
    parts = t.split("_DHA_")[0].split("_") if dha else t.rsplit("_rep", 1)[0].split("_")
    rep = t.rsplit("_rep", 1)[1]
    if dha:
        bg = "_".join(parts[:2]); tm = float(parts[2][:-1])
    else:
        bg = "_".join(parts[:2]); tm = float(parts[2][:-1])
    return bg, tm, dha, rep

meta = pd.DataFrame([parse_sample(t) for t in titles],
                    columns=["bg", "time", "dha", "rep"], index=titles)
print(meta.groupby(["bg", "dha"]).size().to_string(), flush=True)

# probe -> gene, mean multi-probe (log2 space), then weights
g2rows = {}
for pid, g in p2g.items():
    if pid in df.index: g2rows.setdefault(g, []).append(pid)
genes = sorted(set(g2rows) & set(CDF.index))
print("genes with CDS + probe:", len(genes), flush=True)
G = pd.DataFrame({g: df.loc[g2rows[g]].mean(axis=0) for g in genes}).T
G = G.reindex(columns=titles)
W = np.exp2(G.clip(-10, 10))  # relative-abundance proxy vs 3D7 pool

# ---------- 4. per-sample demand shares ----------
FAMS = {a: [c for c in CODONS if STD_CODE[c] == a] for a in AA_FAMS}
dem_rows = []
for s in titles:
    w = W[s].reindex(genes).fillna(0).clip(lower=0)
    tot = (CDF.loc[genes, [f"cod_{c}" for c in CODONS]].multiply(w, axis=0)).sum()
    tot.index = [c.replace("cod_", "") for c in tot.index]
    dc = tot / tot.sum()
    row = {"sample": s}
    for a in AA_FAMS:
        row["AA_" + a] = float(sum(dc[c] for c in FAMS[a]))
    for c in CODONS:
        row["codon_" + c] = float(dc[c])
    row["AAA_AAG_ratio"] = float(dc["AAA"] / max(dc["AAG"], 1e-12))
    dem_rows.append(row)
DEM = pd.DataFrame(dem_rows).set_index("sample")
DEM[["bg", "time", "dha", "rep"]] = meta[["bg", "time", "dha", "rep"]]
DEM.to_csv(f"{OUT}/M3-03_demand_by_sample.tsv", sep="\t")
print("demand written", DEM.shape, flush=True)

# ---------- 5. tests ----------
FEATS = [f"codon_{c}" for a in CAND_CODON for c in FAMS[a]] + ["codon_AAA", "codon_AAG", "AAA_AAG_ratio"] \
    + [f"AA_{a}" for a in CAND_AA]
FEATS = sorted(set(FEATS))
RES_PAIRS = [("Dd2_R539T", "Dd2_WT"), ("Dd2_C580Y", "Dd2_WT"), ("Cam3II_R539T", "Cam3II_revWT")]
rows = []
for (bg, tm), grp in DEM.groupby(["bg", "time"]):
    d = grp[grp.dha]; c = grp[~grp.dha]
    if len(d) == 0 or len(c) == 0: continue
    for f in FEATS:
        t, p = sstats.ttest_ind(d[f].values, c[f].values, equal_var=False)
        rows.append({"bg": bg, "time": tm, "feature": f, "n_dha": len(d), "n_ctrl": len(c),
                     "delta_DHA_CTRL": float(d[f].mean() - c[f].mean()),
                     "p": float(p), "early": bool(tm <= 8)})
eff = pd.DataFrame(rows)
for (bg, tm), idx in eff.groupby(["bg", "time"]).groups.items():
    eff.loc[idx, "q"] = bh(eff.loc[idx, "p"].values)
eff.to_csv(f"{OUT}/M3-03_DHA_effect.tsv", sep="\t", index=False)

# DiD with SE from per-replicate demands (Welch-style z):
did2 = []
for res, wt in RES_PAIRS:
    times = sorted(set(DEM[DEM.bg == res].time) & set(DEM[DEM.bg == wt].time))
    for tm in times:
        for f in FEATS:
            r1 = DEM[(DEM.bg == res) & (DEM.time == tm) & DEM.dha][f].values
            c1 = DEM[(DEM.bg == res) & (DEM.time == tm) & ~DEM.dha][f].values
            r0 = DEM[(DEM.bg == wt) & (DEM.time == tm) & DEM.dha][f].values
            c0 = DEM[(DEM.bg == wt) & (DEM.time == tm) & ~DEM.dha][f].values
            if min(len(r1), len(c1), len(r0), len(c0)) == 0: continue
            d1, d0 = r1.mean() - c1.mean(), r0.mean() - c0.mean()
            se = np.sqrt(r1.var(ddof=1)/len(r1) + c1.var(ddof=1)/len(c1)
                         + r0.var(ddof=1)/len(r0) + c0.var(ddof=1)/len(c0))
            z = (d1 - d0) / max(se, 1e-15)
            p = float(2 * sstats.norm.sf(abs(z)))
            did2.append({"res": res, "wt": wt, "time": tm, "feature": f,
                         "DiD": float(d1 - d0), "se": float(se), "p": p,
                         "early": bool(tm <= 8)})
DID = pd.DataFrame(did2)
for (res, wt, tm), idx in DID.groupby(["res", "wt", "time"]).groups.items():
    DID.loc[idx, "q"] = bh(DID.loc[idx, "p"].values)
DID.to_csv(f"{OUT}/M3-03_K13_DiD.tsv", sep="\t", index=False)
print("DID rows:", len(DID), "early-sig(q<0.05):",
      int(((DID.early) & (DID.q < 0.05)).sum()), flush=True)

# ---------- 6. positive/negative controls ----------
K13 = "PF3D7_1343700"
ctrl_rows = []
if K13 in G.index:
    for (bg, tm), grp in DEM.groupby(["bg", "time"]):
        cols = grp.index
        d = G.loc[K13, cols[grp.dha]].values; c = G.loc[K13, cols[~grp.dha]].values
        if len(d) and len(c):
            t, p = sstats.ttest_ind(d, c, equal_var=False)
            ctrl_rows.append({"gene": "K13", "bg": bg, "time": tm,
                              "delta_log2_DHA_CTRL": float(d.mean() - c.mean()),
                              "p": float(p)})
pd.DataFrame(ctrl_rows).to_csv(f"{OUT}/M3-03_K13_control.tsv", sep="\t", index=False)

# ---------- 7. direction consistency vs M3-01/M3-02 ----------
te = pd.read_csv(f"{ROOT}/data/derived/WP3/M3-01_demand/M3-01_TE_genemodel.tsv", sep="\t")
te_aa = te[(te.trait == "TE_AAfree") & (te.kind == "AA") & (te.q < 0.05)].set_index("feature")["beta"]
cons = []
early_dha = eff[eff.early].groupby("feature").delta_DHA_CTRL.mean()
for a, b in te_aa.items():
    if f"AA_{a}" in early_dha.index:
        cons.append({"AA": a, "TE_beta": float(b),
                     "early_DHA_demand_delta": float(early_dha[f"AA_{a}"]),
                     "sign_match": bool(np.sign(b) == np.sign(early_dha[f"AA_{a}"]))})
pd.DataFrame(cons).to_csv(f"{OUT}/M3-03_consistency_TE.tsv", sep="\t", index=False)

# ---------- 8. QC + deliverables ----------
n_probes_mapped = len(p2g)
qc = {"seed": SEED, "n_samples": len(titles), "n_genes_CDS_probe": len(genes),
      "weights": "2^log2ratio clipped [-10,10], multi-probe mean in log2 space",
      "caveat": "two-color ratio vs 3D7 pool = relative proxy, not absolute TPM; "
                "stage-delay confounding at late times -> early(<=8h) primary",
      "DHA_times_Dd2": sorted(DEM[(DEM.bg.str.startswith('Dd2')) & DEM.dha].time.unique().tolist()),
      "DHA_times_Cam": sorted(DEM[(DEM.bg.str.startswith('Cam')) & DEM.dha].time.unique().tolist())}
pd.Series({k: str(v) for k, v in qc.items()}).to_json(f"{OUT}/M3-03_qc.json", indent=2)

# 保A decision rule (pre-registered D-016)
prim = DID[DID.early & DID.feature.isin(
    ["codon_AAA", "codon_AAG", "AAA_AAG_ratio"] + [f"codon_{c}" for c in FAMS["I"]] + ["AA_K", "AA_I"])]
n_sig = int((prim.q < 0.05).sum())
verdict = ("GO_NARROWED_PF_CANDIDATE" if n_sig >= 1 else "PIVOT_PER_D016")
with open(f"{OUT}/claim_impact.md", "w") as fh:
    fh.write(f"M3-03 early DHAxK13 hits(q<0.05) in LysAAA/AAG+Ile set: {n_sig}. "
             f"Rule D-016: >=1 consistent hit -> {verdict}; else PIVOT. "
             "CLM04 stays HYPOTHESIS until Gate C review.\n")
print(f"primary early hits: {n_sig} -> {verdict}", flush=True)

inputs = [f"{RAW}/geo/GSE151189/GSE151189_series_matrix.txt.gz",
          f"{RAW}/ncbi-datasets/{ACC_PF}/ncbi_dataset/data/{ACC_PF}/cds_from_genomic.fna",
          "/tmp/m3plat.tsv"]
with open(f"{OUT}/input_manifest.tsv", "w") as fh:
    fh.write("file\tsha256\tsize\n")
    for p in inputs:
        fh.write(f"{p}\t{sha(p) if os.path.exists(p) else 'MISSING'}\t{os.path.getsize(p) if os.path.exists(p) else 0}\n")
with open(f"{OUT}/params.yaml", "w") as fh:
    fh.write(f"seed: {SEED}\nprimary: early<=8h DHA response + K13 DiD\nfocus: LysAAA/AAG + Ile + 7+L/13\nweights: exp2_clipped\nfdr: BH_per_bg_time\ndate: {datetime.date.today().isoformat()}\n")
with open(f"{OUT}/README.md", "w") as fh:
    fh.write("# M3-03 DHA x K13\nHOLD info-gain for 保A (D-016): does DHA move candidate-family demand, differently by K13? Early times primary (stage-delay). Two-color log2ratio vs 3D7 pool -> 2^w proxy (approx, documented).\n")
with open(f"{OUT}/checksums.sha256", "w") as fh:
    import glob as _g
    for p in sorted(_g.glob(f"{OUT}/*")):
        if os.path.isfile(p) and not p.endswith("checksums.sha256"):
            fh.write(f"{sha(p)}  {os.path.basename(p)}\n")
print("DONE DELIVERABLES_WRITTEN", flush=True)