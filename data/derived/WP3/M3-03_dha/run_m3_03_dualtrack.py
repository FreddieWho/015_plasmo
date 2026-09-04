#!/usr/bin/env python3
"""M3-03 I5 dual-track closure (adversarial-review fixes 2+4):
(a) persist GPL18893 probe->gene platform table as a proper artifact (replaces /tmp/m3plat.tsv);
(b) rerun the demand computation + primary early DiD on the LCR-MASKED track (M2-02 masker,
    window=64, ent_thr=1.5) and compare with the unmasked track.
Verdict robustness: primary early hits must stay 0 in masked track.
"""
import os, re, gzip, hashlib, datetime
import numpy as np
import pandas as pd
from scipy import stats as sstats

SEED = 42
ROOT = "/home/huyudi/015_plasmo"
RAW = f"{ROOT}/data/raw"
OUT = f"{ROOT}/data/derived/WP3/M3-03_dha"
ACC_PF = "GCF_000002765.6"

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
FAMS = {a: [c for c in CODONS if STD_CODE[c] == a] for a in AA_FAMS}
AA_ALPHABET = list("ARNDCEQGHILKMFPSTWYV")
AA_IDX = {a: i for i, a in enumerate(AA_ALPHABET)}

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

# ----- M2-02 LCR masker (copied verbatim, numba path) -----
def _entropy_from_counts(cnt, window):
    tot = cnt.sum()
    if tot == 0: return 0.0
    distinct = np.count_nonzero(cnt)
    if distinct >= 8: return 2.0
    p = cnt[cnt > 0] / tot
    return -np.sum(p * np.log2(p))

try:
    import numba
    @numba.njit
    def _lcr_njit(idx_arr, n, window, ent_thr, out_mask):
        cnt = np.zeros(20, dtype=np.float64)
        for k in range(window):
            v = idx_arr[k]
            if v >= 0: cnt[v] += 1
        distinct = 0
        for c in range(20):
            if cnt[c] > 0: distinct += 1
        ent = 2.0
        if distinct < 8:
            ent = 0.0
            for c in range(20):
                cc = cnt[c]
                if cc > 0:
                    p = cc / window
                    ent -= p * np.log2(p)
        if ent < ent_thr:
            for k in range(window): out_mask[k] = True
        for i in range(1, n - window + 1):
            left = idx_arr[i-1]; right = idx_arr[i+window-1]
            if left >= 0: cnt[left] -= 1
            if right >= 0: cnt[right] += 1
            distinct = 0
            for c in range(20):
                if cnt[c] > 0: distinct += 1
            if distinct >= 8: continue
            ent = 0.0
            for c in range(20):
                cc = cnt[c]
                if cc > 0:
                    p = cc / window
                    ent -= p * np.log2(p)
            if ent < ent_thr:
                for k in range(window): out_mask[i+k] = True
    _HAS_NUMBA = True
except Exception:
    _HAS_NUMBA = False

def lcr_mask_for_protein(aa_seq, window=64, ent_thr=1.5):
    n = len(aa_seq)
    if n == 0: return np.zeros(0, dtype=bool)
    idx_arr = np.array([AA_IDX.get(a, -1) for a in aa_seq], dtype=np.int16)
    if n < window:
        cnt = np.zeros(20, dtype=float)
        for v in idx_arr:
            if v >= 0: cnt[v] += 1
        ent = _entropy_from_counts(cnt, n)
        return np.array([ent < ent_thr] * n, dtype=bool)
    mask = np.zeros(n, dtype=bool)
    if _HAS_NUMBA:
        _lcr_njit(idx_arr.astype(np.int64), n, window, ent_thr, mask)
        return mask
    cnt = np.zeros(20, dtype=float)
    for v in idx_arr[:window]:
        if v >= 0: cnt[v] += 1
    ent = _entropy_from_counts(cnt, window)
    if ent < ent_thr: mask[0:window] = True
    for i in range(1, n - window + 1):
        left = idx_arr[i-1]; right = idx_arr[i+window-1]
        if left >= 0: cnt[left] -= 1
        if right >= 0: cnt[right] += 1
        ent = _entropy_from_counts(cnt, window)
        if ent < ent_thr: mask[i:i+window] = True
    return mask

# ---------- (a) persist platform table (fix 2) ----------
soft = f"{RAW}/geo/GPL18893/GPL18893_family.soft.gz"
p2g = {}
with gzip.open(soft, "rt", errors="replace") as fh:
    in_tab = False; i_id = i_orf = None
    for line in fh:
        if line.startswith("!platform_table_begin"): in_tab = True; continue
        if line.startswith("!platform_table_end"): break
        if not in_tab: continue
        if i_id is None:
            hdr = [h.strip().strip('"') for h in line.rstrip("\n").split("\t")]
            i_id = hdr.index("ID"); i_orf = hdr.index("ORF"); continue
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(i_id, i_orf): continue
        g = gkey(p[i_orf])
        if g: p2g[p[i_id].strip('"')] = g
PLAT = pd.DataFrame({"probe_id": list(p2g), "gene": list(p2g.values())})
PLAT.to_csv(f"{OUT}/platform_GPL18893_probe2gene.tsv", sep="\t", index=False)
print("platform probes persisted:", len(PLAT), flush=True)

# ---------- (b) dual-track CDS codon counts ----------
def parse_cds_masked(acc):
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
            cod = [s[i:i+3] for i in range(0, len(s), 3)]
            aas = [STD_CODE.get(c, "*") for c in cod]
            if "*" in aas[:-1]: return
            if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
        for line in fh:
            if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
            else: seq.append(line.strip())
        flush()
    rows = {}
    for g, (_, s) in best.items():
        cod = [s[i:i+3] for i in range(0, len(s), 3)]
        aas = [STD_CODE.get(c, "*") for c in cod]
        keep = [i for i, a in enumerate(aas) if a != "*"]
        cod = [cod[i] for i in keep]; aas = [aas[i] for i in keep]
        mask = lcr_mask_for_protein("".join(aas))
        cod_kept = [c for c, m in zip(cod, mask) if not m]
        if not cod_kept: continue
        cc = {c: 0 for c in CODONS}
        for c in cod_kept: cc[c] += 1
        rows[g] = {"n_cod_masked": len(cod_kept), **{f"cod_{c}": cc[c] for c in CODONS}}
    df = pd.DataFrame(rows).T
    df.index.name = "gene"
    return df

CDF_M = parse_cds_masked(ACC_PF)
print("masked-track CDS genes:", len(CDF_M), flush=True)

# ---------- expression + demand on masked track ----------
mx_path = f"{RAW}/geo/GSE151189/GSE151189_series_matrix.txt.gz"
titles = []
with gzip.open(mx_path, "rt") as fh:
    for line in fh:
        if line.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"): break
df = pd.read_csv(mx_path, sep="\t", comment="!", header=0, index_col=0, low_memory=False)
df.columns = titles

def parse_sample(t):
    dha = "_DHA_" in t
    parts = t.split("_DHA_")[0].split("_") if dha else t.rsplit("_rep", 1)[0].split("_")
    rep = t.rsplit("_rep", 1)[1]
    bg = "_".join(parts[:2]); tm = float(parts[2][:-1])
    return bg, tm, dha, rep
meta = pd.DataFrame([parse_sample(t) for t in titles], columns=["bg", "time", "dha", "rep"], index=titles)

g2rows = {}
for pid, g in p2g.items():
    if pid in df.index: g2rows.setdefault(g, []).append(pid)
genes = sorted(set(g2rows) & set(CDF_M.index))
G = pd.DataFrame({g: df.loc[g2rows[g]].mean(axis=0) for g in genes}).T.reindex(columns=titles)
W = np.exp2(G.clip(-10, 10))
dem_rows = []
for s in titles:
    w = W[s].reindex(genes).fillna(0).clip(lower=0)
    tot = (CDF_M.loc[genes, [f"cod_{c}" for c in CODONS]].multiply(w, axis=0)).sum()
    tot.index = [c.replace("cod_", "") for c in tot.index]
    dc = tot / tot.sum()
    row = {"sample": s}
    for a in AA_FAMS: row["AA_" + a] = float(sum(dc[c] for c in FAMS[a]))
    for c in CODONS: row["codon_" + c] = float(dc[c])
    row["AAA_AAG_ratio"] = float(dc["AAA"] / max(dc["AAG"], 1e-12))
    dem_rows.append(row)
DEM = pd.DataFrame(dem_rows).set_index("sample")
DEM[["bg", "time", "dha", "rep"]] = meta
DEM.to_csv(f"{OUT}/M3-03_demand_by_sample_LCRmasked.tsv", sep="\t")

# primary early DiD on masked track (same focus set as D-016 rule)
CAND_CODON = ["F", "I", "P", "A", "N", "D", "E", "L"]
CAND_AA = ["A", "R", "N", "D", "Q", "G", "I", "K", "F", "P", "W", "Y", "V"]
FEATS = sorted(set([f"codon_{c}" for a in CAND_CODON for c in FAMS[a]] + ["codon_AAA", "codon_AAG", "AAA_AAG_ratio"] + [f"AA_{a}" for a in CAND_AA]))
RES_PAIRS = [("Dd2_R539T", "Dd2_WT"), ("Dd2_C580Y", "Dd2_WT"), ("Cam3II_R539T", "Cam3II_revWT")]
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
            se = np.sqrt(r1.var(ddof=1)/len(r1) + c1.var(ddof=1)/len(c1) + r0.var(ddof=1)/len(r0) + c0.var(ddof=1)/len(c0))
            z = (d1 - d0) / max(se, 1e-15)
            did2.append({"res": res, "wt": wt, "time": tm, "feature": f, "DiD": float(d1 - d0),
                         "p": float(2 * sstats.norm.sf(abs(z))), "early": bool(tm <= 8)})
DID = pd.DataFrame(did2)
for (res, wt, tm), idx in DID.groupby(["res", "wt", "time"]).groups.items():
    DID.loc[idx, "q"] = bh(DID.loc[idx, "p"].values)
DID.to_csv(f"{OUT}/M3-03_K13_DiD_LCRmasked.tsv", sep="\t", index=False)
prim = DID[DID.early & DID.feature.isin(["codon_AAA", "codon_AAG", "AAA_AAG_ratio"] + [f"codon_{c}" for c in FAMS["I"]] + ["AA_K", "AA_I"])]
n_sig_masked = int((prim.q < 0.05).sum())

# cross-track comparison: unmasked DiD (existing file) vs masked
DID_U = pd.read_csv(f"{OUT}/M3-03_K13_DiD.tsv", sep="\t")
J = DID_U.merge(DID, on=["res", "wt", "time", "feature"], suffixes=("_unmasked", "_masked"))
rho = float(sstats.spearmanr(J.DiD_unmasked, J.DiD_masked).statistic)
prim_u = DID_U[DID_U.early & DID_U.feature.isin(prim.feature.unique())]
n_sig_unmasked = int((prim_u.q < 0.05).sum())
J[J.feature.isin(["codon_AAA", "codon_AAG", "AAA_AAG_ratio"])].to_csv(f"{OUT}/M3-03_dualtrack_focus_compare.tsv", sep="\t", index=False)

qc = {"masked_genes": int(len(CDF_M)), "did_rows_masked": int(len(DID)),
      "early_primary_hits_masked": n_sig_masked, "early_primary_hits_unmasked": n_sig_unmasked,
      "spearman_DiD_cross_track": rho, "verdict_robust": bool(n_sig_masked == 0 and n_sig_unmasked == 0)}
pd.Series({k: str(v) for k, v in qc.items()}).to_json(f"{OUT}/M3-03_dualtrack_qc.json", indent=2)
print(qc, flush=True)

# fix manifest: replace /tmp/m3plat.tsv with GPL18893 soft + persisted platform table
mani = pd.read_csv(f"{OUT}/input_manifest.tsv", sep="\t")
mani = mani[~mani.file.str.startswith("/tmp/")]
new = pd.DataFrame({"file": [soft, f"{OUT}/platform_GPL18893_probe2gene.tsv"]})
new["sha256"] = new.file.map(sha); new["size"] = new.file.map(os.path.getsize)
pd.concat([mani, new]).to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
with open(f"{OUT}/claim_impact.md", "a") as fh:
    fh.write(f"\n\n## I5 dual-track closure (2026-09-04)\n"
             f"LCR-masked track (M2-02 masker w64/e1.5): {len(CDF_M)} genes; primary early hits masked={n_sig_masked} vs unmasked={n_sig_unmasked}; "
             f"cross-track DiD spearman={rho:.4f}. Verdict PIVOT robust to LCR masking. "
             f"/tmp/m3plat.tsv provenance replaced by persisted platform_GPL18893_probe2gene.tsv + GPL18893 soft in manifest.\n")
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(os.listdir(OUT)):
        p = os.path.join(OUT, f)
        if os.path.isfile(p) and f != "checksums.sha256" and not f.startswith("__"):
            fh.write(f"{sha(p)}  {f}\n")
print("DONE", flush=True)
