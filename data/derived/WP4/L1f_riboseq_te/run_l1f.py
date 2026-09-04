#!/usr/bin/env python
# L1f: third-dataset translation-layer validation (GSE58402, Caro 2014 Ribo-seq, IDC 5 stages)
# chain prediction: Asn-rich genes' TE RISES across IDC (ring -> schizont), tracking Asn-tRNA charging (L1c: R 1.7 -> T 4.8 -> S 3.9)
import gzip, hashlib, json, os, re, tarfile
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L1f_riboseq_te"
os.makedirs(OUT, exist_ok=True)
TAR = f"{ROOT}/data/raw/geo/GSE58402/GSE58402_RAW.tar"
ALIAS = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_GeneAliases.txt"
STAGES = ["2h_ring", "10h_latering", "19h_troph", "31h_latetroph", "46h_schizont"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

# ---------- alias map: old v3 ids -> PF3D7 ----------
a2g = {}
for line in open(ALIAS):
    p = line.rstrip("\n").split("\t")
    if len(p) < 2 or not p[0].startswith("PF3D7"): continue
    for a in p[1:]:
        if a: a2g[a.lower()] = p[0]

# ---------- parse rpkm from tar ----------
def read_rpkm(name):
    d = {}
    with gzip.open(tar.extractfile(name), "rt") as z:
        for line in z:
            p = line.split()
            if len(p) >= 2:
                g = a2g.get(p[0].lower())
                if g:
                    try: d[g] = max(d.get(g, 0), float(p[1]))
                    except ValueError: pass
    return d

tar = tarfile.open(TAR)
mRNA, RPF = {}, {}
for i in range(1, 6):
    mRNA[STAGES[i-1]] = read_rpkm(f"GSM141029{i}_mRNA_{i}_rpkm.txt.gz")
    gsm = 295 + i  # 296..300
    RPF[STAGES[i-1]] = read_rpkm(f"GSM1410{gsm}_ribosome_footprints_{i}_rpkm.txt.gz")
tar.close()
print("mapped mRNA genes per stage:", {k: len(v) for k, v in mRNA.items()}, flush=True)

# ---------- TE per stage ----------
genes = set.intersection(*[set(mRNA[s]) for s in STAGES], *[set(RPF[s]) for s in STAGES])
rows = []
for g in genes:
    rec = {"gene": g}
    ok = True
    for s in STAGES:
        m, r = mRNA[s].get(g, 0), RPF[s].get(g, 0)
        if m < 1 or r < 0.5: ok = False; break
        rec[f"TE_{s}"] = r / m
    if ok: rows.append(rec)
TE = pd.DataFrame(rows).set_index("gene")
print("genes with TE at all 5 stages:", len(TE), flush=True)

# ---------- features (same as L1e) ----------
STD = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
               "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
acc = "GCF_000002765.6"
CDS = f"{ROOT}/data/raw/ncbi-datasets/{acc}/ncbi_dataset/data/{acc}/cds_from_genomic.fna"
best = {}
with open(CDS) as fh:
    hid, seq = None, []
    def flush():
        if hid is None: return
        m = re.search(r"\[locus_tag=([^\]]+)\]", hid)
        if not m: return
        g = re.search(r"PF3D7_\d+", m.group(1))
        if not g: return
        g = g.group(0)
        s = "".join(seq).upper().replace("U", "T")
        if "[partial" in hid or len(s) < 150 or len(s) % 3 or not set(s) <= set("ATGC"): return
        cod = [s[i:i+3] for i in range(0, len(s), 3)]
        if "*" in [STD.get(c, "*") for c in cod][:-1]: return
        if len(s) > len(best.get(g, ("", ""))[1]): best[g] = (hid, s)
    for line in fh:
        if line.startswith(">"): flush(); hid, seq = line.rstrip(), []
        else: seq.append(line.strip())
    flush()
feat = {}
for g, (_, s) in best.items():
    cod = [s[i:i+3] for i in range(0, len(s), 3) if STD.get(s[i:i+3], "*") != "*"]
    if cod:
        cc = Counter(cod); n = len(cod)
        feat[g] = {"AAT_freq": cc["AAT"] / n, "Asn_content": (cc["AAT"] + cc["AAC"]) / n}
FE = pd.DataFrame(feat).T

M = TE.join(FE, how="inner")
M["TE_late_over_ring"] = np.log2(M["TE_46h_schizont"] / M["TE_2h_ring"])
M["TE_troph_over_ring"] = np.log2(M["TE_19h_troph"] / M["TE_2h_ring"])
q90 = M.Asn_content.quantile(0.9)
hi = M.Asn_content >= q90
res = {
    "n_genes": len(M),
    "spearman_Asn_vs_TE46over2": float(stats.spearmanr(M.Asn_content, M.TE_late_over_ring).statistic),
    "spearman_p": float(stats.spearmanr(M.Asn_content, M.TE_late_over_ring).pvalue),
    "spearman_Asn_vs_TE19over2": float(stats.spearmanr(M.Asn_content, M.TE_troph_over_ring).statistic),
    "spearman19_p": float(stats.spearmanr(M.Asn_content, M.TE_troph_over_ring).pvalue),
    "median_late_over_ring_AsnTop10": float(M.loc[hi, "TE_late_over_ring"].median()),
    "median_late_over_ring_rest": float(M.loc[~hi, "TE_late_over_ring"].median()),
    "MW_p_AsnTop10_rises": float(stats.mannwhitneyu(M.loc[hi, "TE_late_over_ring"], M.loc[~hi, "TE_late_over_ring"], alternative="greater").pvalue),
    "spearman_AAT_vs_TE46over2": float(stats.spearmanr(M.AAT_freq, M.TE_late_over_ring).statistic),
}
prof = {s: {"AsnTop10": float(np.log2(M.loc[hi, f"TE_{s}"].median())), "rest": float(np.log2(M.loc[~hi, f"TE_{s}"].median()))} for s in STAGES}
res["median_log2TE_profile"] = prof
json.dump(res, open(f"{OUT}/L1f_qc.json", "w"), indent=2)
M.to_csv(f"{OUT}/L1f_te_by_stage.tsv", sep="\t")
print(json.dumps(res, indent=1)[:1200], flush=True)

open(f"{OUT}/claim_impact.md", "w").write(f"""# L1f claim impact (third-dataset translation-layer validation; exploratory)

- date: 2026-09-04
- data: GSE58402 Caro 2014 Ribo-seq+mRNA, W2 strain, 5 IDC stages (2/10/19/31/46 hpi), processed RPKM.
- prediction from chain: Asn-rich genes' TE rises ring->schizont, tracking Asn-tRNA charging (L1c: R 1.7 -> T 4.8 -> S 3.9).
- spearman(Asn content, TE 46h/2h) = {res['spearman_Asn_vs_TE46over2']:.3f} (p={res['spearman_p']:.2g});
  Asn-top10% median late/ring TE shift = {res['median_late_over_ring_AsnTop10']:.3f} vs rest {res['median_late_over_ring_rest']:.3f}, MW p={res['MW_p_AsnTop10_rises']:.2g}.
- per-stage profile in L1f_qc.json median_log2TE_profile.
- caveat: W2 strain, 2014-era RPKM quant, no replicates per stage (single timecourse) -> descriptive support level.
- claim impact: independent translation-layer evidence for/against the Asn axis; feeds condition-1 closure.
""")
mani = pd.DataFrame({"file": [TAR, ALIAS, CDS]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(f"{OUT}/{f}")]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
