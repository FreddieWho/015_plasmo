#!/usr/bin/env python
# L1h (ring 3c): Mok 2021 genotype-context protein test - do Asn-rich proteins tilt in ART-resistant K13 mutants?
# chain prediction: ART-resistant state is a chronic stress-adapted state -> same Asn tilt as starvation (L1g +0.167)
import hashlib, json, os, re
from collections import Counter
import numpy as np
import pandas as pd
import openpyxl
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L1h_mok_genotype_context"
os.makedirs(OUT, exist_ok=True)
MOK = f"{ROOT}/data/manual_inbox/M04_NATCOMM2021_MOK/41467_2020_20805_MOESM6_ESM.xlsx"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

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
q90 = FE.Asn_content.quantile(0.9)

# block map: (name, detected_col, fc_col, p_col) ; fc/p positions within each mutant sub-block
BLOCKS = [("Rings_E1_R539T", 11, 14, 16), ("Rings_E1_C580Y", 11, 18, 20),
          ("Rings_E2_R539T", 23, 26, 28), ("Rings_E2_C580Y", 23, 30, 32),
          ("Troph_E1_R539T", 35, 38, 40), ("Troph_E1_C580Y", 35, 42, 44),
          ("Troph_E2_R539T", 47, 50, 52)]
wb = openpyxl.load_workbook(MOK, read_only=True)
ws = wb["Supplmentary Data 3"]
data = list(ws.iter_rows(min_row=6, values_only=True))
wb.close()

def num(v):
    try: return float(v)
    except (TypeError, ValueError): return np.nan

res = []
long_recs = []
for name, dcol, fcol, pcol in BLOCKS:
    recs = []
    for r in data:
        g = r[0]
        if not g or not str(g).startswith("PF3D7"): continue
        if str(r[dcol]) != "Y": continue
        fc = num(r[fcol]); p = num(r[pcol])
        if not np.isfinite(fc): continue
        recs.append({"gene": str(g), "log2FC": fc, "p": p})
    D = pd.DataFrame(recs).merge(FE, left_on="gene", right_index=True, how="inner")
    if len(D) < 50:
        print(name, "too few:", len(D)); continue
    D["AsnTop10"] = D.Asn_content >= q90
    hi, lo = D[D.AsnTop10], D[~D.AsnTop10]
    sp = stats.spearmanr(D.Asn_content, D.log2FC)
    sig = D[(D.p < 0.05)]
    sig_hi = sig[sig.AsnTop10]
    res.append({"contrast": name, "n": len(D), "n_AsnTop10": int(D.AsnTop10.sum()),
                "median_FC_AsnTop10": float(hi.log2FC.median()), "median_FC_rest": float(lo.log2FC.median()),
                "MW_p": float(stats.mannwhitneyu(hi.log2FC, lo.log2FC).pvalue),
                "spearman_Asn_FC": float(sp.statistic), "spearman_p": float(sp.pvalue),
                "n_sig": len(sig), "frac_AsnTop10_in_sig": float(len(sig_hi) / max(len(sig), 1))})
    for _, row in D.iterrows():
        long_recs.append({"contrast": name, "gene": row.gene, "log2FC": row.log2FC, "p": row.p,
                          "Asn_content": row.Asn_content})
R = pd.DataFrame(res)
R.to_csv(f"{OUT}/L1h_tests.tsv", sep="\t", index=False)
pd.DataFrame(long_recs).to_csv(f"{OUT}/L1h_long.tsv", sep="\t", index=False)
print(R.to_string(), flush=True)

# meta: does the sign agree across blocks/mutants?
signs = [r["spearman_Asn_FC"] for r in res]
meta = {"n_contrasts": len(res), "n_positive": int(sum(1 for x in signs if x > 0)),
        "median_spearman": float(np.median(signs))}
json.dump({"tests": res, "meta": meta}, open(f"{OUT}/L1h_qc.json", "w"), indent=2)
open(f"{OUT}/claim_impact.md", "w").write(f"""# L1h claim impact (ring 3c: Mok 2021 genotype context; exploratory)

- date: 2026-09-04
- data: Mok 2021 Nat Commun Suppl Data 3 (PXD019612 processed DE): Cam3.II R539T/C580Y vs WT,
  rings x2 expts + trophs x2 expts, TMT.
- question: do Asn-rich proteins tilt in ART-resistant K13 mutants the way they do under starvation (L1g +0.167)?
- result: {meta['n_positive']}/{meta['n_contrasts']} contrasts positive spearman; median rho={meta['median_spearman']:.3f}.
  Per-contrast stats in L1h_tests.tsv.
- interpretation: genotype-context (chronic adapted state), NOT an acute stress response; sign-consistent
  blocks support the axis being engaged in the resistant phenotype; sign-flipped blocks weaken it.
- claim impact: context evidence for/against Asn axis in ART resistance; no claim change.
""")
mani = pd.DataFrame({"file": [MOK, CDS]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(f"{OUT}/{f}")]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
