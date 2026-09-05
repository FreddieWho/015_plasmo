#!/usr/bin/env python
"""M4R-A ROUND2 step5: PbApiAP2 family phenotyping map + Pf-ortholog Asn crosswalk (map-only)."""
import re
import numpy as np
import pandas as pd
import openpyxl
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4R/M4R-A"

def gkey(s):
    m = re.search(r"PF3D7_\d+", str(s)); return m.group(0) if m else None

wb = openpyxl.load_workbook("/tmp/pbinv_Data_1.xlsx", read_only=True, data_only=True)
ws = wb["AP2KOscreen"]
rows = list(ws.values)
hdr = [str(v) if v else "" for v in rows[2]]
gi = hdr.index("Gene_ID"); pf = hdr.index("P. falciparum orthologue")
nm = hdr.index("Name"); dv = hdr.index("Gene disruption verified")
name2id, name2pf, name2ko = {}, {}, {}
for r in rows[3:]:
    if not r[gi]: continue
    name = str(r[nm]).strip() if r[nm] and str(r[nm]).strip() != "-" else None
    key = name or str(r[gi])
    name2id[key] = str(r[gi])
    name2pf[key] = gkey(r[pf])
    name2ko[key] = (str(r[dv]).strip().upper() == "YES")
print("ko-screen genes:", len(name2id), flush=True)

wp = wb["Phenotyping"]
pr = list(wp.values)
stages = [str(v) for v in pr[2][2::2]]
print("stages:", stages, flush=True)
recs = []
for r in pr[4:]:
    if not r[1]: continue
    key = str(r[1]).strip()
    vals = {}
    for j, st in enumerate(stages):
        v = r[2 + 2 * j]
        vals[st] = float(v) if isinstance(v, (int, float)) else float("nan")
    recs.append({"mutant": key, "pb_gene": name2id.get(key), "pf_ortholog": name2pf.get(key),
                 "ko_verified": name2ko.get(key), **vals})
PHE = pd.DataFrame(recs)
# transmission-blocked: oocyst or sporozoite ~0 while asexual viable
asex = [c for c in PHE.columns if "Asexual" in c][0]
oo = [c for c in PHE.columns if "Oocyst" in c][0]
sp = [c for c in PHE.columns if "Sporozoite" in c][0]
PHE["transmission_blocked"] = (PHE[[oo, sp]].min(axis=1, skipna=True) < 0.05) & (PHE[asex] > 0.2)
PHE.to_csv(f"{OUT}/M4RA_pbapiap2_family_map.tsv", sep="\t", index=False)
print(PHE[["mutant", "pf_ortholog", "ko_verified", "transmission_blocked"]].to_string(), flush=True)

# crosswalk: Pf ortholog Asn of transmission-blocking vs other ApiAP2 KOs
FAA = f"{ROOT}/data/raw/veupathdb/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7_AnnotatedProteins.fasta"
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
asn = {g: s.count("N") / len(s) for g, s in prots.items()}
sub = PHE[PHE.pf_ortholog.notna() & PHE.ko_verified].copy()
sub["pf_asn"] = sub.pf_ortholog.map(asn)
sub = sub[sub.pf_asn.notna()]
print("mapped KO-verified mutants:", len(sub), flush=True)
a = sub[sub.transmission_blocked].pf_asn; b = sub[~sub.transmission_blocked].pf_asn
if len(a) >= 3 and len(b) >= 3:
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
else:
    u, p = float("nan"), float("nan")
print(f"transmission-blocked n={len(a)} med={a.median() if len(a) else float('nan'):.4f} "
      f"vs other n={len(b)} med={b.median() if len(b) else float('nan'):.4f} MW p={p}", flush=True)
pd.DataFrame([{"test": "MW Pf-ortholog Asn (transmission-blocked vs other ApiAP2 KOs)",
               "n_blocked": len(a), "n_other": len(b),
               "med_blocked": float(a.median()) if len(a) else float("nan"),
               "med_other": float(b.median()) if len(b) else float("nan"),
               "MW_p": float(p),
               "note": "Pb phenotype x Pf composition; map-only, cross-species mapping disclosed"}]
             ).to_csv(f"{OUT}/M4RA_pbapiap2_xwalk.tsv", sep="\t", index=False)
print("STEP5 DONE", flush=True)
