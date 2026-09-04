#!/usr/bin/env python
# L1c (ring 1): tRNA supply side from GSE226632
# (a) tRNA abundance: AA-depletion vs control vs +halofuginone, per AA-anticodon family
# (b) charging proxy per stage: family share in periodate-treated / untreated (R/T/S x 3rep)
# caveat: HS_ namespace (adversarial review B3) -> results labeled QUARANTINED-REF, family-level only
import gzip, hashlib, json, os, re, tarfile
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/huyudi/015_plasmo"
OUT = f"{ROOT}/data/derived/WP4/L1c_trna_charging"
os.makedirs(OUT, exist_ok=True)

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()

def fam(name):
    m = re.match(r"HS_([A-Za-z]{3})-([ACGT]{3})", name)
    return f"{m.group(1)}-{m.group(2)}" if m else None

tar_path = f"{ROOT}/data/raw/geo/GSE226632/GSE226632_RAW.tar"
tar = tarfile.open(tar_path)
COND, STAGE, PERI = {}, {}, {}
for m in tar.getnames():
    b = m.split("/")[-1]
    if "tRNA" not in b: continue
    fh = tar.extractfile(m); d = {}
    with gzip.open(fh, "rt") as z:
        for line in z:
            p = line.split()
            if len(p) >= 2 and p[0].startswith("HS_"): d[p[0]] = float(p[1])
    m0 = re.search(r"tRNA-(control|AA|HF)_(\d)_-count", b)
    m1 = re.search(r"tRNA_([RTS])(\d)_count", b)
    m2 = re.search(r"tRNA-([RTS])(\d)-periodate-count", b)
    if m0: COND[(m0.group(1), m0.group(2))] = d
    elif m2: PERI[(m2.group(1), m2.group(2))] = d
    elif m1: STAGE[(m1.group(1), m1.group(2))] = d
tar.close()
print("cond:", sorted(COND), "\nstage:", sorted(STAGE), "\nperi:", sorted(PERI), flush=True)

def share(d):
    fams = {}
    for k, v in d.items():
        f = fam(k)
        if f: fams[f] = fams.get(f, 0) + v
    t = sum(fams.values())
    return {f: v / t for f, v in fams.items()}

# (a) abundance under stress
fams = sorted({fam(k) for d in COND.values() for k in d} - {None})
rows = []
for f in fams:
    c = np.mean([share(COND[("control", r)]).get(f, 0) for r in "123"])
    a = np.mean([share(COND[("AA", r)]).get(f, 0) for r in "123"])
    h = np.mean([share(COND[("HF", r)]).get(f, 0) for r in "123"])
    if c == 0 and a == 0: continue
    rows.append({"family": f, "share_control": c, "share_AAdepl": a, "share_HF": h,
                 "log2_AA_vs_ctrl": np.log2((a + 1e-6) / (c + 1e-6)),
                 "log2_HF_vs_ctrl": np.log2((h + 1e-6) / (c + 1e-6))})
AB = pd.DataFrame(rows).sort_values("log2_AA_vs_ctrl")
AB["rank_AAdepl"] = AB.log2_AA_vs_ctrl.rank()
AB.to_csv(f"{OUT}/L1c_abundance_by_family.tsv", sep="\t", index=False)
asn = AB[AB.family.str.startswith("Asn")]
print("Asn families:\n", asn.to_string(), flush=True)

# (b) charging proxy per stage: share(periodate)/share(untreated)
rows = []
for f in fams:
    for st in "RTS":
        rats = []
        for r in "123":
            u = share(STAGE[(st, r)]).get(f, 0); p = share(PERI[(st, r)]).get(f, 0)
            if u > 0: rats.append(p / u)
        if rats:
            rows.append({"family": f, "stage": st, "charging_proxy": float(np.mean(rats)),
                         "cv_reps": float(np.std(rats) / np.mean(rats)) if len(rats) > 1 else np.nan})
CH = pd.DataFrame(rows)
CH["rank_within_stage"] = CH.groupby("stage").charging_proxy.rank(ascending=False)
CH.to_csv(f"{OUT}/L1c_charging_by_stage.tsv", sep="\t", index=False)
asn_ch = CH[CH.family.str.startswith("Asn")]
print("Asn charging:\n", asn_ch.to_string(), flush=True)

# verdict: Asn-GTT abundance change rank + charging landscape
ab_rank = float(asn[asn.family == "Asn-GTT"].rank_AAdepl.iloc[0]) if len(asn[asn.family == "Asn-GTT"]) else np.nan
n_fam = len(AB)
verdict = {"Asn_GTT_abundance_rank_low_is_depleted": f"{ab_rank}/{n_fam}",
           "Asn_charging_by_stage": asn_ch[["stage", "charging_proxy", "rank_within_stage"]].to_dict("records")}
json.dump(verdict, open(f"{OUT}/L1c_qc.json", "w"), indent=2)
open(f"{OUT}/claim_impact.md", "w").write(f"""# L1c claim impact (ring 1 of AAT mechanism chain; exploratory)

- date: 2026-09-04
- Asn-GTT (reads AAT) abundance rank under AA depletion: {verdict['Asn_GTT_abundance_rank_low_is_depleted']} (1 = most depleted)
- Asn charging proxy by stage: {json.dumps(verdict['Asn_charging_by_stage'])}
- HF arm = halofuginone (ProRS inhibitor) - Pro-family abundance shift is built-in positive control.
- QUARANTINE: HS_ reference namespace (adversarial review B3) - family-level descriptive only, NOT for formal claims.
- interpretation for AAT chain: see TODO/STATUS after user review.
""")
mani = pd.DataFrame({"file": [tar_path]})
mani["sha256"] = mani.file.map(sha); mani["size"] = mani.file.map(os.path.getsize)
mani.to_csv(f"{OUT}/input_manifest.tsv", sep="\t", index=False)
outs = [f for f in os.listdir(OUT) if f != "checksums.sha256" and os.path.isfile(os.path.join(OUT, f))]
with open(f"{OUT}/checksums.sha256", "w") as fh:
    for f in sorted(outs):
        fh.write(f"{sha(os.path.join(OUT, f))}  {f}\n")
print("DONE", flush=True)
