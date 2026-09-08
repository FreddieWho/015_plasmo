#!/usr/bin/env python
"""MS Fig4 (Route A, bounded): (a) formal-vs-regex enrichment forest;
(b) counterfactual strata; (c) GCN5 + Pb map + tripole honest box.
Plot-only from in-hand TSVs. English labels, colorblind-safe + shape double-encoding.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/huyudi/015_plasmo"
A = f"{ROOT}/data/derived/WP4R/M4R-A"
MS = f"{ROOT}/manuscript/figures"
os.makedirs(f"{MS}/scripts", exist_ok=True)
os.makedirs(f"{MS}/panels", exist_ok=True)
os.makedirs(f"{MS}/legends", exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
C_FORMAL, C_REGEX = "#0072B2", "#E69F00"  # blue vs orange + circle vs square

# ---------------- Fig4a forest ----------------
en = pd.read_csv(f"{A}/M4RA_interpro_enrichment.tsv", sep="\t")
want = ["ApiAP2 [InterPro-formal]", "ApiAP2 [regex-baseline]",
        "PUF_strict [InterPro-formal]", "PUF [regex-baseline]",
        "RNA_binding_broad [InterPro-formal]", "RNA_binding_broad [regex-baseline]",
        "chromatin_reg [InterPro-formal]", "chromatin_reg [regex-baseline]",
        "CCR4_NOT [InterPro-formal]", "CCR4_NOT [regex-baseline]",
        "transcription_reg [InterPro-formal]", "transcription_reg [regex-baseline]",
        "sexual_gametocyte [InterPro-formal]", "sexual_gametocyte [regex-baseline]",
        "proteostasis [InterPro-formal]", "proteostasis [regex-baseline]",
        "PUF_repeatbroad [InterPro-formal]"]
sub = en[en["set"].isin(want)].copy()
labels = {"ApiAP2 [InterPro-formal]": "ApiAP2 (formal, n=24)",
          "ApiAP2 [regex-baseline]": "ApiAP2 (regex, n=28)",
          "PUF_strict [InterPro-formal]": "PUF strict (IPR, n=2)*",
          "PUF [regex-baseline]": "PUF strict (regex, n=2)*",
          "RNA_binding_broad [InterPro-formal]": "RNA-binding broad (formal, n=467)",
          "RNA_binding_broad [regex-baseline]": "RNA-binding broad (regex, n=201)",
          "chromatin_reg [InterPro-formal]": "chromatin (formal, n=176)",
          "chromatin_reg [regex-baseline]": "chromatin (regex, n=131)",
          "CCR4_NOT [InterPro-formal]": "CCR4-NOT (formal, n=8)",
          "CCR4_NOT [regex-baseline]": "CCR4-NOT (regex, n=12)",
          "transcription_reg [InterPro-formal]": "transcription (formal, n=4, 0 hits)",
          "transcription_reg [regex-baseline]": "transcription (regex, n=89)",
          "sexual_gametocyte [InterPro-formal]": "sexual/gametocyte (formal == AP2 set)",
          "sexual_gametocyte [regex-baseline]": "sexual/gametocyte (regex, n=45)",
          "proteostasis [InterPro-formal]": "proteostasis (formal, n=318)",
          "proteostasis [regex-baseline]": "proteostasis (regex, n=166)",
          "PUF_repeatbroad [InterPro-formal]": "PUF repeat-broad (n=79)"}
go_path = f"{A}/M4RA_l010_enrichment.tsv"
if os.path.exists(go_path):
    gl = pd.read_csv(go_path, sep="\t")
    grow = gl[(gl["set"] == "TF_DNAbinding") & (gl["stratum"] == "nonIEA")].iloc[0]
    go_row = pd.DataFrame([{"set": "TF [GO-curated nonIEA]", "n_set": int(grow["n"]),
                             "n_in_top10": int(grow["in_top"]), "OR": float(grow["OR"]),
                             "lo": float(grow["lo"]), "hi": float(grow["hi"]),
                             "p": float(grow["p"]), "q_bh": float(grow["q_bh"])}])
    sub = pd.concat([sub, go_row], ignore_index=True)
    labels["TF [GO-curated nonIEA]"] = f"TF (GO curated non-IEA, n={int(grow['n'])})*"
    order = [s for s in want if s in set(sub["set"])] + ["TF [GO-curated nonIEA]"]
else:
    order = [s for s in want if s in set(sub["set"])]
sub["ylab"] = sub["set"].map(labels)
sub.to_csv(f"{MS}/panels/fig4a_forest.tsv", sep="\t", index=False)

fig, ax = plt.subplots(figsize=(8.2, 6.4))
ys = np.arange(len(order))
for i, s in enumerate(order):
    r = sub[sub["set"] == s].iloc[0]
    formal = "regex-baseline" not in s and "repeatbroad" not in s
    col = "#009E73" if "GO-curated" in s else (C_FORMAL if ("[InterPro-formal]" in s or "repeatbroad" in s) else C_REGEX)
    mk = "^" if "GO-curated" in s else ("o" if "[InterPro-formal]" in s or "repeatbroad" in s else "s")
    if r["OR"] > 0 and np.isfinite(r["lo"]) and np.isfinite(r["hi"]) and r["lo"] > 0:
        lo = max(r["lo"], 0.12)
        ax.errorbar(r["OR"], ys[i], xerr=[[r["OR"] - lo], [min(r["hi"], 160) - r["OR"]]],
                    fmt=mk, color=col, ecolor=col, capsize=3, ms=5)
    else:
        ax.plot(0.14, ys[i], marker="x", color="grey", ms=6)
    ax.text(165, ys[i], f"q={r['q_bh']:.1e}" if r["q_bh"] < 0.01 else f"q={r['q_bh']:.2f}",
            va="center", fontsize=7)
ax.set_xscale("log")
ax.set_yticks(ys)
ax.set_yticklabels([labels[s] for s in order], fontsize=8)
ax.axvline(1, color="k", lw=0.8, ls="--")
ax.set_xlabel("Odds ratio for Asn top-decile membership (95% CI, log scale)")
ax.set_title("Fig4a. Regulator-set Asn enrichment: InterPro-formal vs regex baseline")
from matplotlib.lines import Line2D
ax.get_legend = lambda *a, **k: None  # no box: row labels carry formal/regex/GO meaning
fig.text(0.01, 0.01, "*PUF n=2 undetermined. sexual-formal==AP2: no double-count. " +
         "Green triangle = GO-curated TF (n=10): null, underpowered gap, not refutation.", fontsize=6.5)
fig.tight_layout()
fig.savefig(f"{MS}/Fig4a.png", dpi=300)

# ---------------- Fig4b counterfactual strata ----------------
st = pd.read_csv(f"{A}/M4RA_counterfactual_stratified.tsv", sep="\t")
st.to_csv(f"{MS}/panels/fig4b_strata.tsv", sep="\t", index=False)
lg = pd.read_csv(f"{A}/M4RA_counterfactual_logit.tsv", sep="\t")
is_reg = lg[lg.iloc[:, 0] == "is_reg"].iloc[0]
fig, ax = plt.subplots(figsize=(6.6, 3.6))
x = np.arange(len(st))
w = 0.35
ax.bar(x - w/2, st["med_reg"], w, label="regulator", color=C_FORMAL)
ax.bar(x + w/2, st["med_rest"], w, label="rest", color="#999999")
for i, r in st.iterrows():
    ax.text(i, max(r["med_reg"], r["med_rest"]) + 0.004, f"MW p={r['MW_p']:.2g}",
            ha="center", fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels([f"{t} (n_reg={n})" for t, n in zip(st["lcr_tertile"], st["n_reg"])])
ax.set_ylabel("Median Asn fraction")
ax.set_title("Fig4b. Asn excess holds only in high-LCR stratum")
ax.legend(fontsize=8)
fig.text(0.01, 0.01, f"Length/LCR-adjusted logit OR={is_reg['OR']:.2f} "
         f"({is_reg['OR_lo']:.2f}-{is_reg['OR_hi']:.2f}), p={is_reg['p']:.1e}. "
         "Counterfactual is partial, not universal.", fontsize=7)
fig.tight_layout()
fig.savefig(f"{MS}/Fig4b.png", dpi=300)

# ---------------- Fig4c: GCN5 + Pb map + tripole box ----------------
# GCN5: parse day rows from ragged source-data snapshot; condition blocks per merged header
days, neg, pos = [], [], []
with open(f"{A}/M4RA_gcn5_fig2_snapshot.tsv") as h:
    lines = [l.rstrip("\n").split("\t") for l in h]
datarows = [l for l in lines if l and l[0].strip().isdigit()]
for l in datarows:
    vals = [float(v) for v in l[1:13] if v.strip() not in ("", "-")]
    if len(vals) >= 12:
        days.append(int(l[0])); neg.append(vals[:6]); pos.append(vals[6:12])
days = np.array(days); neg = np.array(neg); pos = np.array(pos)
pd.DataFrame({"day": days, "neg_mean": neg.mean(1), "neg_sd": neg.std(1),
              "pos_mean": pos.mean(1), "pos_sd": pos.std(1)}).to_csv(
    f"{MS}/panels/fig4c_gcn5_curves.tsv", sep="\t", index=False)

pb = pd.read_csv(f"{A}/M4RA_pbapiap2_family_map.tsv", sep="\t")
pbv = pb[pb["ko_verified"] == True]
n_blk = int((pbv["transmission_blocked"] == True).sum())
pbv[["mutant", "pf_ortholog", "transmission_blocked"]].to_csv(
    f"{MS}/panels/fig4c_pbmap.tsv", sep="\t", index=False)
tri = pd.read_csv(f"{A}/M4RA_lifecycle_tripole.tsv", sep="\t")
tribox = tri[(tri["test"].str.contains("spearman")) | (tri["contrast"] == "GFP-vs-NoGFP(all-days)")]
tribox.to_csv(f"{MS}/panels/fig4c_tripole_box.tsv", sep="\t", index=False)

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4),
                         gridspec_kw={"width_ratios": [1.2, 0.7, 1.1]})
ax = axes[0]
ax.errorbar(days, neg.mean(1), yerr=neg.std(1), fmt="o-", color=C_FORMAL,
            ecolor=C_FORMAL, capsize=3, label="-RAP (n=6 pts)")
ax.errorbar(days, pos.mean(1), yerr=pos.std(1), fmt="s-", color=C_REGEX,
            ecolor=C_REGEX, capsize=3, label="+RAP (n=6 pts)")
ax.set_xlabel("Day"); ax.set_ylabel("Parasitaemia (%)")
ax.set_title("GCN5 single-locus\nperturbation (source data)")
ax.legend(fontsize=7)
d6r = neg[-1].mean() / pos[-1].mean()
ax.text(0.95, 0.12, f"day-6 ratio ~{d6r:.0f}x", transform=ax.transAxes, fontsize=8,
        ha="right", va="bottom")

ax = axes[1]
ax.bar(["blocked", "not blocked"], [n_blk, len(pbv) - n_blk],
       color=[C_FORMAL, "#999999"])
ax.set_ylabel("# KO-verified Pb ApiAP2 mutants")
ax.set_title("Pb ApiAP2 family map\n(descriptive, n=11)")

ax = axes[2]
ax.axis("off")
box = ("State-consequence tripole (honest box)\n"
       "gametocyte: weak dir., MW q=0.081 n.s. (exploratory)\n"
       "zygote: null (best q=0.073)\n"
       "liver: EXCLUDED \u2014 detection-driven\n"
       "  (top10 detect 7.3% vs 3.4%; controlled n=84 null)\n"
       "=> no shortlist; H-A stays hypothesis\n"
       "L-010 GO check: LANDED\n(curated TF null, n=10,\nunderpowered; T01 superseded)")
ax.text(0.02, 0.98, box, va="top", ha="left", fontsize=8,
        bbox=dict(boxstyle="round", fc="#FFF8DC", ec="grey"))
axes[2].set_title("Tripole verdict")
fig.suptitle("Fig4c. Single-locus proof + family map + negative-report box")
fig.tight_layout()
fig.savefig(f"{MS}/Fig4c.png", dpi=300)
print("Fig4 done. d6 ratio:", round(float(d6r), 1), "| pb n:", len(pbv), "| l010:", os.path.exists(f"{A}/M4RA_l010_enrichment.tsv"))
