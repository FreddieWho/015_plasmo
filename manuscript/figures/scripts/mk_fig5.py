#!/usr/bin/env python
"""MS Fig5 (Route C, second context): (a) acute-vs-chronic partition forest;
(b) dTE rule-2 box + M4R-X01 gap box. Plot-only from in-hand TSVs/JSON.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/huyudi/015_plasmo"
C = f"{ROOT}/data/derived/WP4R/M4R-C"
MS = f"{ROOT}/manuscript/figures"
import os
os.makedirs(f"{MS}/scripts", exist_ok=True)
os.makedirs(f"{MS}/panels", exist_ok=True)
os.makedirs(f"{MS}/legends", exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
C_ACUTE, C_CHRONIC, C_OTHER = "#0072B2", "#D55E00", "#999999"

part = pd.read_csv(f"{C}/M4RC_acute_chronic_partition.tsv", sep="\t")
part = part[part["feature"] == "asn_frac"].copy()
val = pd.read_csv(f"{C}/M4RC_validation_table.tsv", sep="\t")[["contrast", "spearman_q"]]
part = part.merge(val, on="contrast", how="left")

def base_key(c):
    return c.replace("_noSchiz", "").replace("_STAGEADJ", "")

# Fig5a rows: STAGEADJ preferred for GSE151189 (honest estimand), base rows otherwise
rows = []
seen = set()
for _, r in part.iterrows():
    c = r["contrast"]
    if "noSchiz" in c:
        continue
    if c.startswith("GSE151189") and not c.endswith("STAGEADJ"):
        continue
    rows.append(r)
sel = pd.DataFrame(rows).copy()
# Fisher-z 95% CI for rho
z = np.arctanh(np.clip(sel["spearman_rho"], -0.999, 0.999))
se = 1.0 / np.sqrt(sel["n"] - 3)
sel["rho_lo"] = np.tanh(z - 1.96 * se)
sel["rho_hi"] = np.tanh(z + 1.96 * se)
arm_order = {"acute": 0, "chronic": 1, "other": 2}
sel["arm_rank"] = sel["arm"].map(arm_order)
sel = sel.sort_values(["arm_rank", "contrast"]).reset_index(drop=True)
sel.to_csv(f"{MS}/panels/fig5a_partition.tsv", sep="\t", index=False)
n_sig = int(((sel["spearman_q"] < 0.05) & (sel["spearman_rho"] > 0)).sum())

def short(c):
    c = c.replace("GSE151189_DHA-CTRL_", "").replace("_STAGEADJ", " [stg-adj]")
    c = c.replace("scRNA_", "").replace("_DHAvsDMSO", "")
    c = c.replace("GSE225340_", "").replace("Mok2021_", "").replace("GSE59099_", "")
    c = c.replace("GSE226632_dTEinteraction", "dTE interaction")
    return c

# ---- Fig5a MAIN: summary partition (one row per arm; full forest -> Extended) ----
def arm_group(c):
    if c.startswith("GSE151189"):
        return "acute microarray (GSE151189)"
    if c.startswith("scRNA_") and "DiD" not in c:
        return "acute scRNA (K13 Seq-Well)"
    if c.startswith("scRNA_DiD"):
        return "acute scRNA DiD (genotype x drug)"
    if c.startswith("GSE225340"):
        return "persistence (dormancy)"
    if c.startswith("Mok2021"):
        return "chronic protein (Mok 2021)"
    if c.startswith("GSE59099"):
        return "chronic transcript (GSE59099)"
    if "dTE" in c:
        return "dTE interaction"
    return "other"
sel["group"] = sel["contrast"].map(arm_group)
order_groups = ["acute microarray (GSE151189)", "acute scRNA (K13 Seq-Well)",
                "acute scRNA DiD (genotype x drug)", "persistence (dormancy)",
                "chronic protein (Mok 2021)", "chronic transcript (GSE59099)",
                "dTE interaction"]
order_groups = [g for g in order_groups if (sel["group"] == g).any()]
exp_dir = {g: 1 for g in order_groups}
exp_dir["chronic protein (Mok 2021)"] = -1  # partition hypothesis: chronic-protein signal is negative
def dir_sig(sub):
    d = exp_dir[sub.name]
    return int((((sel.loc[sub.index, "spearman_rho"] * d) > 0) & (sel.loc[sub.index, "spearman_q"] < 0.05)).sum())
summary = sel.groupby("group").agg(med_rho=("spearman_rho", "median"),
                                     min_rho=("spearman_rho", "min"),
                                     max_rho=("spearman_rho", "max"),
                                     n=("spearman_rho", "size"),
                                     n_dir_sig=("spearman_q", dir_sig))
summary = summary.reindex(order_groups)
summary.to_csv(f"{MS}/panels/fig5a_summary.tsv", sep="\t")
group_arm = {"acute microarray (GSE151189)": "acute", "acute scRNA (K13 Seq-Well)": "acute",
             "acute scRNA DiD (genotype x drug)": "acute", "persistence (dormancy)": "acute",
             "chronic protein (Mok 2021)": "chronic", "chronic transcript (GSE59099)": "chronic",
             "dTE interaction": "other"}
fig, ax = plt.subplots(figsize=(8.4, 3.6))
ys = np.arange(len(summary))
for i, (g, r) in enumerate(summary.iterrows()):
    arm = group_arm[g]
    col = C_ACUTE if arm == "acute" else (C_CHRONIC if arm == "chronic" else C_OTHER)
    mk = "o" if arm == "acute" else ("s" if arm == "chronic" else "^")
    ax.errorbar(r["med_rho"], ys[i],
                xerr=[[r["med_rho"] - r["min_rho"]], [r["max_rho"] - r["med_rho"]]],
                fmt=mk, color=col, ecolor=col, capsize=3, ms=5)
    ax.text(r["max_rho"] + 0.006, ys[i],
            f"med {r['med_rho']:+.2f} (n={int(r['n'])}, dir-sig {int(r['n_dir_sig'])})",
            va="center", fontsize=7)
ax.set_yticks(ys)
ax.set_yticklabels(list(summary.index), fontsize=8)
ax.axvline(0, color="k", lw=0.8, ls="--")
ax.set_xlabel("Spearman rho median [min-max] per arm (Asn fraction vs response)")
ax.set_title("Fig5a. Acute-vs-chronic partition, summary (full 50-contrast forest: Extended Data)")
from matplotlib.lines import Line2D
ax.legend([Line2D([0], [0], marker="o", color=C_ACUTE, ls=""),
           Line2D([0], [0], marker="s", color=C_CHRONIC, ls="")],
          ["acute", "chronic"], fontsize=8, loc="best")
fig.subplots_adjust(bottom=0.24, left=0.28, right=0.97, top=0.90)
ax.text(0.5, -0.52,
         "16/27 survive stage-adj & positive; weak stage ref (R2~1e-4-3%); length atten. large but nonzero; dir-sig = expected-dir q<0.05.",
         fontsize=6.5, ha="center", va="top", color="#333333", transform=ax.transAxes)
fig.tight_layout()
fig.savefig(f"{MS}/Fig5a.png", dpi=300)
# ---- Fig5a FULL forest -> Extended Data ----
fig2, ax2 = plt.subplots(figsize=(8.4, max(5, 0.32 * len(sel) + 1.5)))
ys2 = np.arange(len(sel))
cols2 = sel["arm"].map({"acute": C_ACUTE, "chronic": C_CHRONIC}).fillna(C_OTHER)
mks2 = sel["arm"].map({"acute": "o", "chronic": "s"}).fillna("^")
for i, r in sel.iterrows():
    ax2.errorbar(r["spearman_rho"], ys2[i],
                 xerr=[[r["spearman_rho"] - r["rho_lo"]], [r["rho_hi"] - r["spearman_rho"]]],
                 fmt=mks2[i], color=cols2[i], ecolor=cols2[i], capsize=2, ms=4)
    if r["spearman_q"] < 0.05:
        ax2.text(r["rho_hi"] + 0.004, ys2[i], "*", va="center", fontsize=9)
ax2.set_yticks(ys2)
ax2.set_yticklabels([short(c) for c in sel["contrast"]], fontsize=6.5)
ax2.axvline(0, color="k", lw=0.8, ls="--")
ax2.set_xlabel("Spearman rho (Asn fraction vs response, 95% CI)  * q<0.05")
ax2.set_title("Extended Data: full acute-vs-chronic forest (all 50 STAGEADJ-preferred contrasts)")
fig2.tight_layout()
fig2.savefig(f"{MS}/Fig5a_full_forest.png", dpi=300)
sel.to_csv(f"{MS}/panels/fig5a_full_forest.tsv", sep="\t", index=False)

# ---------------- Fig5b boxes ----------------
dte = json.load(open(f"{C}/M4RC_dTE_upgrade.json"))
fig, axes = plt.subplots(1, 2, figsize=(9.5, 2.4))
ax = axes[0]
ax.axis("off")
ax.text(0.02, 0.95,
        "dTE interaction upgrade (rule 2)\n"
        f"sign agreement ratio-vs-interaction: {dte['sign_agreement']:.3f}\n"
        f"spearman: {dte['spearman']:.3f} (p" +
        (f"={dte['spearman_p']:.2g}" if dte['spearman_p'] > 0 else "<1e-300, underflow to 0") + ")\n"
        f"q<0.05 genes: {dte['n_q05']}  -> WEAK SUPPORT retained\n"
        "ratio-based dTE stays DISCOVERY/SUPPORTING",
        va="top", ha="left", fontsize=9,
        bbox=dict(boxstyle="round", fc="#EAF2FF", ec="grey"))
ax.set_title("dTE")
ax = axes[1]
ax.axis("off")
ax.text(0.02, 0.95,
        "M4R-X01 gap (standing)\n"
        "No suitable public acute DHA x K13\n"
        "translation-layer data (Ribo/polysome/\n"
        "proteome/phospho/charging). scRNA is\n"
        "transcript-only and does NOT close it.\n"
        "C stays partition figure, not MOA.",
        va="top", ha="left", fontsize=9,
        bbox=dict(boxstyle="round", fc="#FFF8DC", ec="grey"))
ax.set_title("Gap")
fig.suptitle("Fig5b. dTE rule-2 + standing translation-layer gap")
fig.tight_layout()
fig.savefig(f"{MS}/Fig5b.png", dpi=300)
print(f"Fig5 done. rows={len(sel)} pos-sig={n_sig} dTE_agree={dte['sign_agreement']:.3f}")
