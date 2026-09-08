#!/usr/bin/env python3
"""Fig2 panels (a,b,c). Plot only from in-hand TSVs. No new analysis.
Reads: M2-01_background_explained.tsv, M2-02_background_explained.tsv,
       M4RD_D1_concordance_summary.tsv, M4RD_D1_concordance_aux.tsv, M4RD_D1_LOO.tsv,
       M4RD_D3_matched_controls.tsv, M4RD_D3_multidim.tsv
Writes: manuscript/figures/Fig2{a,b,c}.png (300dpi) + panels/Fig2{a,b,c}.tsv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/huyudi/015_plasmo"
M2C = f"{ROOT}/data/derived/WP2/M2-01_codon"
M2A = f"{ROOT}/data/derived/WP2/M2-02_aa"
M4D = f"{ROOT}/data/derived/WP4R/M4R-D"
OUT = f"{ROOT}/manuscript/figures"
os.makedirs(f"{OUT}/panels", exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, RED = "#0072B2", "#E69F00", "#009E73", "#D55E00"

# ---------- Fig2a: conduction bars ----------
m2c = pd.read_csv(f"{M2C}/M2-01_background_explained.tsv", sep="\t")
m2a = pd.read_csv(f"{M2A}/M2-02_background_explained.tsv", sep="\t")
codon_med = m2c.prop_explained_GC3_masked.median()
aa_med = m2a.r2_masked.median()
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.barh([1, 0], [codon_med * 100, aa_med * 100], color=[BLUE, ORANGE],
        edgecolor="black", linewidth=0.6, height=0.55)
ax.set_yticks([1, 0])
ax.set_yticklabels([f"synonymous codon usage\n(GC3 explains, median, n={len(m2c)} fam-x-spp)",
                    f"amino-acid CLR\n(genome GC explains, median, n={len(m2a)} AAs)"])
ax.set_xlabel("Background-explained variation (%)")
for v, lab in [(codon_med * 100, 1), (aa_med * 100, 0)]:
    ax.text(v + 0.6, lab, f"{v:.1f}%", va="center", fontsize=9, fontweight="bold")
ax.set_xlim(0, 100)
ax.set_title("(a) Composition explains most codon/AA variation\n(residuals are the minority)",
             fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig2a.png", dpi=300)
pd.DataFrame([{"layer": "codon_GC3", "median_explained": codon_med, "n": len(m2c)},
              {"layer": "aa_genomeGC", "median_explained": aa_med, "n": len(m2a)}]).to_csv(
    f"{OUT}/panels/Fig2a.tsv", sep="\t", index=False)
print(f"Fig2a: codon_med={codon_med:.4f} aa_med={aa_med:.4f}")

# ---------- Fig2b: branch-aware forest ----------
cs = pd.read_csv(f"{M4D}/M4RD_D1_concordance_summary.tsv", sep="\t").set_index("stratum")
aux = pd.read_csv(f"{M4D}/M4RD_D1_concordance_aux.tsv", sep="\t")
loo = pd.read_csv(f"{M4D}/M4RD_D1_LOO.tsv", sep="\t")
or_full = float(aux.comparison.map(lambda s: s == "coupled_vs_uncoupled_concordance").pipe(
    lambda m: aux.loc[m, "OR"].iloc[0]))
or_lo, or_hi = loo.Fisher_OR.min(), loo.Fisher_OR.max()
rows = [("coupled sites", "informative_coupled", "o"),
        ("uncoupled sites", "informative_uncoupled", "s"),
        ("all informative branches", "informative_branches", "^")]
fig, ax = plt.subplots(figsize=(7.5, 3.6))
for i, (lab, key, mk) in enumerate(rows):
    r = cs.loc[key]
    ax.errorbar(r.frac, i, xerr=[[r.frac - r.ci_lo], [r.ci_hi - r.frac]],
                fmt=mk, color="black", ecolor="black", elinewidth=1.2, capsize=4, ms=7,
                markerfacecolor=BLUE if "coupled" == key.split("_")[-1] else (ORANGE if "uncoupled" in key else GREEN))
    ax.text(r.ci_hi + 0.012, i, f'{r.frac:.3f} [{r.ci_lo:.3f}-{r.ci_hi:.3f}] n={int(r.n_events)}',
            va="center", fontsize=8)
ax.axvline(0.5, color="grey", ls="--", lw=0.8)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows])
ax.set_xlabel("Branch-concordance fraction (AA replacement follows branch GC/AT change)")
ax.set_title(f"(b) Coupled sites follow branch composition change\nOR={or_full:.2f}; LOO OR range "
             f"{or_lo:.2f}-{or_hi:.2f}, 7/7 significant", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig2b.png", dpi=300)
cs.reset_index()[["stratum", "n_events", "n_concordant", "frac", "ci_lo", "ci_hi", "binom_p"]].to_csv(
    f"{OUT}/panels/Fig2b.tsv", sep="\t", index=False)
print(f"Fig2b: OR={or_full:.3f} LOO range={or_lo:.2f}-{or_hi:.2f}")

# ---------- Fig2c: adaptive-site separation ----------
mc = pd.read_csv(f"{M4D}/M4RD_D3_matched_controls.tsv", sep="\t")
md = pd.read_csv(f"{M4D}/M4RD_D3_multidim.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(7.5, 4.2))
x = np.arange(len(mc))
ax.bar(x, mc.frac_bg_coupled * 100, color=ORANGE, edgecolor="black", linewidth=0.5,
       label="matched-background coupled fraction")
ax.scatter(x, np.zeros(len(mc)), s=70, marker="x", color=RED, linewidths=1.6,
           label="known resistance site: uncoupled (14/14)", zorder=3)
ax.set_xticks(x)
labels = [r.gene.split("_")[0] + f" {r.pf_pos}\ncons={r.cons_known:.2f}" for _, r in mc.iterrows()]
ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
ax.set_ylabel("Background coupled fraction (%)")
ax.set_title("(c) Known resistance sites sit at conserved, composition-uncoupled positions;\n"
             "matched backgrounds couple only at lower conservation (weak reading: distinct sequence spaces)",
             fontsize=10, loc="left")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig2c.png", dpi=300)
mc.to_csv(f"{OUT}/panels/Fig2c.tsv", sep="\t", index=False)
print(f"Fig2c: known sites={len(mc)} all_uncoupled={bool((mc.known_coupled == False).all())}")
print("DONE Fig2")
