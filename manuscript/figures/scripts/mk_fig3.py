#!/usr/bin/env python3
"""Fig3 panels (a,b,c). Plot only from in-hand TSVs. No new analysis.
Reads: M4RD_D2b_matched.tsv, M4RD_D2b_framebias.tsv, M4RD_D2b_verdict.txt,
       M4RD_D6formal_summary.tsv, M4RD_D6formal_pertransition.tsv,
       M4RD_D4_Pf8boundary.tsv, M4RD_D4_Pf8cnv.tsv
Writes: manuscript/figures/Fig3{a,b,c}.png (300dpi) + panels/Fig3{a,b,c}.tsv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/huyudi/015_plasmo"
M4D = f"{ROOT}/data/derived/WP4R/M4R-D"
OUT = f"{ROOT}/manuscript/figures"
os.makedirs(f"{OUT}/panels", exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, RED, PURPLE = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"

verdict = open(f"{M4D}/M4RD_D2b_verdict.txt").read().strip()

# ---------- Fig3a: D2b harbor vs matched ----------
m = pd.read_csv(f"{M4D}/M4RD_D2b_matched.tsv", sep="\t")
fb = pd.read_csv(f"{M4D}/M4RD_D2b_framebias.tsv", sep="\t").set_index("feature")
order = ["lcr_frac", "asn_frac", "polyN_max", "n_ipr"]
pretty = {"lcr_frac": "LCR fraction", "asn_frac": "Asn fraction",
          "polyN_max": "poly-Asn max run", "n_ipr": "#InterPro domains"}
fig = plt.figure(figsize=(11, 3.8))
gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 1, 1.6])
mm = m.set_index("feature")
for j, f in enumerate(order):
    ax = fig.add_subplot(gs[0, j])
    hv, bv = mm.loc[f, "median_harbor"], mm.loc[f, "median_bg"]
    ax.bar([0, 1], [hv, bv], color=[BLUE, ORANGE], edgecolor="black", linewidth=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["harbor\n(n=91)", "matched\nbg"], fontsize=7)
    ax.set_title(pretty[f], fontsize=8)
    top = max(hv, bv)
    ax.set_ylim(0, top * 1.35 if top > 0 else 1)
    ax.text(0.5, top * 1.12, f"p={mm.loc[f, 'MW_p']:.1g}", ha="center", fontsize=7,
            transform=ax.transData)
    if j == 0:
        ax.set_ylabel("median (per-feature scale)")
fig.suptitle("(a) Harbor genes: no LCR excess; Asn-depleted; domain-rich (independent scales per feature)",
             fontsize=10, x=0.02, ha="left")
ax = fig.add_subplot(gs[0, 4])
ax.axis("off")
ax.text(0.02, 0.95,
        "Sampling-frame bias (core-181 vs rest of proteome):\n"
        f"length median {fb.loc['length','median_core']:.0f} vs {fb.loc['length','median_noncore']:.0f} aa "
        f"(p={fb.loc['length','MW_p']:.1e})\n"
        f"Asn median {fb.loc['asn_frac','median_core']:.3f} vs "
        f"{fb.loc['asn_frac','median_noncore']:.3f} (p={fb.loc['asn_frac','MW_p']:.1e})\n"
        f"D2b verdict: {verdict} (frame-internal)\n"
        "Claim ceiling: CANDIDATE/C2; no repeat-first language.",
        fontsize=8, va="top", ha="left",
        bbox=dict(boxstyle="round", facecolor="#f0f0f0", edgecolor="grey"))
fig.tight_layout()
fig.savefig(f"{OUT}/Fig3a.png", dpi=300)
m.to_csv(f"{OUT}/panels/Fig3a.tsv", sep="\t", index=False)
print(f"Fig3a: verdict={verdict} harbor_n=91")

# ---------- Fig3b: D6 transitions ----------
s = pd.read_csv(f"{M4D}/M4RD_D6formal_summary.tsv", sep="\t")
pt = pd.read_csv(f"{M4D}/M4RD_D6formal_pertransition.tsv", sep="\t")
trans = ["T1_AT_Lav", "T2_GC_vivax", "T3_rodent"]
tlabel = {"T1_AT_Lav": "T1 AT (Laverania)\nregistered NEGATIVE",
          "T2_GC_vivax": "T2 GC (vivax clade)",
          "T3_rodent": "T3 rodent"}
sets = ["AP2", "PUF_RNA", "CHROM"]
markers = {"AP2": "o", "PUF_RNA": "s", "CHROM": "^"}
colors = {"AP2": BLUE, "PUF_RNA": ORANGE, "CHROM": GREEN}
fig, ax = plt.subplots(figsize=(8.5, 4.0))
x = np.arange(len(trans))
for j, st in enumerate(sets):
    sub = s[s.set == st].set_index("transition")
    fracs = [sub.loc[t, "n_expected_sign"] / sub.loc[t, "n_genes"] for t in trans]
    ax.scatter(x + (j - 1) * 0.14, fracs, s=90, marker=markers[st], color=colors[st],
               edgecolors="black", linewidths=0.6, label=st, zorder=3)
    for i, t in enumerate(trans):
        r = sub.loc[t]
        ax.text(x[i] + (j - 1) * 0.14 + 0.05, fracs[i] + 0.02,
                f'{int(r.n_expected_sign)}/{int(r.n_genes)} p={r.signtest_p:.3g}', fontsize=6.5)
ax.axhline(0.5, color="grey", ls="--", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels([tlabel[t] for t in trans], fontsize=8)
ax.set_ylabel("fraction with expected sign")
ax.set_ylim(0, 1.12)
ax.set_title("(b) Natural experiment: T2+T3 composition-consistent (rule met); "
             "T1 negative kept; CHROM null throughout", fontsize=10, loc="left")
ax.legend(fontsize=8, title="gene set")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig3b.png", dpi=300)
s.to_csv(f"{OUT}/panels/Fig3b.tsv", sep="\t", index=False)
print(f"Fig3b: T2_AP2={s[(s.transition=='T2_GC_vivax')&(s.set=='AP2')].iloc[0].n_expected_sign}/4 "
      f"T3_PUF p={s[(s.transition=='T3_rodent')&(s.set=='PUF_RNA')].iloc[0].signtest_p:.3g}")

# ---------- Fig3c: Pf8 population layer ----------
bnd = pd.read_csv(f"{M4D}/M4RD_D4_Pf8boundary.tsv", sep="\t")
cnv = pd.read_csv(f"{M4D}/M4RD_D4_Pf8cnv.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(10, 4.0), gridspec_kw={"width_ratios": [3, 2]})
ax = axes[0]
seg = bnd.segregating.fillna(False).astype(bool) if bnd.segregating.dtype == object else bnd.segregating
n_seg, n_tot = int(seg.sum()), len(bnd)
ax.bar([f"segregating\n{n_seg}/{n_tot}"], [n_seg], color=GREEN,
       edgecolor="black", linewidth=0.6)
ax.set_xlim(-0.6, 0.6)
ax.set_ylabel("# resistance markers (Pf8)")
ax.text(0, n_seg + 0.4, f"{n_seg}/{n_tot} segregate; 0 fixed", ha="center", fontsize=8)
ax.set_title(f"(c) {n_seg}/{n_tot} known resistance markers segregate in present-day populations",
             fontsize=10, loc="left")
ax = axes[1]
ax.barh(np.arange(len(cnv)), cnv.freq * 100, color=PURPLE, edgecolor="black", linewidth=0.5)
ax.set_yticks(np.arange(len(cnv)))
ax.set_yticklabels([c.replace("_final_amplification_call", "").replace("_final_deletion_call", "") for c in cnv.locus_call], fontsize=8)
is_del = cnv.locus_call.str.contains("deletion")
for patch, dd in zip(ax.patches, is_del):
    if dd:
        patch.set_facecolor(RED)
ax.set_xlabel("alteration frequency (%) — purple: amplification; red: deletion (HRP2/3)")
for i, r in cnv.iterrows():
    ax.text(r.freq * 100 + 0.4, i, f"{r.freq*100:.1f}% (n={int(r.n_amp)})", va="center", fontsize=7.5)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig3c.png", dpi=300)
bnd.to_csv(f"{OUT}/panels/Fig3c_markers.tsv", sep="\t", index=False)
cnv.to_csv(f"{OUT}/panels/Fig3c_cnv.tsv", sep="\t", index=False)
print(f"Fig3c: segregating={n_seg}/{n_tot}")
print("DONE Fig3")
