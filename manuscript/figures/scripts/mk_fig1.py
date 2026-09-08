#!/usr/bin/env python3
"""Fig1 panels (a,b,c). Plot only from in-hand TSVs. No new analysis.
Reads: M1-01_composition_table.tsv, M1-01_compartment_breakdown.tsv,
       M1-02_species_tree.nwk, M1-02_ancestral_gc.tsv,
       M1-02_independence_grade.tsv, M1-02_leave_one_clade_out.tsv
Writes: manuscript/figures/Fig1{a,b,c}.png (300dpi) + panels/Fig1{a,b,c}.tsv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT = "/home/huyudi/015_plasmo"
M1C = f"{ROOT}/data/derived/WP1/M1-01_composition"
M1P = f"{ROOT}/data/derived/WP1/M1-02_phylogeny"
OUT = f"{ROOT}/manuscript/figures"
os.makedirs(f"{OUT}/panels", exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, RED, PURPLE = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"

SHORT = {"Plasmodium": "P.", "Toxoplasma": "T.", "Babesia": "B.", "Theileria": "Th.",
         "Hepatocystis": "H.", "Haemoproteus": "Hm.", "Leucocytozoon": "L."}

def short(s):
    for g, ab in SHORT.items():
        if s.startswith(g):
            return s.replace(g, ab, 1)
    return s

# ---------- Fig1a: genome GC bar + compartment layering ----------
comp = pd.read_csv(f"{M1C}/M1-01_composition_table.tsv", sep="\t")
bd = pd.read_csv(f"{M1C}/M1-01_compartment_breakdown.tsv", sep="\t")
comp = comp.sort_values("genome_gc").reset_index(drop=True)
nuc = bd[bd.compartment == "nuclear_chromosome"].set_index("species_id")["gc"]

fig, ax = plt.subplots(figsize=(9, 5.2))
y = np.arange(len(comp))
ax.barh(y, comp.genome_gc * 100, color=BLUE, edgecolor="black", linewidth=0.4, label="genome GC%")
nucv = np.array([nuc.get(s, np.nan) for s in comp.species_id]) * 100
ax.scatter(nucv, y, s=42, marker="D", color=ORANGE, edgecolors="black",
           linewidths=0.5, zorder=3, label="nuclear-chromosome GC%")
ax.set_yticks(y)
ax.set_yticklabels([short(s) for s in comp.species], fontsize=8)
ax.set_xlabel("Genome GC (%)")
ax.set_title("(a) Genome GC spans 18.2-52.3% across 18 species; nuclear-chromosome GC tracks genome GC",
             fontsize=10, loc="left")
for i, r in comp.iterrows():
    ax.text(r.genome_gc * 100 + 0.4, i, f"{r.genome_gc*100:.1f}", va="center", fontsize=7)
ax.set_xlim(0, comp.genome_gc.max() * 100 + 8)
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig1a.png", dpi=300)
comp["nuclear_gc"] = [nuc.get(s, np.nan) for s in comp.species_id]
comp[["species_id", "species", "genome_gc", "gc3", "nuclear_chromosome_len", "nuclear_gc"]].to_csv(
    f"{OUT}/panels/Fig1a.tsv", sep="\t", index=False)
print("Fig1a: n_species =", len(comp), "gc_range =",
      round(comp.genome_gc.min()*100, 1), "-", round(comp.genome_gc.max()*100, 1))

# ---------- Newick parser ----------
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = [0]
    def parse_node():
        children = []
        name = None
        length = 0.0
        if s[pos[0]] == "(":
            pos[0] += 1
            children.append(parse_node())
            while s[pos[0]] == ",":
                pos[0] += 1
                children.append(parse_node())
            assert s[pos[0]] == ")", s[pos[0]]
            pos[0] += 1
        nm = []
        while pos[0] < len(s) and s[pos[0]] not in ",():;":
            nm.append(s[pos[0]]); pos[0] += 1
        nm = "".join(nm)
        if s[pos[0]:pos[0]+1] == ":":
            pos[0] += 1
            ln = []
            while pos[0] < len(s) and s[pos[0]] not in ",();":
                ln.append(s[pos[0]]); pos[0] += 1
            length = float("".join(ln)) if "".join(ln) else 0.0
        return {"name": nm or None, "length": length, "children": children}
    return parse_node()

nwk = open(f"{M1P}/M1-02_species_tree.nwk").read()
tree = parse_newick(nwk)
anc = pd.read_csv(f"{M1P}/M1-02_ancestral_gc.tsv", sep="\t")
anc_map = {}
for _, r in anc.iterrows():
    key = frozenset([x.strip() for x in str(r.descendants).split(",")])
    anc_map[key] = (r.node, r.gc)
gc_tip = comp.set_index("species_id")["genome_gc"].to_dict()
sp_name = comp.set_index("species_id")["species"].apply(short).to_dict()

leaves, internals = [], []
yy = [0]
def layout(node, depth):
    if not node["children"]:
        node["_y"] = float(yy[0]); yy[0] += 1
        node["_x"] = float(depth)
        leaves.append(node)
    else:
        for ch in node["children"]:
            layout(ch, depth + 1)
        node["_y"] = float(np.mean([ch["_y"] for ch in node["children"]]))
        node["_x"] = float(depth)
        internals.append(node)
layout(tree, 0)

def leafset(node):
    if not node["children"]:
        return frozenset([node["name"]])
    out = frozenset()
    for ch in node["children"]:
        out = out | leafset(ch)
    return out

fig, ax = plt.subplots(figsize=(9, 6.5))
def draw(node):
    for ch in node["children"]:
        ax.plot([node["_x"], node["_x"]], [node["_y"], ch["_y"]], color="black", lw=0.8)
        ax.plot([node["_x"], ch["_x"]], [ch["_y"], ch["_y"]], color="black", lw=0.8)
        draw(ch)
draw(tree)
maxx = max(n["_x"] for n in leaves)
for n in leaves:
    g = gc_tip.get(n["name"], np.nan)
    ax.scatter([maxx + 0.6], [n["_y"]], s=60, c=[g], cmap="coolwarm", vmin=0.15, vmax=0.55,
               edgecolors="black", linewidths=0.5, zorder=3, marker="o")
    ax.text(maxx + 0.75, n["_y"], f'{sp_name.get(n["name"], n["name"])} ({n["name"]}) {g*100:.1f}%' if not np.isnan(g) else n["name"],
            va="center", fontsize=7)
for n in internals:
    key = leafset(n)
    if key in anc_map:
        nm, g = anc_map[key]
        ax.text(n["_x"] - 0.05, n["_y"] + 0.12, f"{g*100:.1f}%", ha="right", fontsize=6.5,
                color=RED, style="italic")
ax.set_xlim(-0.3, maxx + 2.6)
ax.set_ylim(-1, len(leaves))
ax.set_xlabel("cladogram depth (branch lengths schematic; k-mer distances saturated)")
ax.set_title("(b) Species tree (k-mer NJ) with tip genome GC (dots) and ancestral GC (red italics)",
             fontsize=10, loc="left")
ax.set_yticks([])
fig.tight_layout()
fig.savefig(f"{OUT}/Fig1b.png", dpi=300)
pd.DataFrame([{"leaf": n["name"], "species": sp_name.get(n["name"]), "genome_gc": gc_tip.get(n["name"])} for n in leaves]).to_csv(
    f"{OUT}/panels/Fig1b.tsv", sep="\t", index=False)
print("Fig1b: tips =", len(leaves), "internal_annotated =",
      sum(1 for n in internals if leafset(n) in anc_map))

# ---------- Fig1c: LOO robustness ----------
loo = pd.read_csv(f"{M1P}/M1-02_leave_one_clade_out.tsv", sep="\t")
gr = pd.read_csv(f"{M1P}/M1-02_independence_grade.tsv", sep="\t", header=None, names=["metric", "value"])
grd = dict(zip(gr.metric, gr.value))
fig, ax = plt.subplots(figsize=(9, 3.2))
y = np.arange(len(loo))
ax.barh(y, loo.gc_range_remaining * 100, color=GREEN, edgecolor="black", linewidth=0.4)
ax.set_yticks(y)
ax.set_yticklabels(loo.removed_clade, fontsize=8)
ax.set_xlabel("Retained GC range after clade removal (pp)")
for i, r in loo.iterrows():
    ax.text(r.gc_range_remaining * 100 + 0.3, i,
            f'{r.gc_range_remaining*100:.1f}pp retain={r.retains_gradient}', va="center", fontsize=7)
ax.set_title(f"(c) Leave-one-clade-out: gradient retained {loo.retains_gradient.sum()}/{len(loo)}\n"
             f"independence grade {grd.get('grade')} (low: {grd.get('independent_low_clades')}; "
             f"high: {grd.get('independent_high_clades')})", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig1c.png", dpi=300)
loo.to_csv(f"{OUT}/panels/Fig1c.tsv", sep="\t", index=False)
print("Fig1c: retains =", int(loo.retains_gradient.sum()), "/", len(loo))
print("DONE Fig1")
