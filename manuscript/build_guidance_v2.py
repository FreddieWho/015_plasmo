#!/usr/bin/env python3
"""Assemble single-file manuscript HTML: draft text + figures (embedded) + legends + tables.
Usage: build_manuscript_html.py  -> manuscript/PLASMODIUM-C2F_guidance_doc_v2.html
"""
import base64, html, os, subprocess

ROOT = "/home/huyudi/015_plasmo"
MS = f"{ROOT}/manuscript"
FIGS = ["Fig1a", "Fig1b", "Fig1c", "Fig2a", "Fig2b", "Fig2c",
        "Fig3a", "Fig3b", "Fig3c", "Fig4a", "Fig4b", "Fig4c",
        "Fig5a", "Fig5b", "Fig5a_full_forest"]
LEGEND = {"Fig1": "legends/Fig1_legend.txt", "Fig2": "legends/Fig2_legend.txt",
          "Fig3": "legends/Fig3_legend.txt", "Fig4": "legends/Fig4_legend.txt",
          "Fig5": "legends/Fig5_legend.txt"}
TABLES = [("TableS1_D1_concordance.tsv", "D1 branch-aware concordance"),
          ("TableS1b_D1_LOO.tsv", "D1 branch leave-one-out"),
          ("TableS2_D2b_matched.tsv", "D2b harbor vs matched background"),
          ("TableS2b_D2b_framebias.tsv", "D2b sampling-frame bias quant"),
          ("TableS3_D6formal_summary.tsv", "D6 formal summary"),
          ("TableS3b_D6formal_pertransition.tsv", "D6 per-transition tests"),
          ("TableS3c_D3_multidim.tsv", "D3 multidim matched controls"),
          ("TableS3d_Pf8boundary.tsv", "Pf8 population boundary"),
          ("TableS4_interpro_enrichment.tsv", "InterPro-formal vs regex enrichment"),
          ("TableS5_lifecycle_tripole.tsv", "Lifecycle tripole"),
          ("TableS6_C_partition_asnfrac.tsv", "C acute-vs-chronic partition")]


def tsv_to_md(path, max_rows=250):
    lines = open(path).read().rstrip("\n").split("\n")
    head = lines[0].split("\t")
    out = ["| " + " | ".join(head) + " |", "| " + " | ".join(["---"] * len(head)) + " |"]
    rows = lines[1:]
    trunc = len(rows) > max_rows
    for r in rows[:max_rows]:
        cells = [c.replace("|", "\\|") for c in r.split("\t")]
        out.append("| " + " | ".join(cells) + " |")
    if trunc:
        out.append(f"\n*Table truncated to {max_rows} rows; full table in `{os.path.basename(path)}`.*")
    return "\n".join(out)


parts = []
draft = open(f"{MS}/text/guidance_doc_v02.md").read()
parts.append(draft)
parts.append("\n\n# Figures\n")
for f in FIGS:
    p = f"{MS}/figures/{f}.png"
    assert os.path.exists(p), f"missing {p}"
    parts.append(f"\n## {f}\n\n![]({p})\n")
    key = f[:4]
    leg = f"{MS}/figures/{LEGEND[key]}"
    if os.path.exists(leg):
        parts.append("\n**Legend.** " + open(leg).read().strip().replace("\n", " ") + "\n")
parts.append("\n\n# Supplementary Tables\n")
for fname, title in TABLES:
    p = f"{MS}/tables/{fname}"
    assert os.path.exists(p), f"missing {p}"
    parts.append(f"\n## {fname} — {title}\n\n" + tsv_to_md(p) + "\n")

combined = f"{MS}/_combined_guidance_v2.md"
open(combined, "w").write("\n".join(parts))
out = f"{MS}/PLASMODIUM-C2F_guidance_doc_v2.html"
r = subprocess.run(["pandoc", "--standalone", "--self-contained", "--metadata",
                    "title=PLASMODIUM-C2F guidance doc v2",
                    "--to", "html5", combined, "-o", out],
                   capture_output=True, text=True, timeout=300)
print("pandoc rc:", r.returncode, r.stderr[:500])
h = open(out, encoding="utf-8").read()
n_img = h.count("data:image/png;base64")
print(f"html bytes={len(h)} embedded_png={n_img} (expect {len(FIGS)})")
assert n_img == len(FIGS), "image count mismatch"
print("BUILD_OK:", out)
