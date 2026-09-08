# FIG_FIXLOG — multimodal visual review fixes (8 items)

Scope: scripts under `manuscript/figures/scripts/` only; PNGs re-rendered in place; `data/derived` untouched. Env: `~/.conda/envs/sc/bin/python` + LD_PRELOAD libstdc++. Every fixed PNG read back visually.

## (1) Fig1c + Fig2a truncated titles — FIXED, verified
- mk_fig1.py Fig1c title → two lines; mk_fig2.py Fig2a title → shortened two-line, fontsize 9.
- y-label spelling checked in mk_fig2.py: "synonymous codon usage" correct.
- Verified: Fig1c.png (two-line title fits), Fig2a.png (title fits).

## (2) Fig1a nuclear-track claim without nuclear marks — FIXED, verified
- Compartment table HAS nuclear GC for all 18 spp (fraction scale, 0.18–0.52). Added orange-diamond nuclear markers + legend entry; panel TSV gains `nuclear_gc` column.
- Verified: Fig1a.png shows diamonds + legend; title claim now evidenced.

## (3) Fig3a mixed-unit y-axis — FIXED (redesigned), verified
- Rebuilt as 1×4 faceted bars (LCR fraction / Asn fraction / poly-Asn max run / #InterPro domains), each with independent y-scale + MW p label; sampling-frame box kept as 5th column.
- Verified: Fig3a.png readable, no unit mixing.

## (4) Fig3c left empty category + wrong x-label — FIXED, verified
- Left: single "segregating 34/34" bar + "0 fixed" annotation (dropped empty tick).
- Right: bars recolored by call type (purple amplification / RED HRP2-HRP3 deletions); x-label now "alteration frequency (%) — purple: amplification; red: deletion (HRP2/3)".
- Verified: Fig3c.png.

## (5) Fig4a legend overlap — FIXED, verified
- Legend moved lower-right → upper-left (empty region). Verified: ApiAP2-formal point + q=1.2e-14 label clear.

## (6) Fig4c annotation overlap — FIXED, verified
- "day-6 ratio ~17x" moved to lower-right empty area (ha=right). Verified: legend + annotation both legible.

## (7) Fig5a 50-row forest too dense — FIXED (restructured), verified
- Main Fig5a.png now 7-row arm summary (median rho [min-max] + n/+sig counts; 16/27 + weak-stage + length notes kept).
- Full 50-contrast forest → NEW manuscript/figures/Fig5a_full_forest.png (Extended) + panels/fig5a_full_forest.tsv; main legend references Extended.
- Footnote rewritten short (fits); legend loc="best". Verified both PNGs.

## (8) Fig5b p=0.0 bug — FIXED, verified
- Now reads "p<1e-300, underflow to 0" (JSON value exactly 0.0); figsize height 3.2→2.4. Verified: Fig5b.png.

## Residual notes
- Fig4a/Fig4c footnotes still say "L-010 GO check pending" — correct until L-010 lands; L-010 row goes into Fig4a only when M4RA_l010_enrichment.tsv exists.
- Fig3c deletion red (#D55E00) vs amplification purple: distinguishable + xlabel text backup (not color-only).
- No data files modified; rerun logs deterministic (same inputs → same numbers; Fig5 rows=48 pos-sig=28).
