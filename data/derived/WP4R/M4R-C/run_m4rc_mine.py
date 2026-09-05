#!/usr/bin/env python3
"""M4R-C ROUND1 main-agent driver: bounded steps, reuses salvage helpers.
Usage: run_m4rc_mine.py STEP   (STEP in C2 D1 D2 E1 E2 F FINAL)
Each STEP saves state incrementally; safe to re-run (skips computed contrasts).
"""
import sys
import numpy as np
from scipy import stats
import functools as _ft
_orig_mwu = stats.mannwhitneyu
stats.mannwhitneyu = _ft.partial(_orig_mwu, method="asymptotic")
print("MWU patched to asymptotic", flush=True)
src = open("data/derived/WP4R/M4R-C/run_m4rc_salvage.py").read().split("STEP = ")[0]
exec(src)
import time
t0 = time.time()
STEP = sys.argv[1]
ALL = load_state()
print(f"MINE STEP={STEP} loaded_rows={len(ALL)}", flush=True)
have = set(r["contrast"] for r in ALL)

if STEP == "C2":
    RESP = pd.read_csv("data/derived/WP4/L1_codon_x_stress/L1_response_by_gene.tsv", sep="\t")
    sv = pd.read_csv(f"{OUT}/M4RC_stagevec.tsv", sep="\t", index_col=0).iloc[:, 0].dropna()
    n = 0
    for (bg, tm), j in RESP.groupby(["bg", "time"]):
        y = j.set_index("gene")["resp"].dropna()
        y = y[y.index.isin(FEAT.index)]
        if len(y) < 100:
            continue
        if y.nunique() < 10:
            print(f"  skip degenerate {bg}_{tm}h nunique={y.nunique()}", flush=True)
            continue
        lab = f"GSE151189_DHA-CTRL_{bg}_{tm}h"
        if lab not in have:
            ALL += ftests(y, lab); ALL += ftests(no_schiz(y), lab + "_noSchiz"); n += 2
        common = y.index.intersection(sv.index)
        if len(common) > 100:
            yy = y.loc[common].values; xx = sv.loc[common].values
            X = np.column_stack([np.ones(len(common)), xx])
            beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
            resid = pd.Series(yy - X @ beta, index=common)
            if lab + "_STAGEADJ" not in have:
                ALL += ftests(resid, lab + "_STAGEADJ")
                ALL += ftests(no_schiz(resid), lab + "_STAGEADJ_noSchiz"); n += 2
    save_state(ALL); checkpoint("C2_DHA", ALL)
    print(f"MINE C2 DONE new={n} T={round(time.time()-t0,1)}", flush=True)

elif STEP == "D1":
    with gzip.open("data/raw/geo/GSE225340/GSE225340_all_sample_fpkm_table.txt.gz", "rt") as h:
        fpkm = pd.read_csv(h, sep="\t", index_col=0)
    fpkm.index = fpkm.index.str.replace(r"\.\d+$", "", regex=True)
    ring_cols = [c for c in fpkm.columns if "Rings" in c]
    late_cols = [c for c in fpkm.columns if re.search(r"D([5-9]|1[0-2])$", c)]
    early_cols = [c for c in fpkm.columns if re.search(r"D[1-4]$", c)]
    print(f"  ring={len(ring_cols)} early={len(early_cols)} late={len(late_cols)} cols={fpkm.shape}", flush=True)
    R = fpkm[ring_cols]
    n = 0
    for nm, cols in [("DORM_LATE_vs_Rings", late_cols), ("DORM_EARLY_vs_Rings", early_cols)]:
        out = {}
        for g in fpkm.index:
            r = R.loc[g].dropna().values; d = fpkm.loc[g, cols].dropna().values
            if len(r) >= 2 and len(d) >= 3:
                out[g] = float(np.log2(np.median(d) + 1) - np.log2(np.median(r) + 1))
        y = pd.Series(out)
        print(f"  {nm} n={len(y)} med={y.median():.3f}", flush=True)
        if f"GSE225340_{nm}" not in have:
            ALL += ftests(y, f"GSE225340_{nm}")
            ALL += ftests(no_schiz(y), f"GSE225340_{nm}_noSchiz"); n += 2
    save_state(ALL); checkpoint("D1_persist", ALL)
    print(f"MINE D1 DONE new={n} T={round(time.time()-t0,1)}", flush=True)

if STEP == "D2":
    print("D2: K13 scRNA...", flush=True)
    dge_dir = "data/raw/k13scrna/M4R_K13_SCRNA_2026/K13_mt_v_wt_sc_v1/data/raw/dge"
    smap = pd.read_csv("data/raw/k13scrna/M4R_K13_SCRNA_2026/K13_mt_v_wt_sc_v1/data/raw/sample_map.csv")
    pb = {}
    for _, r in smap.iterrows():
        f = f"{dge_dir}/{r['sample']}_gene_exon_dge_n20000.txt.gz"
        try:
            with gzip.open(f, "rt") as h:
                m = pd.read_csv(h, sep="\t", index_col=0)
        except FileNotFoundError:
            print(f"  MISSING {f}", flush=True); continue
        m.index = m.index.str.replace(r"\.\d+$", "", regex=True)
        s = m.sum(axis=1); s = s[s >= 10]
        pb[r["sample"]] = s
        print(f"  {r['sample']}: cells={m.shape[1]} genes_ok={len(s)}", flush=True)
    S = pd.DataFrame(pb).fillna(0)
    print(f"  pseudobulk {S.shape}", flush=True)
    S.to_csv(f"{OUT}/M4RC_scRNA_pseudobulk.tsv", sep="\t")
    cpm = np.log2(S.div(S.sum(axis=0), axis=1) * 1e6 + 1)
    def sc_contrast(a, b, lab):
        global ALL
        if lab in have:
            return
        if a not in cpm.columns or b not in cpm.columns:
            print(f"  SKIP {lab} (missing cols)", flush=True); return
        y = (cpm[a] - cpm[b]).dropna()
        y = y[y.index.isin(FEAT.index)]
        if len(y) < 100 or y.nunique() < 10:
            print(f"  SKIP {lab} degenerate n={len(y)}", flush=True); return
        print(f"  {lab} n={len(y)} med={y.median():.3f}", flush=True)
        ALL += ftests(y, lab)
        ALL += ftests(no_schiz(y), lab + "_noSchiz")
    nn = 0
    for t in ["2h", "4h", "6h"]:
        sc_contrast(f"MRA1250_{t}DHA", f"MRA1250_{t}DMSO", f"scRNA_WT_{t}_DHAvsDMSO")
        sc_contrast(f"MRA1251_{t}DHA", f"MRA1251_{t}DMSO", f"scRNA_580Y_{t}_DHAvsDMSO")
        lab = f"scRNA_DiD_580YxDHA_{t}"
        if lab not in have and all(c in cpm.columns for c in
                [f"MRA1251_{t}DHA", f"MRA1251_{t}DMSO", f"MRA1250_{t}DHA", f"MRA1250_{t}DMSO"]):
            y = ((cpm[f"MRA1251_{t}DHA"] - cpm[f"MRA1251_{t}DMSO"]) -
                 (cpm[f"MRA1250_{t}DHA"] - cpm[f"MRA1250_{t}DMSO"])).dropna()
            y = y[y.index.isin(FEAT.index)]
            if len(y) >= 100 and y.nunique() >= 10:
                print(f"  {lab} n={len(y)} med={y.median():.3f}", flush=True)
                ALL += ftests(y, lab)
                ALL += ftests(no_schiz(y), lab + "_noSchiz")
    save_state(ALL); checkpoint("D2_scRNA", ALL)
    print(f"MINE D2 DONE T={round(time.time()-t0,1)}", flush=True)

if STEP == "E1E2":
    print("E1: Mok...", flush=True)
    L1h = pd.read_csv("data/derived/WP4/L1h_mok_genotype_context/L1h_long.tsv", sep="\t")
    for ct, j in L1h.groupby("contrast"):
        y = j.set_index("gene")["log2FC"].dropna()
        y = y[y.index.isin(FEAT.index)]
        lab = f"Mok2021_{ct}"
        if lab not in have and len(y) >= 100 and y.nunique() >= 10:
            ALL += ftests(y, lab)
    save_state(ALL); checkpoint("E1_Mok", ALL)
    print("E2: GSE59099...", flush=True)
    try:
        pmap = pd.read_csv("data/derived/WP3/M3-03_dha/platform_GPL18893_probe2gene.tsv", sep="\t")
        import subprocess
        tbl = subprocess.run(["zcat", "data/raw/geo/GSE59099/GSE59099_series_matrix.txt.gz"],
                             capture_output=True, text=True, timeout=300).stdout.splitlines()
        beg = next(i for i, l in enumerate(tbl) if l.startswith("!series_matrix_table_begin"))
        hdr = tbl[beg + 1].replace('"', "").split("\t")
        dat = [l.replace('"', "").split("\t") for l in tbl[beg + 2:]
               if l and not l.startswith("!series_matrix_table_end")]
        G = pd.DataFrame(dat, columns=hdr).set_index("ID_REF")
        meta_rows = [l.replace('"', "").split("\t") for l in tbl if l.startswith("!Sample_characteristics_ch1")]
        chars = {}
        for row in meta_rows:
            for v in row[1:]:
                if ":" in v:
                    k, vv = v.split(":", 1); chars.setdefault(k.strip(), []).append(vv.strip())
        nsamp = len(hdr) - 1
        chars = {k: v for k, v in chars.items() if len(v) == nsamp}
        print(f"  usable chars keys={[k for k in chars]}", flush=True)
        if not any("half" in k.lower() for k in chars):
            hl_pos = [None] * nsamp
            for l in tbl:
                if l.startswith("!Sample_characteristics_ch1") and "half" in l.lower():
                    row = l.replace('"', "").split("\t")
                    cells = (row[1:] + [""] * nsamp)[:nsamp]
                    for i, v in enumerate(cells):
                        if ":" in v and hl_pos[i] is None:
                            hl_pos[i] = v.split(":", 1)[1].strip().rstrip("h").strip()
            chars["halflife_merged"] = hl_pos
            print(f"  halflife_merged filled={sum(1 for v in hl_pos if v)}", flush=True)
        CH = pd.DataFrame(chars)
        CH.index = [c.strip('"') for c in hdr[1:]]
        hl_key = next((k for k in CH.columns if "half" in k.lower()), None)
        hl = pd.to_numeric(CH[hl_key].str.extract(r"([\d.]+)")[0] if hl_key else pd.Series(dtype=float), errors="coerce") if hl_key else pd.Series(np.nan, index=CH.index)
        tp = CH.get("timepoint"); geo = CH.get("geographic origin")
        pre = (tp == "prior to artemisinin combination therapy (ACT)")
        mek = geo.isin(["Pailin, Cambodia", "Mae Sot, Thailand", "Preah Vihear, Cambodia",
                        "Binh Phuoc, Vietnam", "Rattanakiri, Cambodia", "Pursat, Cambodia",
                        "Attapeu, Laos", "Shwe Kyin, Myanmar", "Sisakhet, Thailand",
                        "Ranong, Thailand", "Khun Han, Thailand"])
        R = (hl >= 5) & pre; S = (hl <= 3) & pre
        print(f"  pre-ACT n={int(pre.sum())} R={int(R.sum())} S={int(S.sum())}", flush=True)
        pc = pmap.columns
        probe_col = pc[0]
        gene_col = [c for c in pc if "gene" in c.lower() or "PF3D7" in c.lower()][0]
        pg = dict(zip(pmap[probe_col].astype(str), pmap[gene_col].astype(str)))
        Gnum = G.apply(pd.to_numeric, errors="coerce")
        def gene_level(cols):
            sub = Gnum[cols]
            by = {}
            for pr in sub.index:
                g = pg.get(str(pr), "")
                if isinstance(g, str) and g.startswith("PF3D7"):
                    by.setdefault(g.split(".")[0], []).append(pr)
            rows = {g: sub.loc[prs].median(axis=0).values for g, prs in by.items() if prs}
            return pd.DataFrame(rows, index=cols)
        def g99_test(colsR, colsS, lab):
            global ALL
            if lab in have:
                return
            GR = gene_level(colsR); GS = gene_level(colsS)
            common = GR.columns.intersection(GS.columns).intersection(FEAT.index)
            y = pd.Series(GR[common].median(axis=0) - GS[common].median(axis=0), index=common)
            if len(y) < 100 or y.nunique() < 10:
                print(f"  SKIP {lab} degenerate", flush=True); return
            print(f"  {lab} genes={len(y)} med={y.median():.4f}", flush=True)
            ALL += ftests(y, lab)
            ALL += ftests(no_schiz(y), lab + "_noSchiz")
        Rcols = [c for c in G.columns if R.loc[c]]; Scols = [c for c in G.columns if S.loc[c]]
        g99_test(Rcols, Scols, "GSE59099_RvsS_preACT")
        g99_test([c for c in Rcols if mek.loc[c]], [c for c in Scols if mek.loc[c]], "GSE59099_RvsS_Mekong")
    except Exception as e:
        print(f"  E2 FAILED: {type(e).__name__}: {e}", flush=True)
        with open(f"{OUT}/M4RC_E2_NOTE.txt", "w") as h:
            h.write(f"GSE59099 blocked: {type(e).__name__}: {e}\nStatus: UNRESOLVED-gap.\n")
    save_state(ALL); checkpoint("E1E2", ALL)
    print(f"MINE E1E2 DONE T={round(time.time()-t0,1)}", flush=True)

if STEP == "F":
    print("F: dTE interaction...", flush=True)
    try:
        import tarfile, glob
        os.makedirs("/tmp/m4rc/dte", exist_ok=True)
        if not glob.glob("/tmp/m4rc/dte/*.tabular.txt.gz"):
            tar = tarfile.open("data/raw/geo/GSE226632/GSE226632_RAW.tar")
            tar.extractall("/tmp/m4rc/dte")
        fmap = {}
        for f in glob.glob("/tmp/m4rc/dte/*.tabular.txt.gz"):
            b = os.path.basename(f)
            m = re.search(r"htseq-count-(CM|AA_free)-(total|polysome)-replicate(\d)", b)
            if m:
                fmap[(m.group(1), m.group(2), int(m.group(3)))] = f
        print(f"  files={len(fmap)}", flush=True)
        CNT = {}
        for k, f in fmap.items():
            d = pd.read_csv(f, sep="\t", header=None, index_col=0, names=["count"])
            d.index = d.index.str.replace(r"\.\d+$", "", regex=True)
            CNT[k] = d["count"]
        C = pd.DataFrame(CNT).fillna(0)
        C = C.loc[:, sorted(C.columns)]
        gm = np.exp(np.log(C.where(C > 0).replace(0, np.nan)).mean(axis=1, skipna=True))
        sf = (C.div(gm, axis=0)).median(axis=0)
        Y = np.log2(C.div(sf, axis=1) + 1)
        cols = list(Y.columns)
        X = pd.DataFrame({"icept": 1}, index=cols)
        X["cond"] = [1 if c[0] == "AA_free" else 0 for c in cols]
        X["frac"] = [1 if c[1] == "polysome" else 0 for c in cols]
        X["CxF"] = X.cond * X.frac
        for r in [2, 3]:
            X[f"rep{r}"] = [1 if c[2] == r else 0 for c in cols]
        Xm = X.values; XtXinv = np.linalg.inv(Xm.T @ Xm)
        df_res = len(cols) - Xm.shape[1]
        inter, pse = {}, {}
        for g in Y.index:
            yv = Y.loc[g].values
            if yv.max() < 1:
                continue
            beta, *_ = np.linalg.lstsq(Xm, yv, rcond=None)
            resid = yv - Xm @ beta
            s2 = (resid @ resid) / df_res
            se = math.sqrt(s2 * XtXinv[3, 3])
            inter[g] = beta[3]; pse[g] = beta[3] / se if se > 0 else 0.0
        I = pd.Series(inter); T = pd.Series(pse)
        p_int = pd.Series({g: 2 * stats.t.sf(abs(t), df_res) for g, t in T.items()})
        q_int = pd.Series(bh_fdr(p_int.values), index=p_int.index)
        DTE = pd.DataFrame({"dTE_interaction": I, "t": T, "p": p_int, "q": q_int})
        RATIO = pd.read_csv("data/derived/WP4/L1_codon_x_stress/L1_GSE226632_dTE.tsv",
                            sep="\t").set_index("gene")["dTE"]
        common = DTE.index.intersection(RATIO.index)
        agree = np.sign(DTE.loc[common, "dTE_interaction"]) == np.sign(RATIO.loc[common])
        sp = stats.spearmanr(DTE.loc[common, "dTE_interaction"], RATIO.loc[common], nan_policy="omit")
        DTE.to_csv(f"{OUT}/M4RC_dTE_interaction.tsv", sep="\t")
        print(f"  genes={len(DTE)} common={len(common)} sign_agree={float(agree.mean()):.3f} "
              f"spearman={float(sp.statistic):.3f} (p={sp.pvalue:.2g})", flush=True)
        print(f"  q<0.05 n={int((q_int < 0.05).sum())}", flush=True)
        if "GSE226632_dTEinteraction" not in have and DTE.dTE_interaction.nunique() >= 10:
            ALL += ftests(DTE.dTE_interaction, "GSE226632_dTEinteraction")
        dte_verdict = dict(sign_agreement=float(agree.mean()), spearman=float(sp.statistic),
                           spearman_p=float(sp.pvalue), n_q05=int((q_int < 0.05).sum()))
        with open(f"{OUT}/M4RC_dTE_upgrade.json", "w") as h:
            json.dump(dte_verdict, h, indent=2)
    except Exception as e:
        print(f"  F FAILED: {type(e).__name__}: {e}", flush=True)
        with open(f"{OUT}/M4RC_dTE_NOTE.txt", "w") as h:
            h.write(f"dTE interaction blocked: {type(e).__name__}: {e}\n"
                    "Status: existing ratio-based dTE retained as DISCOVERY/SUPPORTING per packet C4.\n")
    save_state(ALL); checkpoint("F_dTE", ALL)
    print(f"MINE F DONE T={round(time.time()-t0,1)}", flush=True)

if STEP == "FINAL":
    RES = pd.DataFrame(ALL)
    RES.to_pickle(f"{OUT}/M4RC_feature_tests_all.pkl")
    sub = RES[RES.feature == "asn_frac"].copy()
    def arm(c):
        if c.startswith("GSE151189") or c.startswith("scRNA") or c.startswith("GSE225340"):
            return "acute"
        if c.startswith("Mok2021") or c.startswith("GSE59099"):
            return "chronic"
        return "other"
    sub["arm"] = sub.contrast.map(arm)
    sub.to_csv(f"{OUT}/M4RC_acute_chronic_partition.tsv", sep="\t", index=False)
    RES.to_csv(f"{OUT}/M4RC_feature_decomposition.tsv", sep="\t", index=False)
    print(f"FINAL contrasts={RES.contrast.nunique()} rows={len(RES)}", flush=True)
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(sub[["contrast", "arm", "n", "spearman_rho", "spearman_q"]].to_string(), flush=True)
    print("MINE FINAL DONE", flush=True)
