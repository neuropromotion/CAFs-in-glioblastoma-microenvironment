#!/usr/bin/env python3
"""
GSE190504 — rank-based hCAF score (Python port of gse190504_hCAF_within_full.R).

Reproduces the full pipeline from raw inputs:
    * GSE190504_Processed_Data_Spreadsheet_Glioma_Study.xlsx (expression matrix)
    * GSE190504_series_matrix.txt                            (sample metadata)
    * hCAF.txt                                               (CAF gene signature)

Produces (in within_full_hCAF/):
    * GSE190504_hCAF_scores_metadata.csv      per-sample rank/z scores + metadata
    * GSE190504_hCAF_longitudinal_pairs.csv   TP1/TP2 paired table
    * GSE190504_hCAF_test_table.csv           KW / Wilcoxon p-values (with BH)
    * GSE190504_hCAF_within_stats.txt         human-readable summary

For each sample the score is the mean within-sample percentile rank of hCAF
genes (rank/N over all measured genes). This is robust to scale/normalisation
differences and is what the article tables use. We also report a complementary
z-score: mean of per-gene standardised log2(expr+1) across the hCAF set.

Run:
    python3 gse190504_hCAF_rank_score.py
"""

from __future__ import annotations

import io
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, wilcoxon


def bh_adjust(pvals) -> np.ndarray:
    """Benjamini-Hochberg FDR adjusted p-values (no statsmodels dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
BASE = Path("/Users/neuropromotion/Desktop/CAF/BULK/final/GSE190504")
ROOT = BASE.parent

EXPR_XLSX = BASE / "GSE190504_Processed_Data_Spreadsheet_Glioma_Study.xlsx"
META_TXT = BASE / "GSE190504_series_matrix.txt"
HCAF_FILE = ROOT / "hCAF.txt"

OUT_DIR = BASE / "within_full_hCAF"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_SCORES = OUT_DIR / "GSE190504_hCAF_scores_metadata.csv"
OUT_PAIRS = OUT_DIR / "GSE190504_hCAF_longitudinal_pairs.csv"
OUT_PVALS = OUT_DIR / "GSE190504_hCAF_test_table.csv"
OUT_STATS = OUT_DIR / "GSE190504_hCAF_within_stats.txt"


# ---------------------------------------------------------------------------
# self-contained XLSX reader (no openpyxl required)
# ---------------------------------------------------------------------------
_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _col_letters_to_index(letters: str) -> int:
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch.upper()) - ord("A") + 1)
    return idx - 1


def read_xlsx_sheet(xlsx_path: Path, sheet_index: int = 0) -> list[list[str | None]]:
    """Return a dense 2-D list of cell strings for the requested sheet."""
    with zipfile.ZipFile(xlsx_path) as zf:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            tree = ET.parse(zf.open("xl/sharedStrings.xml"))
            for si in tree.getroot().findall("s:si", _NS):
                pieces = [t.text or "" for t in si.iter("{%s}t" % _NS["s"])]
                shared.append("".join(pieces))

        sheet_names = sorted(
            n for n in zf.namelist() if n.startswith("xl/worksheets/sheet")
        )
        sheet_file = sheet_names[sheet_index]
        tree = ET.parse(zf.open(sheet_file))
        rows = tree.getroot().find("s:sheetData", _NS)

        out: list[list[str | None]] = []
        max_cols = 0
        for row in rows.findall("s:row", _NS):
            cells: dict[int, str | None] = {}
            for c in row.findall("s:c", _NS):
                ref = c.attrib["r"]
                col_letters = re.match(r"[A-Z]+", ref).group(0)
                col_idx = _col_letters_to_index(col_letters)
                t = c.attrib.get("t", "n")
                v_node = c.find("s:v", _NS)
                if t == "s":
                    cells[col_idx] = shared[int(v_node.text)] if v_node is not None else None
                elif t == "inlineStr":
                    is_node = c.find("s:is", _NS)
                    txt = "".join((tn.text or "") for tn in (is_node.iter("{%s}t" % _NS["s"]) if is_node is not None else []))
                    cells[col_idx] = txt
                else:
                    cells[col_idx] = v_node.text if v_node is not None else None
            if cells:
                max_cols = max(max_cols, max(cells) + 1)
            out.append([cells.get(i) for i in range(max_cols)])

        # pad ragged rows so the matrix is rectangular
        for r in out:
            if len(r) < max_cols:
                r.extend([None] * (max_cols - len(r)))
        return out


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p) or np.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.4e}"
    return f"{p:.4f}"


def clean_chr(values) -> list[str]:
    out = []
    for v in values:
        s = "" if v is None else str(v)
        s = s.strip()
        s = re.sub(r'^"|"$', "", s)
        out.append(s.strip())
    return out


def parse_series_metadata(path: Path) -> pd.DataFrame:
    """Pull sample-level rows out of a GEO series_matrix.txt."""
    lines = path.read_text(errors="replace").splitlines()

    def get_row(prefix: str) -> list[str]:
        for ln in lines:
            if ln.startswith(prefix):
                return clean_chr(ln.split("\t")[1:])
        return []

    def get_rows(prefix: str) -> list[list[str]]:
        return [clean_chr(ln.split("\t")[1:]) for ln in lines if ln.startswith(prefix)]

    titles = get_row("!Sample_title")
    chars = get_rows("!Sample_characteristics_ch1")
    if len(chars) < 6:
        raise RuntimeError("Expected >=6 !Sample_characteristics_ch1 lines.")

    diagnosis = [re.sub(r"^diagnosis:\s*", "", v) for v in chars[0]]
    tissue = [re.sub(r"^tissue:\s*", "", v) for v in chars[1]]
    histology = [re.sub(r"^histology:\s*", "", v) for v in chars[2]]
    genotype = [re.sub(r"^genotype:\s*", "", v) for v in chars[3]]
    treatment = [re.sub(r"^treatment:\s*", "", v) for v in chars[4]]
    paired_info = chars[5]

    meta = pd.DataFrame({
        "sample_title": titles,
        "diagnosis": diagnosis,
        "tissue": tissue,
        "histology": histology,
        "genotype": genotype,
        "treatment": treatment,
        "paired_info": paired_info,
    })
    meta["idh1"] = np.where(meta["genotype"].str.contains("IDH1-mut", regex=False), "IDH1-mut", "IDH1-wt")
    meta["treatment_simple"] = np.where(meta["treatment"] == "Untreated", "Untreated", "Treated")
    meta["histology"] = pd.Categorical(meta["histology"], categories=["Oligo", "Astro", "Glio"], ordered=True)

    pair_re = re.compile(r"paired patien[t]?:\s*(\S+)")
    tp_re = re.compile(r"time point\s*([12])")
    meta["pair_ref"] = meta["paired_info"].apply(lambda s: (pair_re.search(s).group(1) if pair_re.search(s) else ""))
    meta["timepoint_ref"] = meta["paired_info"].apply(
        lambda s: int(tp_re.search(s).group(1)) if tp_re.search(s) else np.nan
    )

    # propagate own timepoint: if sample X says "paired with Y, time point K",
    # then Y is at TP=K and X is at the other timepoint.
    self_tp = pd.Series(np.nan, index=meta.index, dtype=float)
    title_to_idx = {t: i for i, t in enumerate(meta["sample_title"])}
    for i, row in meta.iterrows():
        ref = row["pair_ref"]
        tpr = row["timepoint_ref"]
        if ref and not np.isnan(tpr) and ref in title_to_idx:
            j = title_to_idx[ref]
            self_tp.iloc[j] = tpr
            self_tp.iloc[i] = 1.0 if tpr == 2 else 2.0
    meta["timepoint_self"] = self_tp
    return meta


def load_expression(xlsx_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Parse the processed expression spreadsheet.

    Layout (1-based as in the R script):
        row 1, col 3..N  -> ENSG IDs
        row 2, col 3..N  -> gene symbols
        row 6..M, col 1  -> sample short id
        row 6..M, col 2  -> sample title (matches series_matrix.txt)
        row 6..M, col 3.. -> expression values (CPM-like, linear scale)

    Returns gene x sample DataFrame (deduped/aggregated) and sample short ids.
    """
    sheet = read_xlsx_sheet(xlsx_path, sheet_index=0)
    if not sheet:
        raise RuntimeError("XLSX appears empty")

    ensg = [s if s is not None else "" for s in sheet[0][2:]]
    sym = [(s.upper() if s is not None else "") for s in sheet[1][2:]]
    for k in range(len(sym)):
        if not sym[k]:
            sym[k] = ensg[k].upper()

    data_rows = sheet[5:]
    sample_short = [r[0] if len(r) > 0 else None for r in data_rows]
    sample_title = [r[1] if len(r) > 1 else None for r in data_rows]
    sample_short = [s if s is not None else "" for s in sample_short]
    sample_title = [s if s is not None else "" for s in sample_title]

    values = np.full((len(data_rows), len(sym)), np.nan, dtype=float)
    for i, r in enumerate(data_rows):
        for j in range(len(sym)):
            v = r[2 + j] if 2 + j < len(r) else None
            if v is None or v == "":
                continue
            try:
                values[i, j] = float(v)
            except ValueError:
                continue
    values = np.nan_to_num(values, nan=0.0)

    expr = pd.DataFrame(values.T, index=pd.Index(sym, name="symbol"), columns=sample_title)
    expr.index = expr.index.astype(str).str.upper().str.replace(r"\..*$", "", regex=True)
    expr = expr.loc[(expr.index != "") & (~expr.index.isna())]
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0, sort=False).sum()
    return expr, sample_short


# ---------------------------------------------------------------------------
# scoring + statistics
# ---------------------------------------------------------------------------
def rank_and_z_scores(expr: pd.DataFrame, sig_genes: list[str]) -> tuple[pd.Series, pd.Series, list[str]]:
    """Return (rank_score, z_score, present_signature_genes)."""
    present = [g for g in sig_genes if g in expr.index]
    if len(present) < 5:
        raise RuntimeError(f"Too few signature genes overlap with expression: {len(present)}")

    n_genes = expr.shape[0]
    ranks = expr.rank(axis=0, method="average") / n_genes
    rank_score = ranks.loc[present].mean(axis=0)

    log_expr = np.log2(expr.clip(lower=0) + 1.0)
    sub = log_expr.loc[present]
    means = sub.mean(axis=1)
    sds = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(means, axis=0).div(sds, axis=0).fillna(0.0)
    z_score = z.mean(axis=0)

    return rank_score, z_score, present


def add_test_row(rows: list[dict], analysis: str, metric: str, comparison: str,
                 n1: int, n2: int | None, test: str, p_value: float, p_bh: float | None = None) -> None:
    rows.append({
        "analysis": analysis,
        "metric": metric,
        "group_or_comparison": comparison,
        "n1": int(n1) if n1 is not None else None,
        "n2": int(n2) if n2 is not None else None,
        "test": test,
        "p_value": p_value,
        "p_bh": p_bh,
    })


def multi_group_tests(rows: list[dict], df: pd.DataFrame, metric: str, grp: str, analysis_name: str) -> None:
    dd = df[[metric, grp]].dropna()
    dd = dd[dd[grp].astype(str) != ""]
    levels = [lv for lv in dd[grp].astype(str).unique() if lv]
    if len(levels) < 2:
        return

    groups_vals = [dd.loc[dd[grp].astype(str) == lv, metric].to_numpy() for lv in levels]
    kw_p = kruskal(*groups_vals).pvalue
    add_test_row(rows, analysis_name, metric, "all groups", len(dd), None, "Kruskal-Wallis", kw_p)

    if len(levels) > 2:
        pairs = list(combinations(levels, 2))
        pvs: list[float] = []
        for a, b in pairs:
            va = dd.loc[dd[grp].astype(str) == a, metric].to_numpy()
            vb = dd.loc[dd[grp].astype(str) == b, metric].to_numpy()
            pvs.append(mannwhitneyu(va, vb, alternative="two-sided").pvalue)
        adj = bh_adjust(pvs)
        for (a, b), p, q in zip(pairs, pvs, adj):
            n1 = int((dd[grp].astype(str) == a).sum())
            n2 = int((dd[grp].astype(str) == b).sum())
            add_test_row(rows, analysis_name, metric, f"{a} vs {b}", n1, n2, "Wilcoxon", p, q)
    else:
        a, b = levels[0], levels[1]
        va = dd.loc[dd[grp].astype(str) == a, metric].to_numpy()
        vb = dd.loc[dd[grp].astype(str) == b, metric].to_numpy()
        p = mannwhitneyu(va, vb, alternative="two-sided").pvalue
        add_test_row(rows, analysis_name, metric, f"{a} vs {b}", len(va), len(vb), "Wilcoxon", p)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"[info] reading expression: {EXPR_XLSX}")
    expr, sample_short = load_expression(EXPR_XLSX)
    print(f"[info] expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    print(f"[info] reading metadata: {META_TXT}")
    meta = parse_series_metadata(META_TXT)

    # align expression columns to metadata sample titles where possible
    keep_cols = [c for c in expr.columns if c in set(meta["sample_title"])]
    expr = expr.loc[:, keep_cols]

    print(f"[info] reading hCAF signature: {HCAF_FILE}")
    hcaf = [ln.strip().upper() for ln in HCAF_FILE.read_text().splitlines() if ln.strip()]
    hcaf = sorted(set(hcaf))

    rank_score, z_score, present = rank_and_z_scores(expr, hcaf)
    print(f"[info] hCAF overlap: {len(present)} / {len(hcaf)}")

    title_to_short = dict(zip(expr.columns, [sample_short[i] for i, t in enumerate(expr.columns) if i < len(sample_short)]))
    short_lookup = {t: s for t, s in zip(meta["sample_title"], sample_short[: len(meta)])}
    scores = pd.DataFrame({
        "sample_title": rank_score.index,
        "sample_short_id": [short_lookup.get(t, "") for t in rank_score.index],
        "CAF_rank_score": rank_score.values,
        "CAF_z_score": z_score.reindex(rank_score.index).values,
    })
    dat = scores.merge(meta, on="sample_title", how="left")
    dat.to_csv(OUT_SCORES, index=False)
    print(f"[saved] {OUT_SCORES}")

    # ---- longitudinal pair table
    edges = meta.loc[meta["pair_ref"].astype(bool), ["sample_title", "pair_ref"]].copy()
    edges["key"] = edges.apply(
        lambda r: "|".join(sorted([r["sample_title"], r["pair_ref"]])), axis=1
    )
    edges = edges.drop_duplicates("key")

    tp_lookup = dict(zip(dat["sample_title"], dat["timepoint_self"]))
    pair_rows = []
    seen_keys: set[str] = set()
    for _, r in edges.iterrows():
        a, b = r["sample_title"], r["pair_ref"]
        ta = tp_lookup.get(a, np.nan)
        tb = tp_lookup.get(b, np.nan)
        if pd.isna(ta) or pd.isna(tb) or ta == tb:
            continue
        tp1 = a if ta == 1 else b
        tp2 = a if ta == 2 else b
        k = "|".join(sorted([tp1, tp2]))
        if k in seen_keys:
            continue
        seen_keys.add(k)
        pair_rows.append({"tp1": tp1, "tp2": tp2})
    pairs = pd.DataFrame(pair_rows)

    if not pairs.empty:
        get = lambda s, col: dat.set_index("sample_title")[col].reindex(s).to_numpy()
        pairs["CAF_rank_tp1"] = get(pairs["tp1"], "CAF_rank_score")
        pairs["CAF_rank_tp2"] = get(pairs["tp2"], "CAF_rank_score")
        pairs["CAF_z_tp1"] = get(pairs["tp1"], "CAF_z_score")
        pairs["CAF_z_tp2"] = get(pairs["tp2"], "CAF_z_score")
        pairs["delta_rank_tp2_minus_tp1"] = pairs["CAF_rank_tp2"] - pairs["CAF_rank_tp1"]
        pairs["delta_z_tp2_minus_tp1"] = pairs["CAF_z_tp2"] - pairs["CAF_z_tp1"]
        pairs["histology_tp1"] = get(pairs["tp1"], "histology").astype(str)
        pairs["histology_tp2"] = get(pairs["tp2"], "histology").astype(str)
        pairs["idh1_tp1"] = get(pairs["tp1"], "idh1")
        pairs["idh1_tp2"] = get(pairs["tp2"], "idh1")
    pairs.to_csv(OUT_PAIRS, index=False)
    print(f"[saved] {OUT_PAIRS} ({len(pairs)} pairs)")

    # ---- tests
    rows: list[dict] = []

    multi_group_tests(rows, dat, "CAF_rank_score", "histology", "within_main")
    multi_group_tests(rows, dat, "CAF_rank_score", "idh1", "within_main")
    multi_group_tests(rows, dat, "CAF_rank_score", "treatment_simple", "within_main")
    multi_group_tests(rows, dat, "CAF_z_score", "histology", "within_main")
    multi_group_tests(rows, dat, "CAF_z_score", "idh1", "within_main")
    multi_group_tests(rows, dat, "CAF_z_score", "treatment_simple", "within_main")

    for idh in ("IDH1-mut", "IDH1-wt"):
        dsub = dat[dat["idh1"] == idh]
        if len(dsub) > 8 and dsub["histology"].dropna().astype(str).nunique() >= 2:
            multi_group_tests(rows, dsub, "CAF_rank_score", "histology", f"histology_within_{idh}")

    for h in ("Oligo", "Astro", "Glio"):
        dsub = dat[dat["histology"].astype(str) == h]
        if len(dsub) > 8 and dsub["idh1"].dropna().nunique() == 2:
            multi_group_tests(rows, dsub, "CAF_rank_score", "idh1", f"idh_within_{h}")

    if len(pairs) >= 3:
        wr = wilcoxon(pairs["CAF_rank_tp2"], pairs["CAF_rank_tp1"], alternative="two-sided")
        wz = wilcoxon(pairs["CAF_z_tp2"], pairs["CAF_z_tp1"], alternative="two-sided")
        wr1 = wilcoxon(pairs["CAF_rank_tp2"], pairs["CAF_rank_tp1"], alternative="greater")
        add_test_row(rows, "longitudinal_paired", "CAF_rank_score", "TP2 vs TP1", len(pairs), len(pairs), "paired Wilcoxon", wr.pvalue)
        add_test_row(rows, "longitudinal_paired", "CAF_z_score", "TP2 vs TP1", len(pairs), len(pairs), "paired Wilcoxon", wz.pvalue)
        add_test_row(rows, "longitudinal_paired", "CAF_rank_score", "TP2 > TP1 (one-sided)", len(pairs), len(pairs), "paired Wilcoxon", wr1.pvalue)

    tests = pd.DataFrame(rows)
    if not tests.empty:
        tests["p_bh_final"] = tests["p_bh"]
        key = tests["analysis"] + "||" + tests["metric"]
        for k, idx in tests.groupby(key).groups.items():
            mask = tests.loc[idx, "p_bh_final"].isna()
            need = tests.loc[idx][mask]
            if len(need) > 1:
                adj = bh_adjust(need["p_value"].to_numpy())
                tests.loc[need.index, "p_bh_final"] = adj
            elif len(need) == 1:
                tests.loc[need.index, "p_bh_final"] = need["p_value"].to_numpy()
    tests.to_csv(OUT_PVALS, index=False)
    print(f"[saved] {OUT_PVALS}")

    # ---- human-readable summary
    lines = [
        f"Run: {pd.Timestamp.now()}",
        f"Samples: {len(dat)}",
        f"hCAF overlap: {len(present)} / {len(hcaf)}",
        f"Inferred longitudinal pairs: {len(pairs)}",
        "",
    ]
    for _, r in tests.iterrows():
        lines.append(
            f"[{r['analysis']}] {r['metric']} | {r['group_or_comparison']} | "
            f"{r['test']} | p={fmt_p(r['p_value'])}, BH={fmt_p(r['p_bh_final'])}"
        )
    OUT_STATS.write_text("\n".join(lines) + "\n")
    print(f"[saved] {OUT_STATS}")


if __name__ == "__main__":
    main()
