from __future__ import annotations

import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, pearsonr

# ---------------------------------------------------------------------------
# paths  
ROOT = Path(__file__).resolve().parent
BASE = ROOT / "GSE190504"

EXPR_XLSX = BASE / "GSE190504_Processed_Data_Spreadsheet_Glioma_Study.xlsx"
META_TXT = BASE / "GSE190504_series_matrix.txt"
REF_CSV = ROOT / "bulk_ref" / "reference_celltype_mean.csv"

OUT_DIR = BASE / "glio_deconv_bulkref_cellanneal"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PROP = OUT_DIR / "GSE190504_glio_cellanneal_proportions.csv"
OUT_STATS = OUT_DIR / "GSE190504_glio_cellanneal_stats.txt"

CELLANNEAL_REPO = ROOT / "cellanneal_repo"
sys.path.insert(0, str(CELLANNEAL_REPO))
 
from scipy.optimize import dual_annealing as _scipy_dual_annealing  
from cellanneal import general as _ca_general   
_ca_general.dual_annealing = _scipy_dual_annealing

from cellanneal.pipelines import run_cellanneal  


# cellanneal hyper-parameter
DISP_MIN = 0.5
BULK_MIN = 1e-8
BULK_MAX = 1.0
MAXITER = 500

# ---------------------------------------------------------------------------
# self-contained XLSX reader (avoids the openpyxl dependency)
# ---------------------------------------------------------------------------
_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _col_letters_to_index(letters: str) -> int:
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch.upper()) - ord("A") + 1)
    return idx - 1


def read_xlsx_sheet(xlsx_path: Path, sheet_index: int = 0) -> list[list[str | None]]:
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
        tree = ET.parse(zf.open(sheet_names[sheet_index]))
        rows = tree.getroot().find("s:sheetData", _NS)

        out: list[list[str | None]] = []
        max_cols = 0
        for row in rows.findall("s:row", _NS):
            cells: dict[int, str | None] = {}
            for c in row.findall("s:c", _NS):
                col_letters = re.match(r"[A-Z]+", c.attrib["r"]).group(0)
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


def sanitize_ct(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name))


def clean_chr(values) -> list[str]:
    out = []
    for v in values:
        s = "" if v is None else str(v)
        s = s.strip()
        s = re.sub(r'^"|"$', "", s)
        out.append(s.strip())
    return out


def parse_series_metadata(path: Path) -> pd.DataFrame:
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
    if len(chars) < 5:
        raise RuntimeError("Expected >=5 !Sample_characteristics_ch1 lines.")

    histology = [re.sub(r"^histology:\s*", "", v) for v in chars[2]]
    genotype = [re.sub(r"^genotype:\s*", "", v) for v in chars[3]]
    treatment = [re.sub(r"^treatment:\s*", "", v) for v in chars[4]]

    meta = pd.DataFrame({
        "sample_id": titles,
        "histology": histology,
        "genotype": genotype,
        "treatment": treatment,
    })
    meta["idh1"] = np.where(meta["genotype"].str.contains("IDH1-mut", regex=False), "IDH1-mut", "IDH1-wt")
    meta["treatment_simple"] = np.where(meta["treatment"] == "Untreated", "Untreated", "Treated")
    return meta


def load_expression(xlsx_path: Path) -> pd.DataFrame:
    """Return gene x sample DataFrame, gene symbols as index, linear scale."""
    sheet = read_xlsx_sheet(xlsx_path, sheet_index=0)
    if not sheet:
        raise RuntimeError("XLSX appears empty")

    ensg = [s if s is not None else "" for s in sheet[0][2:]]
    sym = [(s.upper() if s is not None else "") for s in sheet[1][2:]]
    for k in range(len(sym)):
        if not sym[k]:
            sym[k] = ensg[k].upper()

    data_rows = sheet[5:]
    sample_title = [r[1] if len(r) > 1 else None for r in data_rows]
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
    return expr


def load_reference(ref_csv: Path) -> pd.DataFrame:
    ref = pd.read_csv(ref_csv, index_col=0)
    ref.index = ref.index.astype(str).str.upper().str.replace(r"\..*$", "", regex=True)
    ref = ref.loc[(ref.index != "") & (~ref.index.isna())]
    ref = ref[~ref.index.duplicated(keep="first")]
    ref.columns = [sanitize_ct(c) for c in ref.columns]
    ref = ref.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return ref


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"[info] reading expression: {EXPR_XLSX}")
    expr = load_expression(EXPR_XLSX)
    print(f"[info] expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    print(f"[info] reading metadata: {META_TXT}")
    meta = parse_series_metadata(META_TXT)

    glio_samples = [s for s in expr.columns if s in set(meta.loc[meta["histology"] == "Glio", "sample_id"])]
    bulk = expr.loc[:, glio_samples].copy()
    print(f"[info] Glio samples: {bulk.shape[1]}")

    print(f"[info] reading reference: {REF_CSV}")
    ref = load_reference(REF_CSV)
    print(f"[info] reference: {ref.shape[0]} genes x {ref.shape[1]} cell types")

    common = bulk.index.intersection(ref.index)
    print(f"[info] common genes (raw): {len(common)}")
    bulk2 = bulk.loc[common].copy()
    ref2 = ref.loc[common].copy()

    keep = (bulk2.sum(axis=1) > 0) & (ref2.sum(axis=1) > 0)
    bulk2 = bulk2.loc[keep]
    ref2 = ref2.loc[keep]
    print(f"[info] common genes after non-zero filter: {bulk2.shape[0]}")
    if bulk2.shape[0] < 1000:
        raise RuntimeError("Too few common genes after filtering; aborting.")

    # composition normalisation per sample / per cell type, matching the
    # GSE222515 cellanneal_bulkref pipeline.
    bulk2 = bulk2.div(bulk2.sum(axis=0), axis=1).fillna(0.0)
    ref2 = ref2.div(ref2.sum(axis=0), axis=1).fillna(0.0)

    print(f"[info] running cellanneal (maxiter={MAXITER}) ...")
    mix = run_cellanneal(
        celltype_df=ref2,
        bulk_df=bulk2,
        disp_min=DISP_MIN,
        bulk_min=BULK_MIN,
        bulk_max=BULK_MAX,
        maxiter=MAXITER,
    )

    res = mix.copy().reset_index().rename(columns={"index": "sample_id"})
    # standardise CAF column name
    caf_col = next((c for c in res.columns if c.lower() == "cafs"), None)
    if caf_col is None:
        raise RuntimeError("CAFs column not found in cellanneal output")
    if caf_col != "CAFs":
        res = res.rename(columns={caf_col: "CAFs"})

    res = res.merge(meta[["sample_id", "idh1", "treatment_simple", "histology"]], on="sample_id", how="left")
    res.to_csv(OUT_PROP, index=False)
    print(f"[saved] {OUT_PROP}")

    # ---- statistics
    caf_wt = res.loc[res["idh1"] == "IDH1-wt", "CAFs"].to_numpy()
    caf_mut = res.loc[res["idh1"] == "IDH1-mut", "CAFs"].to_numpy()
    p_idh = mannwhitneyu(caf_wt, caf_mut, alternative="two-sided").pvalue if len(caf_wt) and len(caf_mut) else np.nan

    caf_unt = res.loc[res["treatment_simple"] == "Untreated", "CAFs"].to_numpy()
    caf_trt = res.loc[res["treatment_simple"] == "Treated", "CAFs"].to_numpy()
    p_trt = mannwhitneyu(caf_unt, caf_trt, alternative="two-sided").pvalue if len(caf_unt) and len(caf_trt) else np.nan

    # report cellanneal fit QC if present
    qc_cols = [c for c in res.columns if c.lower() in {"rho_spearman", "rho_pearson"}]
    qc_summary = ""
    if qc_cols:
        rho_s = res["rho_Spearman"].median() if "rho_Spearman" in res.columns else np.nan
        rho_p = res["rho_Pearson"].median() if "rho_Pearson" in res.columns else np.nan
        qc_summary = f"Median rho_Spearman: {rho_s:.4f} | Median rho_Pearson: {rho_p:.4f}"

    lines = [
        f"Run: {pd.Timestamp.now()}",
        f"Glio samples: {len(res)}",
        f"Common genes used: {len(common)} (after nonzero filter: {bulk2.shape[0]})",
        f"Reference cell types: {', '.join(ref2.columns)}",
        f"Median CAF proportion: {res['CAFs'].median():.6f}",
        f"Nonzero CAF count: {int((res['CAFs'] > 0).sum())}",
        f"IDH1 wt vs mut p (MWU): {fmt_p(p_idh)} (n_wt={len(caf_wt)}, n_mut={len(caf_mut)})",
        f"Untreated vs Treated p (MWU): {fmt_p(p_trt)} (n_unt={len(caf_unt)}, n_trt={len(caf_trt)})",
    ]
    if qc_summary:
        lines.append(qc_summary)

    OUT_STATS.write_text("\n".join(lines) + "\n")
    print(f"[saved] {OUT_STATS}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
