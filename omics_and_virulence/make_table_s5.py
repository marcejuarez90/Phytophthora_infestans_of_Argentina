#!/usr/bin/env python3
# ══════════════════════════════════════════════════════════════════════════════
#  make_table_s5.py — Supplementary Table S5, gene expression of PiNSL-19
#
#  One workbook with the two quantities the manuscript reports for every
#  annotated T30-4 gene:
#
#    DE   log2 fold change and p value of each infection timepoint relative to
#         in vitro mycelium (inputs/tabla_LFC_pvalor.tsv). The source names its
#         timepoints T1-T5; only T3, T4 and T5 are the 24, 48 and 72 hpi
#         contrasts the manuscript uses, so T1 and T2 are dropped here.
#    TPM  the StringTie TPM matrix the expression panels are drawn from
#         (inputs/TPM_matrix.tsv), one column per library, plus the per-
#         condition means that Figure 3C plots as log2(TPM + 1).
#
#  Library names encode the design: the letter is the experiment (A, B, C) and
#  the digit the timepoint (3 = 24 hpi, 4 = 48 hpi, 5 = 72 hpi); Ctrl* are the
#  mycelium controls. Nine libraries: two of mycelium, two at 24 hpi, three at
#  48 hpi and two at 72 hpi.
#
#  Run from this folder:   python make_table_s5.py
#  Writes ../Sup_material/TableS5_gene_expression_PiNSL19.xlsx
# ══════════════════════════════════════════════════════════════════════════════
import os

import pandas as pd

DE_TSV = "inputs/tabla_LFC_pvalor.tsv"
TPM_TSV = "inputs/TPM_matrix.tsv"
OUT = "../Sup_material/TableS5_gene_expression_PiNSL19.xlsx"

# the three contrasts the manuscript reports, in the source's own naming
TIMEPOINTS = {"T3": "24hpi", "T4": "48hpi", "T5": "72hpi"}

# library -> condition, from the design encoded in the column names
CONDITION = {"Ctrl1": "Myc", "Ctrl2": "Myc",
             "B3": "24hpi", "C3": "24hpi",
             "A4": "48hpi", "B4": "48hpi", "C4": "48hpi",
             "B5": "72hpi", "C5": "72hpi"}
COND_ORDER = ["Myc", "24hpi", "48hpi", "72hpi"]


def differential_expression():
    de = pd.read_csv(DE_TSV, sep="\t")
    de.index.name = None
    de = de.reset_index().rename(columns={"index": "gene_id"})
    if de.columns[0] != "gene_id":                 # the ids are the first column
        de = de.rename(columns={de.columns[0]: "gene_id"})
    keep = ["gene_id"]
    rename = {}
    for src, label in TIMEPOINTS.items():
        for suffix, new in (("LFC", "log2FC"), ("pval", "pvalue")):
            col = f"{src}{suffix}"
            assert col in de.columns, f"{col} is not a column of {DE_TSV}"
            keep.append(col)
            rename[col] = f"{label}_{new}"
    dropped = [c for c in de.columns if c not in keep]
    print(f"DE: {len(de)} genes; dropped {dropped}")
    return de[keep].rename(columns=rename)


def tpm_matrix():
    tpm = pd.read_csv(TPM_TSV, sep="\t")
    libs = [c for c in tpm.columns if c != "gene_id"]
    unknown = sorted(set(libs) - set(CONDITION))
    assert not unknown, f"libraries with no condition: {unknown}"

    ordered = sorted(libs, key=lambda c: (COND_ORDER.index(CONDITION[c]), c))
    out = tpm[["gene_id"] + ordered].rename(
        columns={c: f"{CONDITION[c]}_{c}" for c in ordered})
    for cond in COND_ORDER:                        # what Figure 3C averages
        cols = [f"{CONDITION[c]}_{c}" for c in ordered if CONDITION[c] == cond]
        out[f"mean_{cond}"] = out[cols].mean(axis=1)
    n = {c: sum(1 for x in libs if CONDITION[x] == c) for c in COND_ORDER}
    print(f"TPM: {len(out)} genes; libraries per condition {n}")
    return out


def main():
    de, tpm = differential_expression(), tpm_matrix()
    assert set(de["gene_id"]) == set(tpm["gene_id"]), \
        "the two sheets do not cover the same genes"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with pd.ExcelWriter(OUT, engine="openpyxl") as xls:
        de.to_excel(xls, sheet_name="DE", index=False)
        tpm.to_excel(xls, sheet_name="TPM", index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
