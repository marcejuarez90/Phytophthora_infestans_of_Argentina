# RT-qPCR validation of the PiNSL-19 RNA-seq

Compares the log2 fold change of four genes measured by RT-qPCR with the value
obtained from the RNA-seq of the same experiment, at 24, 48 and 72 hpi.

## Files

| file | content |
|---|---|
| `qPCR_vs_RNAseq.csv`       | tidy input: one row per gene x timepoint |
| `qPCR_vs_RNAseq.R`         | statistics and figure, reproducible from the CSV alone |
| `qPCR_vs_RNAseq.png`      | preview of the regression panel (the PDF is `qPCR_vs_RNAseq.pdf`) |
| `qPCR_vs_RNAseq_stats.txt` | printed output of the script |

## Provenance of the two columns

**`qpcr_log2FC`** — from the RT-qPCR workbook `qPCR_deltadelta.xlsx` (raw Ct values; not included here), as
`log2FC = -DDCt`. DDCt is the value in the `deltadelta` cell of the left-hand
block of each gene sheet: three technical replicates per infection timepoint
normalised to the geometric mean of the reference gene PITG_10941 (EF2a) and
calibrated against the mean DCt of the in vitro mycelium samples (`avM`).
Note that the column labelled `lfc` in the workbook is in fact the fold change
2^-DDCt, not a log2 fold change; it is not used here.
PITG_09173 is taken from the `final` sheet, which already stores log2FC.

**`rnaseq_log2FC`** — columns `T3LFC`, `T4LFC`, `T5LFC` of
`../inputs/tabla_LFC_pvalor.tsv`, the differential-expression table that the DE
sheet of Supplementary Table S5 is built from (`../make_table_s5.py`)
(T3 = 24 hpi, T4 = 48 hpi, T5 = 72 hpi), i.e. each infection timepoint against
the in vitro mycelium control, the same contrast as the qPCR. `rnaseq_p` is the
matching `T*pval` column.

## Genes

The workbook holds sheets for PITG_08943, PITG_13567, PITG_04339 and PITG_07059
plus a `final` sheet for PITG_09173. PITG_04339 is excluded here. The sheet
named `13567` carries the title "PITG_16567" in its first cell, which is not a
gene of the T30-4 annotation; the sheet name PITG_13567 is used.

## Run

    Rscript qPCR_vs_RNAseq.R

Needs only ggplot2 and cowplot.
