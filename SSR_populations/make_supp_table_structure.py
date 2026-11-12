#!/usr/bin/env python3
"""
make_supp_table_structure.py
Builds Supplementary Table S2 (structure and diversity of the EU_2_A1 lineage, dataset S2)
from the CSV files written by haplotype_network_v4.R: Section 11 (outputs/s2/unsupervised/)
and the poppr diversity / index of association / AMOVA output of dataset S2 (outputs/s2/).

    python make_supp_table_structure.py

Output: ../Sup_material/TableS2_EU2A1_population_structure.xlsx
  Sheets: Clusters x country | K selection | Pairwise PhiST | Cross-validation | Notes |
          Diversity by country | Diversity by locus | Diversity by genotype | Diversity by year |
          Index of association | AMOVA three groups
Requires pandas and openpyxl.
"""
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

D = os.path.join("outputs", "s2", "unsupervised")
D2 = os.path.join("outputs", "s2")
OUT = os.path.join("..", "Sup_material", "TableS2_EU2A1_population_structure.xlsx")
COUNTRIES = ["Argentina", "Uruguay", "Panama", "Bolivia", "Brazil", "Chile", "Peru", "UK"]
DATASETS = {"all_isolates": "All isolates",
            "clone_corrected": "Clone-corrected (one isolate per MLG within each country)"}
rd = lambda f: pd.read_csv(os.path.join(D, f))
B = Font(bold=True)


def put(ws, r, values, bold=False):
    for c, v in enumerate(values, 1):
        cell = ws.cell(row=r, column=c, value=v)
        if bold:
            cell.font = B
    return r + 1


def put_df(ws, r, df, title=None):
    if title:
        ws.cell(row=r, column=1, value=title).font = Font(bold=True, size=12)
        r += 1
    r = put(ws, r, list(df.columns), bold=True)
    for row in df.itertuples(index=False):
        r = put(ws, r, [None if (isinstance(v, float) and pd.isna(v)) else v for v in row])
    return r + 1


def fit(ws, widths=None):
    for i, col in enumerate(ws.columns, 1):
        w = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[get_column_letter(i)].width = min(max(10, w + 2), 60)
    ws.freeze_panes = "A2"


wb = Workbook()

# 1. clusters x country ------------------------------------------------------
ws = wb.active
ws.title = "Clusters x country"
long = rd("clusters_x_country_long.csv")
rows = []
for ds, name in DATASETS.items():
    sub = long[long["dataset"] == ds]
    n_ds = int(sub[sub["K"] == "K2"]["Freq"].sum())
    for k in ["K2", "K3", "K4"]:
        t = sub[sub["K"] == k].pivot(index="Cluster", columns="Country", values="Freq").fillna(0).astype(int)
        t = t.reindex(columns=[c for c in COUNTRIES if c in t.columns])
        for cl, r in t.iterrows():
            rows.append([f"{name} (n = {n_ds})", int(k[1:]), int(cl)] + [int(r.get(c, 0)) for c in COUNTRIES]
                        + [int(r.sum())])
put(ws, 1, ["Dataset", "K", "Cluster"] + COUNTRIES + ["Total"], bold=True)
for i, r in enumerate(rows, 2):
    put(ws, i, r)
fit(ws)

# 2. K selection -------------------------------------------------------------
ws = wb.create_sheet("K selection")
ks = rd("k_chosen_summary.csv").copy()
ks["dataset"] = ks["dataset"].map(DATASETS)
ks.columns = ["Dataset", "K chosen by diffNgroup", "Number of seeds (of 30)"]
r = put_df(ws, 1, ks, "K chosen by find.clusters (diffNgroup) over 30 random seeds")
bic = rd("bic_by_k.csv")
bic["dataset"] = bic["dataset"].map(DATASETS)
bw = bic.pivot(index="dataset", columns="K", values="BIC").round(1).reset_index()
bw.columns = ["Dataset"] + [f"K = {c}" for c in bw.columns[1:]]
r = put_df(ws, r, bw, "BIC by K (seed 42)")
st = rd("cluster_stability_ari.csv").copy()
st["dataset"] = st["dataset"].map(DATASETS)
st = st.rename(columns={"dataset": "Dataset", "n_starts": "Random starts", "ARI_min": "ARI (minimum)",
                        "ARI_median": "ARI (median)"})
put_df(ws, r, st, "Reproducibility of the partition across random starts (adjusted Rand index; 1 = identical)")
ws.column_dimensions["A"].width = 62
for col in "BCDEFGHI":
    ws.column_dimensions[col].width = 20

# 3. pairwise PhiST ----------------------------------------------------------
ws = wb.create_sheet("Pairwise PhiST")
ph = rd("pairwise_phiST_long.csv").copy()
ph["PhiST"] = ph["PhiST"].round(3)
ph["p"] = ph["p"].round(4)
ph["p_holm"] = ph["p_holm"].round(4)
ph.columns = ["Country 1", "Country 2", "n genotypes 1", "n genotypes 2", "PhiST", "p (9999 permutations)",
              "p (Holm-adjusted)"]
r = put_df(ws, 1, ph.sort_values("PhiST"), "Pairwise PhiST between countries (clone-corrected; sorted by PhiST)")
mat = pd.read_csv(os.path.join(D, "pairwise_phiST_matrix.csv"), index_col=0).round(3)
mat.insert(0, "", mat.index)
put_df(ws, r, mat, "PhiST matrix (negative values set to 0)")
fit(ws)

# 4. cross-validation --------------------------------------------------------
ws = wb.create_sheet("Cross-validation")
x = rd("xval_hypothesis_ARUY.csv").copy()
x.columns = ["Dataset", "Hypothesis (a priori, by country)", "Isolates", "Isolates in hypothesis group",
             "Isolates in other group", "MLGs", "MLGs shared between groups", "Best number of PCs",
             "Mean correct assignment (%)", "Replicates"]
put_df(ws, 1, x)
fit(ws)

# 5. notes -------------------------------------------------------------------
ws = wb.create_sheet("Notes")
notes = [
    "Supplementary tables for the structure of the EU_2_A1 lineage (dataset S2).",
    "Clustering: adegenet::find.clusters (k-means on the principal components of allele frequencies, up to 50 PCs, 10 starts, 100,000 iterations);",
    "  K chosen with the diffNgroup criterion; the BIC is only used to compare K. It is not a Bayesian method and does not use country labels.",
    "Two versions of the data: all isolates (641; two groups with fewer than 3 isolates excluded) and clone-corrected (146; one isolate per MLG within each country).",
    "Cluster numbers are arbitrary labels. Clusters are described by their composition; they are not used as groups for tests on the same data.",
    "Pairwise PhiST: AMOVA on clone-corrected genotypes (poppr.amova, ade4), Bruvo distance, the same 11 loci for every pair (Pi63 removed for >5% missing),",
    "  9999 permutations, Holm correction for 28 comparisons. UK = European reference isolates. Countries with fewer than 5 genotypes are not included.",
    "Cross-validation (xvalDapc): hypothesis Argentina + Uruguay vs the rest, defined a priori by country; 100 replicates, 90/10 training/test, up to 50 PCs.",
    "  The clone-corrected row keeps one isolate per MLG within each group.",
    "Source files: outputs/s2/unsupervised/*.csv (written by Section 11 of haplotype_network_v4.R).",
    None,
    None,
    "Sheets added from the poppr output of the same analysis:",
    "  Diversity by country / locus / genotype / year - genotypic diversity and clonality indices.",
    "  Index of association - Ia and rbarD for the whole lineage.",
    "  AMOVA three groups - variance partition of the a priori country-based three-group scheme.",
]
for i, t in enumerate(notes, 1):
    ws.cell(row=i, column=1, value=t)
ws.column_dimensions["A"].width = 150

# 6-11. poppr output of dataset S2 (outputs/s2/*.csv) -------------------------
def poppr_sheet(name, csv, title, first_col, col_a_width, footnotes=()):
    ws = wb.create_sheet(name)
    df = pd.read_csv(os.path.join(D2, csv)).round(4)
    df.columns = [first_col] + list(df.columns[1:])
    df[first_col] = df[first_col].astype(str)
    ws.cell(row=1, column=1, value=title).font = Font(bold=True, size=12)
    r = put_df(ws, 2, df)
    for t in footnotes:
        ws.cell(row=r, column=1, value=t)
        r += 1
    ws.column_dimensions["A"].width = col_a_width
    for i in range(2, df.shape[1] + 1):
        ws.column_dimensions[get_column_letter(i)].width = 10
    ws.freeze_panes = "A3"


poppr_sheet("Diversity by country", "diversity_S2_per_country.csv",
            "Genotypic diversity and clonality of EU_2_A1 by country (poppr)", "Pop", 70, [
                "N, isolates; MLG, multilocus genotypes; eMLG, rarefied MLG richness (SE); H, Shannon index;",
                "G, Stoddart and Taylor's G; lambda, Simpson's 1-lambda; E.5, evenness; Hexp, expected heterozygosity;",
                "Ia and rbarD, index of association and its standardized form (999 permutations).",
                "Rows with N = 1 (ref, United Kingdom(Wales)) are retained for completeness but are not interpretable."])
poppr_sheet("Diversity by locus", "diversity_S2_per_locus.csv",
            "Allelic diversity of EU_2_A1 by SSR locus", "Group", 43)
poppr_sheet("Diversity by genotype", "diversity_S2_per_genotype.csv",
            "Genotypic diversity of EU_2_A1 (whole lineage)", "Group", 48)
poppr_sheet("Diversity by year", "diversity_S2_per_year.csv",
            "Genotypic diversity of EU_2_A1 by sampling year", "Group", 49)
poppr_sheet("Index of association", "ia_S2.csv",
            "Index of association for the whole EU_2_A1 lineage (999 permutations)", None, 70, [
                "Ia, index of association; rbarD, standardized index of association; p.Ia and p.rD, permutation p-values."])
# The percentage of variance, Std.Obs and p come from the randtest printed by the AMOVA
# step of haplotype_network_v4.R (they are not in the CSV).
poppr_sheet("AMOVA three groups", "amova_S2_EU2A1_3subpop_Pop_3pop.csv",
            "AMOVA of the a priori three-group scheme within EU_2_A1 (Argentina+Uruguay, Panama, rest)",
            "Source of variation", 70, [
                "Clone-corrected data, Bruvo's distance, 999 Monte Carlo permutations (poppr.amova, ade4 method).",
                "35.75% of the total molecular variance was explained between the three groups (Std.Obs = 64.1, p = 0.001)."])

wb.save(OUT)
print("wrote", OUT)
