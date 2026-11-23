# SSR_populations — population-genetics panels and tables

R pipeline for the SSR-based population genetic analyses and the panels of
the manuscript: multilocus genotype (MLG) assignment, genotypic diversity,
clonality (index of association), minimum spanning networks (MSN),
ordination (PCA, PCoA), UPGMA trees, AMOVA, DAPC (unsupervised clustering +
supervised cross-validation), the geographic distribution maps, and the
composition of the final Figure 1 panel.

Corresponds to the "Population genetic analysis" subsection of Materials
and Methods. This README documents the inputs, parameters and usage needed
to reproduce those results; the manuscript text itself only summarizes them.

## Contents

```
SSR_populations/
├── haplotype_network_v4.R   all analyses and panel figures (run in R / RStudio)
├── make_supp_table_structure.py   Supplementary Table S2 (../Sup_material/TableS2_EU2A1_population_structure.xlsx)
├── make_sampling_year_panels.R    isolates per sampling year, one panel per file (outputs/s1/)
├── Amplicon_pop/            NJ phylogeny (Newick) from the 98-amplicon SNPs (own README)
├── inputs/
│   └── full_dataset_SSR.csv  analysis-ready isolate table (read by the script)
└── outputs/                  written by the scripts (not in the repository; see ../build_all.sh)
    ├── s1/   dataset 1 (all isolates): Figure 1, S1 analyses, maps
    └── s2/   dataset 2 (EU_2_A1 only): sub-population structure
```

## Requirements

Developed and run with R 4.3.1 (Windows). Packages (CRAN unless noted), with
the versions used:

| Used by | Packages |
|---|---|
| Analyses (`library()` at the top of the script) | poppr 2.9.8, adegenet 2.1.11, ade4 1.7.24, ape 5.8.1, ggplot2 4.0.3, dplyr 1.1.3, tidyr 1.3.2, plotly 4.12.0, knitr 1.51, maps 3.4.3, scatterpie 0.2.6, ggforce 0.5.0, scales 1.4.0 |
| UPGMA tree (`upgma()`) | phangorn 2.12.1 |
| Argentina isolate map (last block of the script) | sf 1.1.0, rnaturalearth 1.2.0, ggspatial 1.1.11, ggrepel 0.9.6, cowplot 1.2.0, viridisLite 0.4.3 |
| Saving | svglite 2.2.2 (`.svg` output of `ggsave()`), gtable 0.3.6 (fixed panel size) |

```r
install.packages(c("poppr", "adegenet", "ade4", "ape", "phangorn", "ggplot2",
                   "dplyr", "tidyr", "plotly", "knitr", "maps", "scatterpie",
                   "ggforce", "scales", "sf", "rnaturalearth", "ggspatial",
                   "ggrepel", "cowplot", "viridisLite", "svglite", "gtable", "png"))
```

`ggtree` is **not** needed: the Figure 1 tree is drawn directly with ggplot2
(the `ggtree` 3.8.x release is not compatible with ggplot2 4.x).

## Input data

`inputs/full_dataset_SSR.csv` — one row per isolate (2,187 isolates), read as
is: the script does no further filtering. Columns:

- `Sample` — isolate ID (duplicates are disambiguated with `make.unique()`
  at load time).
- `Genotype` — MLG/lineage label (e.g. `EU_2_A1`, `EC_1`, `EU_13_A2`, `PE_3`,
  `US_1`, `US_7_A2`, `EU_8_A1`, `BR_1`, `PE_7.1`, ...).
- `Country` — sampling country (Euroblight reference-panel isolates carry
  placeholder codes that the script pools as `EU_controls`).
- `Origin`, `Year`, `Latitude`, `Longitude` — origin label, sampling year
  (integer; used for temporal grouping; 394 of the 2,187 isolates have no
  year, and they are kept in every figure) and
  coordinates.
- `source` — where the isolate comes from (`this work`, `Kreuze 2017`,
  `Patarroyo 2020`, `Izarra 2025`, `Izarra 2026`, `EU_controls`). The script
  collapses it into two provenance classes: **This work** and **Context
  samples** (everything else); see `provenance_category()`, which also accepts
  the legacy label `Flor_Lucca_check` of earlier tables.
- 12 SSR loci, named exactly `D13`, `PinfSSR8`, `PinfSSR4`, `Pi04`, `Pi70`,
  `PinfSSR6`, `Pi63`, `PiG11`, `Pi02`, `PinfSSR11`, `PinfSSR2`, `Pi4B`.
  Genotype calls are triploid, alleles separated by `/`
  (`df2genind(..., sep = "/", ploidy = 3)`).

Bruvo's distance needs the repeat length of each locus (`replen`, Section 1):
ten of the twelve are dinucleotide, while `Pi63` and `Pi70` are trinucleotide
and get 3 bp. `RL()` subsets that vector **by locus name**, because
`missingno()` and `poppr.amova()` drop loci with more than 5% missing data and a
vector indexed by position would then give the wrong repeat length to the loci
that remain.

Genotype labels are folded into one set of categories (`Pop_geno_map`,
`map_genotype_category()`): `EU_2_A1` (incl. `PE_7.1`), `EC1` (incl. `EC_1`),
`EU_13_A2`, `US7` (incl. `US_7`, `US_7_A2`), `US_1`, `PE_3`,
`BR_1 / EU_8_A1`, and `Other` for everything else. Every figure that colors
by genotype uses this categorization and the palette `map_geno_colors`;
countries use `country_colors`; provenance uses `provenance_colors`
(blue / red).

## Usage

1. Set the R working directory to this folder (all reads and writes use
   relative paths).
2. Run `haplotype_network_v4.R` from the top, interactively in R or RStudio
   (several sections `print()` plots for inline inspection, so it is not
   designed to be run unattended end to end). Sections 1-3 (parameters, data
   loading, helper functions) must run first. Sections 7, 8 and 9 are
   self-contained and can be re-run on their own afterwards (each rebuilds
   its own genind and Bruvo-distance objects).

Runtime: the index-of-association test (`ia_test()`) and the AMOVA
(`amova_analysis()`) use 999 permutations and are the slowest steps; the
`run_analysis(s1, ...)` call alone takes several minutes. Random seeds are
set wherever a stochastic step occurs (MSN layout, `find.clusters()`,
`xvalDapc()`, AMOVA permutations), so results are deterministic for a given
input file.

## What each part of the script does

- **Section 1 — Global parameters**: locus list, color and linetype palettes
  (genotype, country, provenance), the genotype categorization and the
  provenance rule, country coordinates for the maps.
- **Section 2 — Data loading**: reads `inputs/full_dataset_SSR.csv` and
  builds two working subsets with `prep_subset()`: `s1` (the whole table) and
  `s2` (`EU_2_A1` + `PE_7.1`). Each carries `Pop_geno` (lineage, `PE_7.1`
  folded into `EU_2_A1`), `Pop_geno_map` (unified categories), `Pop_country`
  and `Period`.
- **Section 3 — Helper functions**: MSN, PCA, PCoA
  (`pcoa_publication_plot()`), UPGMA, AMOVA, index of association, DAPC and
  cluster validation, diversity export, the map functions, bar plots, and
  `run_analysis()`, which chains them for a given subset.
- **Section 4 — Analysis 1** (`run_analysis(s1, ...)`): MSN, PCA, PCoA
  (without ellipses), UPGMA, index of association and diversity on `s1`;
  DAPC and AMOVA are switched off for this dataset.
- **Section 5 — Analysis 2** (`run_analysis(s2, ...)`): the full battery on
  `s2`, plus the supervised Argentina+Uruguay vs. rest DAPC with
  cross-validation.
- **Section 6 — Maps**: two maps over every isolate of the table (no filter
  on sampling year), by genotype (Map A) and by provenance (`This work` vs.
  `Context samples`, Map B). The `post2005` in the file names is historical:
  there is no year filter. Map B labels each country with its oldest-newest
  sampling year (isolates without a year are ignored for that label).
- **Section 7 — Figure 1 panels** (details below): the two PCoA panels, the
  UPGMA tree with color columns, and the flat AMOVA on `Pop_geno_map`.
- **Section 8 — EU_2_A1 sub-structure**: PCoA on `s2` with ellipses over
  three a priori populations (Argentina+Uruguay / Panama / rest), the
  matching AMOVA, and a country-colored MSN.
- **Section 9 — PCoA with the map genotypes**: `add_map_pcoa()` on `s1`
  (EU_2_A1, EC1, EU_13_A2, US7, US_1, PE_3 and the rest as `Other`), colored
  by country and by genotype, without ellipses.
- **Section 10 — Locus-level differentiation**: per-locus allele-frequency
  comparison (Panama vs. rest of EU_2_A1; Argentina+Uruguay vs. rest) and
  DAPC allele loadings for the Argentina+Uruguay contrast.
- **Section 11 — Data-driven structure of EU_2_A1 (S2)**: standalone block that lets
  the structure come out of the genotypes without choosing groups. (1) Unsupervised
  clustering with `adegenet::find.clusters()` (k-means on the principal components of the
  allele frequencies, K compared with the BIC; it is *not* a Bayesian method and never sees
  the country), run on all isolates and on clone-corrected data (one isolate per MLG within
  each country): K chosen by `diffNgroup` over 30 random seeds, BIC curve, reproducibility
  across random starts (adjusted Rand index) and the clusters x country tables for K = 2, 3
  and 4. (2) Pairwise PhiST between countries (AMOVA on clone-corrected genotypes, Bruvo
  distance, the same 11 loci for every pair, 9999 permutations, Holm-adjusted p) and a
  hierarchical clustering of the PhiST matrix. The clusters are described by their
  composition; they should not be used as groups for an AMOVA or cross-validation on the
  same data (k-means maximizes between-cluster separation, so those tests are significant
  by construction). It also repeats the cross-validation of the a priori hypothesis
  Argentina+Uruguay vs the rest (`xvalDapc`) on all isolates and on clone-corrected data.
  `python make_supp_table_structure.py` then writes Supplementary Table S2
  (`../Sup_material/TableS2_EU2A1_population_structure.xlsx`): the Section 11 sheets plus the
  poppr diversity, index-of-association and AMOVA output of dataset S2 (`outputs/s2/*.csv`).
- **Last block — Argentina isolate map**: location of the seven isolates
  used in the infection assays, with their response on cv. Innovator
  (PiNSL-19, Ar21009, Biocca and Cba2: hypersensitive response, avirulent;
  Ar20079, Ar20112 and Tafi: virulent).

## Panels of Figure 1

The panels the manuscript's Figure 1 is made of, each written on its own:

- **A** — `mapa_EU2A1_EC1_post2005.png`: map of all isolates by genotype
  (Section 6, Map A), drawn with the shared genotype palette and no legend.
- **B** — `Fig1b_by_genotype.png` (left) and `Fig1b_by_country.png` (right):
  the same PCoA of `s1` (Bruvo distance, `cmdscale`) colored by genotype and by
  country, placed side by side under a single panel letter. One circle per
  MLG with radius scaled by the number of isolates (`pie_radius_range`), alpha
  0.75, MLGs spanning more than one group drawn as pie charts. No legends and
  no ellipses: colors are explained by the legends panel and the dominant
  lineages already separate along PCo1/PCo2. Both plots share axis limits,
  panel size (4.4 in wide) and canvas size (5.4 × 5 in), so they line up
  exactly and point for point.
- **C** — `Fig1_UPGMA_table.png`: UPGMA tree (Bruvo distance) of all 2,187
  isolates of `s1`, no tip labels and no legends, using the full height of
  the page. Each leaf carries three color columns — **Country**,
  **Genotype** and **Source** (blue = This work, red = Context samples).
  Saved at 4 × 10 in, the exact size of its box, so it is not rescaled.
  Built by `upgma_color_table_plot()`.
- **Legends** — `Fig1_legends.png` (no panel letter): the Country, Genotype
  and Source color keys in a single column, in the free space to the right of
  the map, between A and C. It is the single color reference for the whole
  figure (map, PCoAs and tree share the palettes `map_geno_colors`,
  `country_colors` and `provenance_colors`). Built by
  `upgma_color_legend_plot()`; `upgma_color_table_plot(legends_below = TRUE)`
  can still draw them under the tree instead.

Each image is placed at its native aspect ratio (contain-fit, centered)
inside its box, so some empty margin is expected (e.g. the portrait map in
its box). The tree enters Figure 1 as a raster; the vector version is
`Fig1_UPGMA_table.pdf` / `.svg`.

## Outputs

All figures are saved as `.png`, `.pdf` and `.svg`. The definitive manuscript figures and
Table S2 is written outside this folder, to `../Sup_material/`:

| Manuscript item | Script | File |
|---|---|---|
| Supplementary Table S1 | `inputs/full_dataset_SSR.csv` (copied) | `../Sup_material/TableS1_isolates_and_SSR_genotypes.csv` |
| Supplementary Table S2 | `make_supp_table_structure.py` | `../Sup_material/TableS2_EU2A1_population_structure.xlsx` |

The map `outputs/s2/mapa_argentina_innovator.png` is panel B of Figure 4 (built by
`../omics_and_virulence/omics_and_virulence.ipynb`).

**`outputs/s1/`** (dataset 1)

- `Fig1_UPGMA_table` — Figure 1C tree with the color columns.
- `Fig1_legends` — Figure 1 color legends (Country / Genotype / Source).
- `Fig1b_by_genotype`, `Fig1b_by_country` — the two PCoA plots of Figure 1B.
- `mapa_EU2A1_EC1_post2005`, `mapa_provenance_post2005` — Maps A (genotype,
  no legend, Figure 1A) and B (provenance, with legend).
- `Fig_mapgenotypes_df1_PCoA_mapgenotypes_{country,genotype}` — Section 9.
- `sampling_years_*` (`.pdf`, `.png`) — isolates per sampling year, one file per panel:
  all isolates by genotype and by country, EU_2_A1 by genotype and by country.
- `S1_PCoA_by_{country,geno}`, `UPGMA_S1_by_{country,geno}` (`.pdf`/`.svg`),
  `UPGMA_totpop.svg` — Section 4 exploratory figures.
- `amova_S1_*.csv`, `ia_S1.csv`, `diversity_S1_per_{locus,genotype,country,year}.csv`.

**`outputs/s2/`** (dataset 2, EU_2_A1)

- `Fig_EU2A1_PCoA_3subpops` — Section 8.
- `S2_PCoA_by_{country,geno}`, `UPGMA_S2_by_{country,geno}` (`.pdf`/`.svg`),
  `UPGMA_2A1.svg`.
- `mapa_argentina_innovator` — Argentina isolate map.
- `unsupervised/` (Section 11): `k_by_seed.csv`, `k_chosen_summary.csv`, `bic_by_k.csv` and
  `bic_by_k` (figure), `cluster_stability_ari.csv`, `clusters_x_country_long.csv` and one
  `clusters_x_country_<dataset>_K<k>.csv` per dataset and K, `pairwise_phiST_long.csv`,
  `pairwise_phiST_matrix.csv`, `phiST_heatmap` and `phiST_dendrogram` (figures),
  `xval_hypothesis_ARUY.csv` (cross-validation, all isolates vs clone-corrected); all of them
  are collected in Supplementary Table S2 by `make_supp_table_structure.py`.
- `amova_S2_EU2A1_3subpop_Pop_3pop.csv`, `ia_S2.csv`,
  `diversity_S2_per_{locus,genotype,country,year}.csv`.

Plots that are only printed to the active graphics device (MSN, PCA, DAPC,
bar plots, Section 10 loadings) are not saved to file; add `ggsave()` calls
where a standalone file is needed.
