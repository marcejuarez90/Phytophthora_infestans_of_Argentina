# Phytophthora infestans in Argentina — analysis code

Code and input data behind the figures and supplementary tables of the manuscript on the
EU_2_A1 clonal lineage of *Phytophthora infestans* in Southern South America.

```
SSR_populations/      SSR population genetics (R): Figures 1-2, Supplementary Figures S1-S2, Tables S1-S2
omics_and_virulence/  genome, RNA-seq and infection assays (Jupyter + R): Figures 3-5, Supplementary Figures S3-S5, Table S5
build_all.sh          runs everything in order
```

Each folder has its own README with requirements, inputs and usage.

```
bash build_all.sh     # see the header of the script to point it at specific R / Python / Jupyter installs
```

Every panel is built from the data in this repository. **The page layout of the manuscript —
how panels are composed into a figure, their sizing and typography — is not part of this
repository**: what is here reproduces each panel and each result on its own.

| Item | Built by |
|---|---|
| Figure 1 | `SSR_populations/make_figura1.R`, from the panels of `haplotype_network_v4.R` |
| Figure 2 | panels from `SSR_populations/haplotype_network_v4.R` (the published figure is assembled by hand) |
| Figures 3 and 4 | `omics_and_virulence/omics_and_virulence.ipynb` |
| Figure 5 panels | `omics_and_virulence.ipynb` (lesion assay) and `omics_and_virulence/avrblb2_family/` (AVR-blb2 family tree and coverage; SNP-205 genotyping) |
| Supplementary Figure S1 | `SSR_populations/make_figsup_sampling_years.R` |
| Supplementary Figure S2 | `SSR_populations/make_figsup_amplicon_phist.R` |
| Supplementary Figure S3 | `omics_and_virulence/omics_and_virulence.ipynb` |
| Supplementary Figure S4 | `omics_and_virulence/qPCR_validation/qPCR_vs_RNAseq.R` |
| Supplementary Figure S5 | `omics_and_virulence/avr3a_allele/plot_avr3a.py` |
| Supplementary Table S1 (isolates and SSR genotypes) | copy of `SSR_populations/inputs/full_dataset_SSR.csv` |
| Supplementary Table S2 (EU_2_A1 structure and diversity) | `SSR_populations/make_supp_table_structure.py` |
| Supplementary Table S5 (differential expression of PiNSL-19) | copy of `omics_and_virulence/inputs/tabla_LFC_pvalor.tsv` |

Supplementary Tables S3, S4 and S6 (sequencing and alignment metrics, primer lists, effector
annotation) are not produced by any script and are not in this repository; they are part of the
manuscript's supplementary material. The effector annotation they summarize is here as
`omics_and_virulence/inputs/SD_effector_list.xlsx`.

## Isolate table

`SSR_populations/inputs/full_dataset_SSR.csv` holds the 2,187 isolates the manuscript analyses.
It does not include 55 isolates of this work from Ecuador and Costa Rica that are not yet
cleared for publication.
