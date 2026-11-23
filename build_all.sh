#!/usr/bin/env bash
# Rebuilds every panel and every code-generated supplementary table from the raw
# inputs. Panels land in each folder's outputs/, the tables in Sup_material/.
# Run from any folder:   bash build_all.sh
#
# No figure is assembled here: how panels are arranged into the figures of the
# manuscript is a presentation choice and is not part of this repository.
#
# Interpreters can be overridden, e.g. on Windows/WSL:
#   RSCRIPT="/mnt/c/Program Files/R/R-4.3.1/bin/Rscript.exe" \
#   PYTHON=/mnt/c/Users/<user>/anaconda3/python.exe \
#   JUPYTER=/mnt/c/Users/<user>/anaconda3/Scripts/jupyter-nbconvert.exe bash build_all.sh
#
# SSR_populations/haplotype_network_v4.R runs every SSR analysis and writes the panels and CSVs
# (SSR_populations/outputs/) that every SSR panel comes from; it is the slow step.
set -euo pipefail
RSCRIPT=${RSCRIPT:-Rscript}
PYTHON=${PYTHON:-python}
JUPYTER=${JUPYTER:-jupyter-nbconvert}
cd "$(dirname "$0")"
mkdir -p Sup_material

# ── SSR population genetics ──────────────────────────────────────────────────
cd SSR_populations
"$RSCRIPT" haplotype_network_v4.R          # outputs/s1, outputs/s2 (panels, CSVs, Argentina map)
"$RSCRIPT" Amplicon_pop/plot_amplicon_tree.R Amplicon_pop/results_excl_08_7054A/ar_eu/ar_eu_nj.nwk \
    Amplicon_pop/sample_genotypes.tsv Amplicon_pop/results_excl_08_7054A/ar_eu/ar_eu_nj_tree_panel \
    midpoint 70 6.5 8.4 bottomleft         # Supplementary Figure S2, panel A
"$RSCRIPT" make_sampling_year_panels.R     # outputs/s1/sampling_years_*.pdf
"$PYTHON" make_supp_table_structure.py     # Sup_material/TableS2_EU2A1_population_structure.xlsx
cp inputs/full_dataset_SSR.csv ../Sup_material/TableS1_isolates_and_SSR_genotypes.csv
cd ..

# ── Genomics, transcriptomics and virulence assays ───────────────────────────
cd omics_and_virulence
"$JUPYTER" --to notebook --execute --inplace --ExecutePreprocessor.timeout=1200 \
    omics_and_virulence.ipynb              # every panel of the genome, RNA-seq and assay sections -> outputs/
"$PYTHON" make_table_s5.py                 # Sup_material/TableS5_gene_expression_PiNSL19.xlsx (DE + TPM sheets)
cd qPCR_validation
"$RSCRIPT" qPCR_vs_RNAseq.R                # RT-qPCR vs RNA-seq regression
cd ../avrblb2_family
"$PYTHON" avrblb2_family.py                # outputs/avrblb2_family.csv and .pdf (AVR-blb2 family tree + coverage)
"$PYTHON" plot_snp205.py                   # outputs/fig5B_snp205.pdf (SNP-205 genotyping, Figure 5B)
cd ../avr3a_allele
"$PYTHON" plot_avr3a.py                    # outputs/figS5_avr3a.pdf (Supplementary Figure S5)
cd ../..

echo "Done. Panels in each folder's outputs/, tables in Sup_material/."
echo "Tables S3, S4 and S6 are not code-generated (primer lists, alignment metrics,"
echo "effector annotation) and are not part of this repository; see README.md."
