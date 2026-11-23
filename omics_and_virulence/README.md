# omics_and_virulence — genome, RNA-seq and infection-assay panels

Genome architecture and effector transcriptional landscape of PiNSL-19, the potato leaf
infection assay and the *N. benthamiana* transgenic-line assay. `../SSR_populations/`
covers the population-genetics figures.

| Result | Built by | Written to |
|---|---|---|
| Genome and effector-expression panels (Figure 3, Supplementary Figure S3) | notebook | `outputs/fig3/`, `outputs/fig_sup3/`, `outputs/fig_sup4/` |
| AVR read-depth profiles (Figure 3, panel D) and potato-assay panels (Figure 4) | notebook | `outputs/fig4/` |
| *N. benthamiana* assay (Figure 5, panel C) | notebook | `outputs/fig5/` |
| AVR-blb2 paralog family (Figure 5) | `avrblb2_family/` (own README) | `avrblb2_family/outputs/` |
| RT-qPCR vs RNA-seq regression (Supplementary Figure S4) | `qPCR_validation/` (own README) | `qPCR_validation/` |
| Supplementary Table S5 | `make_table_s5.py` | `../Sup_material/TableS5_gene_expression_PiNSL19.xlsx` |

Every panel is written on its own. How panels are arranged into the figures of the manuscript
is a presentation choice and is not part of this repository.

The notebook reads its read-mapping QC from CSVs in `inputs/`. Those are produced on the compute
cluster, from the BAM files, by `effector_mapping_qc/` and `avrblb2_pipeline/`, which need the
alignments and are not part of this repository; the sections below name their scripts so the
workflow can be followed, and every cell that uses them says what it does when the CSV is
absent.

Figure 3: **A** heterozygous allele-frequency histogram; **B** gene copy number by effector
family; **C** curated effector panel, mean log2(TPM + 1), rows clustered (PITG_15278 and
AVR-cap1b are not shown; AVR-smira2, PITG_07558, is labeled AVR8); **D** WGS read depth over
AVR1, AVR3b, AVR4 (PITG_07387) and AVR3a (shared y axis 0-40).

Figure 4: **A** potato leaf infection assay (Spunta / Innovator photos and % leaf area
affected); **B** map of the Argentine isolates and their phenotype on Innovator
(`../SSR_populations/outputs/s2/mapa_argentina_innovator.png`, drawn by the last section of
`../SSR_populations/haplotype_network_v4.R`); **C** AVR1/AVR2/AVR3a/AVR3b mean log2(TPM + 1).

## Layout

```
omics_and_virulence/
├── omics_and_virulence.ipynb   # the notebook -- run top to bottom
├── README.md
├── inputs/                     # all data the notebook reads
├── qPCR_validation/            # Supplementary Figure S4 (R script + data, own README)
└── outputs/                    # one file per panel, plus PNG/SVG previews
    ├── fig3/  fig4/  fig5/
    ├── fig_sup3/               # Supplementary Figure S3 (folder names kept from an earlier numbering)
    └── fig_sup4/               # mapping-QC standalone outputs, also Supplementary Figure S3
```

## Run

```
jupyter nbconvert --to notebook --execute --inplace omics_and_virulence.ipynb   # ~4 min
```

from this folder (the notebook uses paths relative to it). It refreshes `outputs/`. Run
`../SSR_populations/haplotype_network_v4.R` first if the Argentina map has to be regenerated.

## Notebook sections (details)

The notes below describe each analysis section. Panel letters in them refer to the figures of
the manuscript, for orientation only: the notebook writes each panel on its own.


- **Panel A** — heterozygous SNP allele-frequency distribution (ploidy
  evidence). Reads nQuire (`lrdmodel`) output and plots the bimodal
  allele-frequency histogram with reference lines at 1/3 and 2/3
  (triploid signature). Plotted in vertical/portrait format (allele
  frequency on the Y axis, number of sites on the X axis).
- **Panel B** — gene copy number by effector family (RxLR, Crinkler,
  apoplastic — plus ubiquitin and a genome-wide total, once available; see
  below), restricted to genes with reliable whole-genome-sequencing mapping
  (`mean_mapq > 30`). This is the figure's original copy-number panel,
  brought back once a mean-MAPQ-filtered version made it defensible in the
  main text; the *unfiltered* version (all genes, any mapping quality) now
  lives instead in the "CNV + volcano + Venn" Supplementary Figure below,
  shown right next to this filtered one for direct comparison. Merges
  `gene_ploidy_v2.csv` with the WGS mapping-QC per-gene table and plots the
  same style of copy-number stacked bar as the unfiltered version, with
  each bar's *n* showing genes retained over genes tested.

  **Bars shown are dynamic, not hardcoded**: RXLR/CRINKLER/APOPLASTIC are
  always present; the load cell looks for an *extended* QC file
  (`inputs/effector_mapping_qc_wgs_per_gene_ext.csv`) that additionally
  covers the rest of the genome (built via `effector_mapping_qc/
  effector_mapping_qc.py --all-genes`, which tags every non-effector gene
  `OTHER` straight from the GTF -- no separate CSV to build; see that
  folder's README, "Genome-wide run") and, if found and an `OTHER` category
  is present, adds a `total` bar automatically (`FCNV_BAR_ORDER` is built at
  run time from whichever categories are actually present in the loaded QC
  file; a `UBIQUITIN` category is also supported and gets its own bar if
  present, for anyone using the older `--extra-genes` path instead of
  `--all-genes`); if only the base effectors-only file
  (`effector_mapping_qc_wgs_per_gene.csv`) is present, it falls back to the
  original three bars with a printed note, so this section works either way.
  The `total` bar, when present, is a genuine genome-wide aggregate (built
  the same way as the unfiltered panel's `total` bar) rather than an
  effectors-only stand-in. Figure width and tick-label rotation adjust
  automatically for 3 vs. 4/5 bars, in both the standalone plot and its
  embedding elsewhere.

  **Result (real data, effectors-only run): RXLR (432/557, 78% retained) and
  apoplastic (248/303, 82%) are essentially unchanged by the filter**
  (3-copy fraction 72.2%→77.5% and 73.9%→76.6% respectively) — evidence
  their unfiltered calls are already trustworthy — **while Crinkler retains
  only 30/128 genes (23%), and among those its copy-number skew attenuates**
  (2-copy fraction 49.2%→30.0%; 3-copy fraction 34.4%→50.0%), moving
  noticeably closer to the ~3-copy majority seen in the other two families
  (though n=30 is small, so this reads as attenuation, not proof the skew
  fully disappears). The genome-wide `total` result is pending the
  `--all-genes` cluster run (not yet completed with real data).

  `FCNV_MAPQ_THRESHOLD` and the `FCNV_*` styling variables are exposed as
  config at the top of the section. Self-contained; requires
  `inputs/effector_mapping_qc_wgs_per_gene.csv` (see the "Effector mapping
  QC" sections and `effector_mapping_qc/README.md` for how it's
  produced) alongside `gene_ploidy_v2.csv`, plus optionally
  `inputs/effector_mapping_qc_wgs_per_gene_ext.csv` for the genome-wide
  `total` bar.
- **Panel C** — heatmap of the top 100 most upregulated effectors during
  infection. Two-step criterion: first, restrict to annotated effectors
  significantly induced (log2FC > 2, adjusted p-value < 0.05 — same
  thresholds used for the DE stats reported in Results) at one or more
  of the 24/48/72 hpi timepoints vs. in vitro mycelium (160 effectors
  meet this bar: 105 RxLR, 55 apoplastic, 0 Crinkler); then, of those,
  rank by degree of upregulation (max log2FC across the three
  timepoints) and keep the top 100 (79 RxLR, 21 apoplastic, 0 Crinkler
  — Crinkler drops out entirely at the filtering step, before ranking).
  This replaces an earlier version that picked the "top 100 by lowest
  p-value at any timepoint, regardless of sign" — that criterion let a
  handful of confidently *repressed* genes into the panel and, combined
  with a sequential (viridis) colormap with no true zero point, made
  the panel read as more downregulation-heavy than the underlying data
  actually is; the filter-then-rank-then-cap logic fixed that. Rows are
  hierarchically clustered (dendrogram) by Euclidean distance between
  each gene's 3-timepoint log2FC profile with average linkage; both the
  distance metric (`CLUSTER_METRIC`) and linkage method (`CLUSTER_METHOD`)
  are exposed as config variables at the top of the section if the
  grouping needs to be revisited (e.g. `CLUSTER_METRIC = "correlation"`
  groups by the *shape* of the temporal profile instead of the raw
  magnitude of the fold-change). Columns are fixed in temporal order.
  Color scale is a **continuous viridis colormap**: an earlier version
  used a custom diverging violet-white-green map with a true zero point,
  built for when this panel could contain real downregulation; now that
  the panel is pre-filtered to induced genes only, the few remaining
  slightly-negative cells are negligible (down to about -0.1 log2FC,
  against inductions up to +14.6), so a plain continuous scale
  (`HEATMAP_CMAP`/`HEATMAP_NORM`) over the actual data range is simpler
  and just as informative. Row labels are colored by effector type:
  RXLR = red, Crinkler = blue, every other annotated class unified as
  "Apoplastic" = gold. Includes a printed count of how many induced
  effectors fall into each category.

  (An IGV screenshot of the AVR1 locus, PITG_16663, was considered early
  on as a possible extra panel; it isn't code-generated and isn't part
  of this notebook's A/B/C panel lettering.)
- **Expression-profile clusters among non-Crinkler effectors (exploratory
  — not part of any manuscript figure)** — looks for discrete kinetic
  "programs" (joint activation or silencing) among RXLR/apoplastic
  effectors (**Crinkler excluded** — see the Crinkler mapping-quality
  Supplementary Figure below for why its expression calls aren't trusted
  enough to include in a clustering analysis meant to find genuine
  biological groups). Criterion is looser than Panel C's: any gene with
  |log2FC| > 1 and adjusted p < 0.05 at ≥1 of the three in-planta
  timepoints (not just log2FC > 2 and not just induction — this is
  deliberately symmetric so repression is caught too, since repression
  signatures are typically smaller-magnitude than the strongest induction
  and log2FC > 2 alone would filter them all out). 308 non-Crinkler
  effectors have complete 3-timepoint data; 218 meet this criterion.

  **Two-step, direction + shape clustering** (revised from an earlier
  single-pass version): a first version ran k-means directly on raw
  log2FC values, which groups genes by *how induced/repressed* they are,
  not by *when* they turn on or off — concretely, AVR-blb2 (rising sharply
  24→48 hpi) and AVR2 (already near-plateau at 24 hpi) landed in the same
  "very strong induction" bucket purely by magnitude, while AVR3a (same
  rising shape as AVR-blb2, lower amplitude) was split into a different
  one. The fix clusters on each gene's **shape** instead: **(1)** split
  by overall direction first (induced: mean log2FC > 0, 174 genes;
  repressed: mean log2FC ≤ 0, 44 genes), **(2)** within each direction,
  z-score each gene's own 3-timepoint profile (subtract its own mean,
  divide by its own SD, to strip out amplitude) and run a dependency-free
  k-means on the z-scored shapes (`K_INDUCED = 3`, `K_REPRESSED = 2`, best
  of `N_KMEANS_SEEDS = 15` random initializations by inertia — same
  from-scratch-implementation convention as the AVR-blb2 NJ tree; k chosen
  by inspecting inertia/interpretability at k = 2-4 per direction).
  Clusters within each direction are relabeled in ascending order of the
  shape centroid's first-timepoint (24 hpi) z-score, which corresponds to
  earliest-peaking to latest-rising for the induced side. Produces two
  plots side by side (mean raw log2FC profile per program ± SEM; program
  composition by RXLR vs. Apoplastic), plus a lookup printing which
  program each AVR gene of interest falls into.

  **Line colors are a single viridis scale by overall magnitude**, not two
  separate hue families: each program's mean raw profile (averaged across
  the 3 timepoints) gives it a score, the 5 programs are ranked by that
  score, and `PROGRAM_COLORS` samples 5 evenly-spaced points off viridis in
  that rank order -- darkest violet for the most repressed program,
  yellow for the most induced, with the 3 induced programs' relative
  strength landing them in the teal-to-yellow-green stretch in between.
  This replaces an earlier green-for-induced/red-for-repressed scheme and
  is computed once (after `cluster_centroids_raw` is fully populated) so
  it's shared by this section's own plot, its embedding in both Figure 3
  wrap-up layouts below, and stays a one-line change (`_cl_program_score`)
  if a different ranking is ever wanted.

  **Real result:** five programs. **Induced -- Early peak** (n=82,
  5.94/4.99/3.29, 68% RXLR) is already near-maximal by 24 hpi and declines
  through 72 hpi. **Induced -- Mid-course peak** (n=36, 3.97/5.78/3.96,
  56% RXLR) peaks at 48 hpi and drops back by 72 hpi — **both AVR2 loci
  (PITG_08943, PITG_22870) fall here.** **Induced -- Late, sustained**
  (n=56, 1.15/5.00/4.54, 62% RXLR) is low at 24 hpi, jumps by 48 hpi, and
  stays high through 72 hpi — **AVR3a and both AVR-blb2 paralogs
  (PITG_20300, PITG_04090) fall here**, resolving the exact issue that
  motivated the revision: under magnitude-based clustering these genes
  were split apart by amplitude despite sharing this delayed-onset,
  sustained kinetic. **Repressed -- Early, transient** (n=29,
  -3.39/-1.20/-0.10) is sharply repressed at 24 hpi and largely recovers
  by 72 hpi, 100% apoplastic (degradative enzymes: hydrolases,
  serine/cysteine proteases). **Repressed -- Delayed, progressive** (n=15,
  -1.27/-4.15/-2.74) is mild at 24 hpi, deepest at 48 hpi, still
  substantial at 72 hpi, 93% apoplastic (pectin lyase, CBEL-elicitin,
  cysteine protease, glucanase). Both repression programs are essentially
  RXLR-free, consistent with silencing being an apoplastic-specific
  phenomenon, but splitting confirms that response itself isn't uniform;
  on the induced side, RXLR effectors are spread across all three kinetics
  (56-68%) rather than concentrated in one.

  This section doesn't feed a lettered panel in the current figure set
  (Figure 3 recovered its original allele-frequency / copy-number /
  heatmap layout once the mean-MAPQ-filtered copy-number panel became
  available — see Panel B above); it's kept as a standalone analysis in
  case it's useful for the Discussion or a future figure. `K_INDUCED`,
  `K_REPRESSED`, `N_KMEANS_SEEDS`,
  `CLUSTER_LFC_THRESHOLD`/`CLUSTER_PVAL_THRESHOLD`, and
  `INDUCED_COLORS`/`REPRESSED_COLORS` are exposed as config variables at
  the top of this section (prefixed `CLUSTER_`/`_cl_` to avoid clashing
  with Panel C's own `LFC_THRESHOLD`/`PVAL_THRESHOLD`/etc., since both
  sections reload independent copies of the same input files under
  different names). Self-contained. The composition bar chart's bar
  *heights* are raw gene counts (not normalized to 100%), so a small
  program (e.g. n=15) doesn't visually loom as large as a big one (n=82)
  just because both would otherwise reach the top of the axis — an earlier
  version normalized to % of program instead. Saved as
  `outputs/fig3/cluster_profiles_and_composition.png`.

  (An earlier iteration of this notebook also tried extending this
  clustering genome-wide, tagging non-effector genes `OTHER`, plus a
  per-gene-lines variant of the plot. Both were tried and dropped — the
  per-gene spaghetti-plot version in particular didn't read well with
  ~250 overlapping lines per panel — in favor of the pie-chart contrast
  below, which answers the same "how much do effectors matter genome-wide"
  question more legibly.)
- **Gene category composition: total genome vs. expressed vs. upregulated
  (exploratory — not part of any manuscript figure)** — a 100%-stacked
  horizontal bar chart (RXLR / Crinkler / Apoplastic / Other), same
  category breakdown across three increasingly restrictive gene sets, to
  contrast the whole transcriptome against the genes that actually matter
  during infection: **(1) Total genome** — every gene in
  `tabla_LFC_pvalor.tsv`, regardless of expression. **(2) Expressed** —
  genes with complete (non-missing) log2FC/p-value at all three in-planta
  timepoints (a gene only gets a value if DESeq2 could test it — used as a
  practical proxy for "expressed enough to analyze," not a TPM-based
  cutoff). **(3) Upregulated during infection** — Panel C's criterion
  (log2FC > 2, adjusted p < 0.05, at ≥1 timepoint). Unlike the
  effector-clustering section above, Crinkler is included here (no
  reliability concern — this is composition, not kinetic clustering).

  **Stacked bars instead of three separate pies (an earlier version):**
  all three gene sets share a common 0–100% baseline, so a category's
  share and how it shifts across gene sets are direct length comparisons
  rather than angle judgments across three circles, and the fold-enrichment
  figure (previously reported only in this section's result text) is now
  annotated directly on the plot. Built by
  `plot_category_composition_bars()`, defined once in the load cell and
  reused both by this section's own standalone plot and its embedding in
  the Figure 3 wrap-up layout below, so both stay visually identical.
  Segments below 8% of their bar don't get an inline label (RXLR/Crinkler/
  Apoplastic are all single digits in the Total genome/Expressed rows —
  their smallness next to the Upregulated row *is* the point, and a label
  that size either overflows its own segment or gets clipped by the axis).

  **Result (real data):** total genome (n=18,179): RXLR 3.06% (557),
  Crinkler 0.70% (128), Apoplastic 1.67% (303), Other 94.57% (17,191).
  Upregulated (n=496): RXLR 21.17% (105, a ~6.9x enrichment), Apoplastic
  11.09% (55, ~6.6x), Crinkler 0.00% (0/128, complete depletion, consistent
  with "0/128 Crinkler genes ever induced"), Other 67.74% (336, still the
  majority by raw count but only ~0.72x its genome share). The "expressed"
  bar (n=10,213) is a useful middle step: RXLR actually *drops* to 1.45%
  there (roughly half its genome share) while Apoplastic stays close to
  flat and Crinkler ticks up slightly — so the enrichment is specific to
  the upregulation step, not a byproduct of which genes happen to be
  expressed at all. `PIE_LFC_THRESHOLD`/`PIE_PVAL_THRESHOLD` and
  `PIE_CATEGORY_COLORS` are exposed as config at the top of the section
  (variable names kept `PIE_*`/`pie_*` even though the plot itself is no
  longer a pie chart, to avoid churning every reference across this
  section and its wrap-up embedding). Self-contained. Saved as
  `outputs/fig3/gene_category_composition_total_expressed_upregulated.png`.
- **Supplementary Figure — CNV + volcano + Venn** — consolidates the
  *unfiltered* gene copy-number-by-family view (all genes, any mapping
  quality — Figure 3's Panel B is the same breakdown restricted to
  well-mapped genes only, for direct comparison) together with the
  volcano/Venn overlap analysis, because the CNV skew and the expression
  pattern shown here need the same mapping-quality caveat (next
  Supplementary Figure), and it felt better to keep both readings of the
  copy-number data (raw and filtered) easy to compare rather than only
  ever showing one.
  - **Panel A** — gene copy-number distribution per effector family (RxLR,
    Crinkler, apoplastic, ubiquitin, plus an aggregate total), as a stacked
    bar chart of copy-number proportions (2/3/4/>=5/NA), unfiltered by
    mapping quality. Reads `gene_ploidy_v2.csv`.
  - **Panels B-D** — volcano plots (one per timepoint, 24/48/72 hpi), grey
    background for all non-effector genes, effector families overlaid in
    the same color code as panel A (RXLR = red, Crinkler = blue, apoplastic
    = gold), with reference lines at the log2FC and adjusted-p thresholds.
    Self-contained (reloads its own copies of `tabla_LFC_pvalor.tsv` and
    `tabla_efectores.xlsx`, independent of Panel C of Figure 3, the
    heatmap).
  - **Panel E** — one large 3-set Venn diagram (24/48/72 hpi) for all
    upregulated genes.
  - **Panels F-I** — a 2x2 grid of smaller Venn diagrams applying the same
    24/48/72 hpi comparison separately to all effectors combined (F), RxLR
    (G), Crinkler (H) and apoplastic (I) — all drawn with a hand-rolled,
    dependency-free `venn3` implementation (no `matplotlib_venn` available
    in this environment), consistent with the notebook's from-scratch
    approach elsewhere (NJ tree). `draw_venn3()` takes optional
    `fontsize_count`/`fontsize_label`/`fontsize_title` arguments so the 2x2
    grid's smaller panels can use reduced font sizes proportional to their
    smaller area.
  - The printed core-set sizes (153 genes / 81 effectors upregulated at
    all three timepoints, split 61 RxLR + 20 apoplastic + 0 Crinkler)
    match the values reported in Results and are the same gene sets used
    for the g:Profiler functional-enrichment analysis described in
    Methods. The Crinkler-family Venn (H) is all zeros at every region.

  Final layout: CNV bar chart on top left (**A**), volcano plots stacked
  below it (24/48/72 hpi, panels **B**-**D**); on the right, the
  total-genes Venn on top (panel **E**, large) and the 2x2 per-family grid
  below it (panels **F**-**I**). Saved as
  `outputs/fig_sup4/fig_sup4_final.pdf`/`.svg`/`.png`.
- **Supplementary Figure — Crinkler mapping-quality diagnostics** —
  companion figure to the one above: that figure shows *what* Crinkler's
  copy-number and expression pattern looks like; this one shows the direct
  evidence for *why* those specific numbers need a caveat. RxLR and
  apoplastic are not subject to this concern.
  - **Panels A-B** — genome architecture, the manuscript's original,
    indirect technical-artifact check: contig GC content vs. length (A)
    and "two-speed genome" intergenic distances (B), both for RxLR
    (left column) and Crinkler (right column), colored by each gene's own
    copy number (NA / 2 / 3 / 4 / >=5). These are proxies for
    repetitiveness/mappability, not a direct measurement, and by
    themselves don't distinguish Crinkler from RxLR/apoplastic.
  - **Panels C-D** — direct RNA-seq mapping QC by effector category: mean
    MAPQ (C) and fraction of reads with MAPQ = 0 (D). Requires
    `inputs/effector_mapping_qc_per_gene.csv`, produced by running
    `effector_mapping_qc.py` on the cluster (see
    `effector_mapping_qc/README.md` for the full workflow and
    `run_effector_qc.sh` to launch it); shows "not available" in these two
    panels if the file isn't there yet. **Result (run on real cluster
    data, 2026-09-11): Crinkler maps substantially worse** — median MAPQ
    29.4 vs. 60 for RXLR/apoplastic, median fraction multi-mapped (NH > 1)
    51.5% vs. ~0% (both p < 0.001, this second metric not plotted as its
    own panel here but reported in Methods/Results) — while coverage
    breadth/depth is as high or higher, meaning that coverage is built
    from ambiguously-mapped reads shared across paralogs rather than
    uniquely assigned ones. The "0/128 Crinkler genes ever induced"
    result (Figure 3, Panel C) should be read with this caveat.
  - **Panels E-F** — the same two metrics (mean MAPQ, fraction MAPQ = 0)
    for the whole-genome sequencing (WGS) reads that the copy-number calls
    (Panel A of the Supplementary Figure above; `gene_ploidy_v2.csv`,
    Knaus et al. 2020 method) are built from. Requires
    `inputs/effector_mapping_qc_wgs_per_gene.csv` (same script, run with
    `--only pinsl.bam` — `run_effector_qc_wgs.sh`). **Result (real cluster
    data, 2026-09-11): same pattern as RNA-seq, if anything stronger** —
    median MAPQ 16.3 for Crinkler vs. 58.8 (RXLR) / 60.0 (apoplastic),
    median fraction of MAPQ = 0 reads 48.75% vs. 0% for both others (all
    p < 1e-30); coverage breadth ~1.0 across all three, so Crinkler isn't
    under-covered, that coverage is just built from ambiguous reads. Note:
    the WGS BAM was aligned with BWA-MEM, which doesn't set the NH:i tag
    by default (unlike the STAR/HISAT2 RNA-seq BAMs), so
    `frac_NH1`/`frac_NHgt1` are undefined for this run — `frac_mapq0`
    plays the equivalent multi-mapping-signature role here, and
    `plot_qc_panel()` shows "not available" instead of crashing on the
    undefined NH columns (a real bug hit and fixed during this analysis).
  - **Panels G-H** — within the WGS run, each gene's own mapping quality
    (mean MAPQ, G; fraction of MAPQ = 0 reads, H) plotted against its
    assigned copy number, points colored by category
    (`correlate_with_copy_number()` / `plot_copy_number_scatter()`).
    **Result: within Crinkler specifically, genes with better mapping
    quality tend to be assigned higher copy number** (mean MAPQ vs. copy
    number, Spearman rho = 0.20, p = 0.024; fraction of MAPQ = 0 reads vs.
    copy number, rho = -0.19, p = 0.036) — a dose-response relationship
    not seen in RxLR or apoplastic (p > 0.07 for the same comparisons).
    Modest effect sizes that don't survive strict correction for the
    number of metrics tested, but directionally consistent with the much
    stronger between-family differences in panels C-F.

  Reuses `compare_categories()`/`plot_qc_panel()`/`correlate_with_copy_number()`/
  `plot_copy_number_scatter()` from the standalone
  `effector_mapping_qc_compare.py` script, so this notebook section and
  that script stay in sync; each optional panel shows "not available" if
  its required cluster-side CSV hasn't been generated/copied into
  `inputs/` yet, so this section is safe to leave in the notebook before
  that's been done. Saved as `outputs/fig_sup3/fig_sup3_final.pdf`/`.svg`.
- **AVR1/AVR2/AVR3a/AVR3b mean TPM across timepoints (exploratory, for
  Figure 4)** — line plot of mean log2(TPM + 1) (Myc, 24, 48, 72 hpi) for
  the four effectors tied to the R1/R2/R3a/R3b specificities tested on
  Innovator, one color per gene (viridis), styled to match the
  potato-assay boxplot. Each replicate is log2(TPM + 1)-transformed
  individually (pseudocount = 1) before averaging, so the shaded ± 1 SD
  band is computed in log2 space rather than log-transforming a raw-TPM
  mean/SD. `TPM_efectores.tsv` has two distinct PITG loci both annotated
  "AVR2" (PITG_08943 and PITG_22870); both are plotted separately rather
  than merged, since which one is the "real" AVR2 is an open question
  flagged elsewhere in the manuscript. Lines that are flat at 0 (AVR1 and
  AVR3b — TPM = 0 at every timepoint, so log2(TPM + 1) = 0 for both) are
  nudged apart by a small fixed +0.3 (log2 units) offset purely so they
  don't visually overlap — noted directly on the plot.
- **Figure 4 — potato leaf infection assay (Spunta / Innovator)** —
  boxplot of % leaf area affected per isolate x variety, with
  representative leaf photos placed above (Spunta) and below
  (Innovator) each isolate column. Row labels include the R genes each
  variety expresses ("Spunta / R1", "Innovator / R1 R2 R3a R3b"). Points
  below a 5%-lesion threshold in Innovator are flagged with a red X as a
  hypersensitive response rather than a compatible lesion. Extracted
  from an earlier version of this notebook; only the plotting step is
  reproduced here (the upstream VIA/pliman image-annotation pipeline
  that produces `areas_por_hoja_total.csv` from the raw leaf photos and
  VGG Image Annotator JSON files is not included).
- **Figure 5 — *N. benthamiana* infiltration assay** — reads `WL_Nb.tsv`
  (one row per infiltrated leaf: `variedad` = transgenic line, `aislamiento`
  = *P. infestans* isolate, `healthy`/`disease` = percent leaf area,
  `HR` = hypersensitive response Y/N), drops the `GS4` line (not part of
  this figure), and plots percent leaf area affected (`disease`) as
  boxplots grouped by transgenic line (X axis, labeled "Transgenic line":
  EV, R1, R3a, Rpi-blb2, Rpi-amr1, Rpi-amr3 — the raw `blb2`/`amr1`/`amr3`
  values are relabeled with the `Rpi-` prefix for display, matching the
  labeling already used in `Nbs.png`), stacked one row per isolate
  (PiNSL-19, Cba2, Ar20079, Tafi). Unlike Figure 4, color here doesn't
  encode the categorical grouping variable (transgenic line) — each
  isolate's row gets a single color from `ISOLATE_COLOR_GROUP_5`, which
  groups isolates by virulence phenotype on Innovator rather than
  alternating row to row: PiNSL-19 and Cba2 share one color, Ar20079 and
  Tafi the other. Points are colored dots where `HR = NO`, red crosses
  where `HR = YES` (a real
  per-leaf HR call, unlike Figure 4's distance-based 5% proxy). Saved as
  `outputs/fig5/fig5_lesion_boxplot.pdf`/`.svg`.
- **AVR1/AVR3a/AVR-blb2/AVR-amr1/AVR-amr3 mean log2(TPM + 1) across
  timepoints (exploratory; not a panel of the published Figure 5, whose
  expression data are in Figure 3C)** — same source and log2(TPM + 1) transform as
  the Figure 4 TPM panel (per-replicate transform, then averaged), for the
  five effectors that correspond to the transgenic lines tested in
  Figure 5 (*AVR1*↔R1, *AVR3a*↔R3a, *AVR-blb2*↔blb2, *AVR-amr1*↔amr1,
  *AVR-amr3*↔amr3). *AVR-blb2* here is represented solely by its paralog
  PITG_04090 (annotated "Avrblb2 family, secreted RxLR effector peptide,
  putative" in `tabla_efectores.xlsx`) rather than PITG_20300, since
  full paralog-level detail for the whole AVR-blb2 family (all 19
  members, including PITG_20300 and PITG_04090) is shown separately in
  the AVR-blb2 family analysis. PITG_04090 is not in `TPM_efectores.tsv`,
  so it's pulled from the full genome-wide `TPM_matrix.tsv` instead (see
  Input data below) and reshaped to match. Only *AVR1* is flat at 0
  across all timepoints here, so no visual offset or on-plot note is
  needed. Saved as `outputs/fig5/AVR1_AVR3a_blb2_amr1_amr3_TPM.png`.
- Panel B's standalone plot cell uses a section-local name for its tick
  labels (`fcnv_tick_labels`) rather than the bare `tick_labels` set by
  Panel A, so running the notebook top to bottom cannot clobber it.
