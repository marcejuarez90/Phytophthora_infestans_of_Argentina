# Amplicon_pop — phylogeny from the 98-amplicon SNPs

Builds trees (Newick output) from the SNPs of the 98-amplicon panel (mostly
effector genes), straight from the multi-sample VCF produced with BWA-MEM +
bcftools: a neighbor-joining tree with bootstrap and a per-sample QC /
genotype-label check. Designed to be copied to
the server and run there; the trees are then drawn locally.

## Files

| File | What it is |
|---|---|
| `run_server.sh` | **the script to run on the server** (both analyses, NJ + QC) |
| `amplicon_nj_phylogeny.py` | the analysis script (Python ≥ 3.8 + numpy; bcftools is not needed) |
| `samples_ar_eu.txt` | sample list for the `ar_eu` analysis: 20 Argentine + 35 European reference isolates |
| `samples_exclude.txt` | samples left out of both analyses (currently `08_7054A`, see below) |
| `sample_genotypes.tsv` | genotype assigned to each sample: Argentine isolates (EU_2_A1); European references (from the reference table, with `14D11A`/`14D11B` as EU_8_A1); 54 Peruvian samples of the Kreuze collection (EC1, US_1, PE_3, PE_7, PE_7_Puno, with year, locality and host). Eleven samples of the `all` analysis still have no genotype (`97_28_1_1`, `HL_3b_19`, `HL_7a_40`, `POX067_*`) and are drawn grey as "Not assigned" |
| `plot_amplicon_tree.R` | draws a tree with every leaf colored by genotype (R + ape) |

## Running it on the server

```bash
cd Amplicon_pop
bash run_server.sh                     # defaults: VCF=/u03/mjuarez/marce/Pi_pop/total/totalpop_final_filter.vcf
                                       #           OUTDIR=results_excl_08_7054A
bash run_server.sh /path/to/other.vcf  my_outdir   # other VCF / output folder
NBOOT=100 THREADS=8 bash run_server.sh             # fewer bootstraps / more threads
```

`run_server.sh` checks that python + numpy are available, notes the versions
and the VCF checksum, runs both analyses and saves everything printed in
`<OUTDIR>/run.log`.

### The two analyses

| Analysis | Samples |
|---|---|
| `all` | every sample in the VCF (120), minus the excluded ones |
| `ar_eu` | the 55 samples in `samples_ar_eu.txt` (20 Argentine: `Ar20079` … `Ar21015`, `Biocca`, `Cba2`, `Garcia`, `PiNSL`, `Tafi`; 35 European references: `01_29` … `FAB22012E`), minus the excluded ones |

`97_28_1_1`, `HL_3b_19` and `HL_7a_40` are in the data directory but not in the
requested set; add them to `samples_ar_eu.txt` to include them. Lines starting
with `#` are ignored in both lists.

**Excluded sample `08_7054A`.** It is EU_2_A1 by SSR (`7054_A` in the SSR table),
but its amplicon profile is intermediate between EU_2_A1 and EU_13_A2 (mean
Prevosti distance 0.098 to EU_13_A2 isolates vs 0.130 to EU_2_A1 references),
and keeping it changed the topology of the NJ tree substantially. It is left out
until its origin is clarified; remove its line from `samples_exclude.txt` to
include it again.

Sample names are taken from the VCF header after removing any directory and
`.bam` / `.sorted` suffixes. The script prints the requested samples that are
**not** found in the VCF (and writes them to `ar_eu/ar_eu_samples_missing_from_vcf.txt`).

## Outputs

For each analysis, in `<OUTDIR>/<analysis>/`:

| File | Content |
|---|---|
| `<analysis>_nj.nwk` | **Newick** NJ tree with branch lengths; internal nodes carry the bootstrap support (%) |
| `<analysis>_sample_qc.tsv` | per sample: SNPs called, missing and heterozygous fractions, nearest neighbor, and the **label check** (see below) |
| `<analysis>_report.txt` | what was filtered and why, bootstrap summary, label-check warnings |
| `<analysis>_distance_matrix.tsv` | pairwise Prevosti distances |
| `<analysis>_samples_used.txt` | samples in the tree |

**Label check** (`--genotypes`, here `sample_genotypes.tsv`): for every labeled
sample, the mean distance to the other members of each genotype group; the
sample is flagged (`flag_closer_to_other_group = 1`, and a `LABEL CHECK` line in
the report) when it is closer, on average, to another group than to its own.
This is what would have flagged `08_7054A`. Groups of 2–3 samples give noisy
means, so read the flags together with the tree.

## Main options of `amplicon_nj_phylogeny.py` (`--help` lists all)

| Option | Default | Meaning |
|---|---|---|
| `--vcf` | — | multi-sample VCF (`.vcf` or `.vcf.gz`) |
| `--analyses` | `all,ar_eu` | which analyses to run |
| `--subset-list` | `samples_ar_eu.txt` | list for `ar_eu` |
| `--exclude` | none | sample(s) to leave out of **both** analyses: comma-separated names or a file with one name per line |
| `--genotypes` | none | sample/genotype table; adds the label check to the QC file |
| `--nboot` | 100 | bootstrap replicates (`run_server.sh` uses 1000) |
| `--max-site-missing` | 0.2 | drop sites with more missing calls than this (within each analysis) |
| `--max-sample-missing` | 0.5 | drop samples with more missing calls than this |
| `--min-gq` | 0 (off) | calls whose genotype quality (gap between the two lowest `PL`) is below this become missing, e.g. `--min-gq 20` |
| `--min-dp` | 0 (off) | needs per-sample `DP` (your VCF has `GT:PL` only, so it does nothing there) |
| `--keep-nonpass` | off | by default sites whose FILTER is not `PASS` or `.` are skipped |
| `--outgroup` | none | comma-separated sample(s) to root the NJ tree on; unrooted otherwise |
| `--seed` | 1 | random seed for the bootstrap |

## Method

- **Sites:** biallelic SNPs only (single-base REF and ALT), with FILTER `PASS`
  or `.`. Indels, multi-allelic and non-PASS sites are skipped and counted in
  the log.
- **Genotype:** the ALT-allele fraction computed from `GT`, whatever the ploidy
  (diploid `0/1` → 1/2; triploid `0/0/1` → 1/3, ...). Missing calls (`./.`) are
  missing.
- **Distance:** Prevosti's — mean over the sites called in both samples of
  |ALT fraction in sample *i* − ALT fraction in sample *j*|.
- **Filters,** applied separately in each analysis: sites above
  `--max-site-missing`; then samples above `--max-sample-missing`; then sites
  that are no longer polymorphic in that sample set.
- **NJ tree:** neighbor-joining (Saitou & Nei); negative branch lengths are set
  to 0; unrooted (trifurcating root) unless `--outgroup` is given. Support:
  SNP columns are resampled with replacement `--nboot` times, an NJ tree is
  built for each replicate and every internal edge of the reference tree gets
  the percentage of replicates that contain the same bipartition.

### Things to keep in mind
- **The VCFs are diploid calls** (`bcftools call -mv`, `FORMAT = GT:PL`): a
  triploid isolate shows up as `0/0`, `0/1` or `1/1` and the dosage (1/3 vs 2/3)
  is lost. Recovering it would need the read counts per allele
  (`bcftools mpileup -a FORMAT/AD` on the BAMs), which this script does not do.
- **Few SNPs.** The panel gives ~150–175 usable SNPs, so the deep nodes of NJ
  trees are unstable (removing one sample changed the topology a lot). Read the
  bootstrap before interpreting deep branches, and use the distance matrix for
  group-to-group statements.
- Tested on synthetic VCFs (clades recovered with 100% support, `.vcf.gz`,
  `.bam` names, exclusion lists, label check); the real run is yours.

## Drawing the trees (locally)

1. **Files to bring back** from `<OUTDIR>/` (or re-download the whole folder):
   `ar_eu/ar_eu_nj.nwk`, `ar_eu/ar_eu_sample_qc.tsv`,
   `ar_eu/ar_eu_report.txt`, the same three for `all/`, and `run.log`.
2. **Draw** (from the `Amplicon_pop` folder, with R and `ape`; `phangorn` for the
   midpoint root):

   ```bash
   Rscript plot_amplicon_tree.R results_excl_08_7054A/ar_eu/ar_eu_nj.nwk sample_genotypes.tsv results_excl_08_7054A/ar_eu/ar_eu_nj_tree
   ```

   Output: `<prefix>.pdf`, `.png`, `.svg`. Optional 4th argument: layout
   (`midpoint` by default, `unrooted` for a radial tree, `asis`); 5th: minimum
   support to print on the nodes (default 70).
   6th and 7th arguments: figure width and height in inches (default 9.5 x 0.17 per leaf + 1.8);
   8th: legend position, `right` (outside, default) or `bottomleft` (inside the empty area under the
   root, for small panels). Used to draw the panel of the supplementary figure:
   `Rscript plot_amplicon_tree.R <tree.nwk> sample_genotypes.tsv <prefix>_panel midpoint 70 6.5 8.4 bottomleft`.
3. **What it shows:** a colored dot on each leaf and the sample name in the same
   color, per genotype; Argentine isolates in bold; support on the internal nodes;
   a legend with the number of samples per genotype. EU_2_A1, EU_13_A2 and
   EU_8_A1 use the same colors as the SSR figures. Samples that are not in
   `sample_genotypes.tsv` (e.g. the other 65 samples of the `all` tree) are drawn
   grey as "Not assigned".
