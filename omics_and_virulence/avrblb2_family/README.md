# avrblb2_family — the AVR-blb2 paralog family in PiNSL-19

`AVRblb2` is not a single gene. The reference effector annotation
(`../inputs/SD_effector_list.xlsx`) labels **19** genes "Avrblb2 family, secreted
RxLR effector peptide" and assigns them to RxLRfam5 — the AVR-blb2 *extended*
family. Seven of those are the `Avrblb2` paralogs proper of Oliva et al. (2015,
MPMI 28:901-912): PITG_04085, PITG_04086, PITG_04090, PITG_18683, PITG_20300,
PITG_20301 and PITG_20303. They share >90% amino-acid identity with each other
(92.4% mean here) and 35-58% with the other twelve, and they form a clade of
their own in this panel's tree, which is why the figure bands them off.

`Rpi-blb2` recognition turns on one residue, the position equivalent to 69 in
PITG_20300 (Oh et al. 2009: Ala/Ile/Val trigger the HR, Phe evades it). That call
is only defined inside the banded clade and is **not** extended to the rest: the
two paralogs that carry Phe-69 and therefore break Rpi-blb2 resistance,
PITG_20301 and PITG_20303, are simply marked in bold red.

This folder answers, for each of the 19:

| | question | source |
|---|---|---|
| **residue 69** | which residue does the reference carry at that position (table only) | `../inputs/paralogs_protein{,_aligned}.fasta` |
| **WGS coverage** | is the paralog present in the PiNSL-19 genome | `inputs/avrblb2_gene_coverage_profile_mapq30.csv` |
| **expression** | is it transcribed during infection (table only; the
curve is plotted in the notebook) | `../inputs/TPM_matrix.tsv` |

```
python avrblb2_family.py
```

writes `outputs/avrblb2_family.csv` (one row per paralog),
`outputs/present_paralogs.txt` and `outputs/avrblb2_family.{pdf,png}` — the
neighbor-joining tree with the whole-genome coverage profile.

Expression is **not** drawn here. The infection time course of the paralogs this
panel calls present is plotted in `../omics_and_virulence.ipynb`, next to the
other AVR loci, so that a single axis carries both; that cell re-derives the
presence call from this folder's coverage profile with the same rule, so the two
panels cannot disagree.

## Presence call

`genome_call` in the output table is three-level, from the per-base whole-genome
depth over the CDS: `present` (at least one read over *every* coding base),
`partial` (some of the CDS covered) or `absent` (no read anywhere in the CDS).
On this isolate that is 14 present, 1 partial (PITG_04089, 48% of its CDS) and 4
absent. The bar is deliberately strict, because depth alone says little here:
every member of the family sits below the genomic median gene depth (41x), not
because copies are missing but because near-identical copies share their reads.

## Why the residue is read off the reference and not off the reads

Several of these paralogs are nearly identical, so a read cannot be assigned to
one of them with any confidence. Calling the allele at residue 69 from the
isolate's own reads therefore reads whatever mapped to that locus, including
reads that belong to a relative.

PITG_20301 and PITG_20303 are the two paralogs the reference annotates with the
evasion residue, and they show the problem directly: **no whole-genome read maps
to either of them**, yet they collect 346 and 634 RNA-seq reads, of which 53% and
56% are multi-mapped (`frac_NHgt1` in the output table). The DNA says the genes
are not in this isolate; the RNA signal is cross-mapping from the expressed
paralogs of the family.

So the question the family can actually answer is presence/absence, and the
script keeps the two measurements that decide it apart: total depth and the
depth of uniquely mapping reads (MAPQ >= 30). A paralog is called absent only
when neither returns a read along its entire coding sequence. On this isolate
that is PITG_09632, PITG_18683, PITG_20301 and PITG_20303 — the two evasion
paralogs among them.

Expression is reported as the gene-level mean TPM per condition, the same
quantity as the curated effector panel, not as reads over any single codon.

## Coverage input

`inputs/avrblb2_gene_coverage_profile_mapq30.csv` is produced by
`extract_avrblb2_coverage.py`, which needs the indexed BAM and the T30-4
annotation, so it runs where the alignment lives and is not part of the normal
build:

```
python3 extract_avrblb2_coverage.py --bam pinsl.bam \
    --annotation P_infestans_T30-4.gtf \
    --out inputs/avrblb2_gene_coverage_profile_mapq30.csv
```

Its output is committed here so `avrblb2_family.py` runs without the BAM. The
columns are those of `../inputs/avr_gene_coverage_profile_mapq30.csv` (`Gene_ID,
Gene_name, Category, chrom, strand, genomic_pos, gene_pos_bp, in_cds, depth,
depth_mapq30`): CDS +/- 300 bp, `gene_pos_bp` relative to the CDS start in the
gene's own orientation, and both total and MAPQ >= 30 depth.

## Requirements

Those of `../README.md` (Python with pandas, numpy, matplotlib). The
neighbor-joining is implemented here, so no phylogenetics package is needed; the
protein alignment was built with MAFFT (`--auto`) and is committed.
