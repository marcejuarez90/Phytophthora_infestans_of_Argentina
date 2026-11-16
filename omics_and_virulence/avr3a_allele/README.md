# AVR3a allele of PiNSL-19 (KI vs EM)

`avr3a_ki_em.py` genotypes the three residues that distinguish the two forms of
AVR3a (PITG_14371) from the whole-genome alignment:

| residue | 19 | 80 | 103 | |
|---|---|---|---|---|
| **AVR3a-KI** | C | K | I | recognized by R3a |
| **AVR3a-EM** | S | E | M | evades R3a |

R3a recognition is triggered by the presence of KI, so "this isolate escapes
R3a" is a negative claim. The script therefore reports the floor the claim
rests on as well as the call itself:

1. **codon composition** at each of the three sites, counted from whole reads,
   with and without a MAPQ >= 30 filter (AVR3a is single-copy and maps cleanly,
   so the two should agree — they are reported separately so that can be
   checked rather than assumed);
2. the **phased haplotype** of residues 80 and 103. Their codons start 69 bp
   apart and span 72 bp, so a 150 bp read covers both: this is what separates a
   homozygous EM isolate from a mixture that happens to look EM at each site
   taken alone;
3. the **error floor**, as the fraction of read bases differing from the
   reference over the whole AVR3a coding sequence;
4. a **one-sided 95% upper bound** on the frequency of an unobserved KI allele.

The residue numbering is verified, not trusted: the script translates the gene
model and stops unless residues 19, 80 and 103 each carry one of the two known
variants, printing the neighbourhood so an offset can be found instead of
silently counting the wrong codon.

## Running it

It reuses the exon- and strand-aware CDS machinery of
`../avrblb2_family/avrblb2_snp205.py`, so that coordinate code has one copy and
one place to be right. Copy **both** scripts into the folder that holds the BAM:

```bash
cp avr3a_ki_em.py ../avrblb2_family/avrblb2_snp205.py  /path/to/bam/
cd /path/to/bam
python3 avr3a_ki_em.py --self-test          # arithmetic only, no data needed
python3 avr3a_ki_em.py --bam pinsl.bam --fasta T30.fa --gtf T30.gtf \
    --out avr3a_kiem
```

Writes `avr3a_kiem/codon_pileup.tsv`, `residue_calls.tsv`, `haplotypes.tsv`,
`error_floor.tsv`, `detection_limit.tsv` and `avr3a_ki_em.{pdf,png}`.
Copy the TSVs back here to redraw the panel locally.
