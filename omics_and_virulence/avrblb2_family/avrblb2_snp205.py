#!/usr/bin/env python3
# ══════════════════════════════════════════════════════════════════════════════
#  avrblb2_snp205.py — is the Rpi-blb2 evasion allele anywhere in this isolate?
#
#  Cluster-side companion of avrblb2_family.py. Where that script asks "does a
#  read map to PITG_20301?", this one asks the question that does not depend on
#  assigning a read to a paralog at all:
#
#      does ANY read in the library carry the Phe-69 codon of Avrblb2?
#
#  Why that is the right question
#  ------------------------------
#  The Avrblb2 paralogs are too similar for a read to be assigned to one of them
#  with confidence, and a paralog carried by the isolate but missing from the
#  T30-4 assembly would have its reads pile onto the closest reference relative.
#  Both problems disappear if the reads are pooled: a read from any Phe-69 gene,
#  mapped anywhere in the family, still shows T at the first base of codon 69.
#
#  This is the same measurement Oliva et al. (2015, MPMI 28:901-912) make by
#  pyrosequencing, where they call it SNP 205 and report the proportion of each
#  nucleotide in the genome pool:
#
#      205A -> Ile-69      205G -> Ala-69 or Val-69      205T -> Phe-69 (evasion)
#      205C -> never observed in 352 isolates
#
#  The marker is valid as a single base because the codons for the four natural
#  residue-69 variants are invariant across P. infestans strains, and because the
#  seven Avrblb2 paralogs are all 100 aa and align without gaps: residue 69 is
#  residue 69, and nucleotide 205 is nucleotide 205, in every one of them.
#
#  What the script reports
#  -----------------------
#    1. per-paralog pileup of codon 69 (all three bases), with and without the
#       MAPQ >= 30 filter -- the filter is reported, never applied, because the
#       multi-mapping reads are exactly the informative ones here;
#    2. the pooled SNP-205 proportions and the resulting genotype (205GAT,
#       205GA, 205GT ...), directly comparable with Oliva et al.'s Figure 5;
#    3. the detection limit. The claim is negative, so it is only worth as much
#       as the floor it sits on: at every coding position of the seven paralogs
#       the script counts reads carrying a base other than the reference one,
#       which measures sequencing error and cross-paralog mismatch together, in
#       the same library and at the same depth. "No T at 205" means something
#       only above that floor;
#    4. the linked-site test. Residue 40 is the second position at which the two
#       evasion paralogs differ from every other member of the family (Gly vs
#       Val), and it sits 87 bp from codon 69 -- within one read. A read carrying
#       Gly-40 AND Phe-69 can only come from an evasion paralog, wherever it
#       mapped. Finding none is a stronger statement than two separate zeros.
#
#  Requirements: python3, pysam, numpy, pandas, matplotlib. pysam is needed for
#  the linked-site test (it is the only part that looks at whole reads).
#
#  Usage, from the folder that holds the BAM:
#
#      python3 avrblb2_snp205.py --bam pinsl.bam --fasta T30.fa --gtf T30.gtf \
#          --out avrblb2_snp205
#
#  Writes <out>/codon69_pileup.tsv, <out>/pooled_snp205.tsv,
#  <out>/error_floor.tsv, <out>/linked_reads.tsv and <out>/avrblb2_snp205.{pdf,png}.
#  Copy the TSVs back next to avrblb2_family.py to redraw the figure locally.
#
#  Self-test (no data needed):   python3 avrblb2_snp205.py --self-test
# ══════════════════════════════════════════════════════════════════════════════
import argparse
import collections
import os
import re
import sys

# ── The genes, and the two sites ─────────────────────────────────────────────
# The seven Avrblb2 paralogs proper (Oliva et al. 2015). All 100 aa, ungapped.
CORE_PARALOGS = ["PITG_04085", "PITG_04086", "PITG_04090", "PITG_18683",
                 "PITG_20300", "PITG_20301", "PITG_20303"]

# Expected reference residues, from the T30-4 proteins. The script checks the
# sequence it extracts against these: a mismatch means the GTF, the FASTA or the
# gene model is not the one this table was built from, and it stops rather than
# reporting counts from the wrong codon.
EXPECTED = {           # gene: (residue 40, residue 69)
    "PITG_04085": ("V", "I"), "PITG_04086": ("V", "I"), "PITG_04090": ("V", "A"),
    "PITG_18683": ("V", "I"), "PITG_20300": ("V", "A"), "PITG_20301": ("G", "F"),
    "PITG_20303": ("G", "F"),
}
EVASION_PARALOGS = ["PITG_20301", "PITG_20303"]

RES_69, RES_40 = 69, 40                      # 1-based residue numbers
NT_205 = (RES_69 - 1) * 3 + 1                # = 205, the marker of Oliva et al.
NT_118 = (RES_40 - 1) * 3 + 1                # = 118, the linked diagnostic site

# Single-copy, cleanly mapping loci: a positive control that the pileup machinery
# works and a second estimate of the error floor outside a paralog family.
# AVR3a is one locus at ~34x with mean MAPQ 60 and no MAPQ-0 reads. AVR2 is a bad
# control precisely because it has two annotated loci.
DEFAULT_CONTROLS = ["PITG_14371"]

SNP205_MEANING = {"A": "Ile-69", "G": "Ala-69 / Val-69", "T": "Phe-69 (evasion)",
                  "C": "not observed in natural populations"}

CODON_TABLE = {}
for _i, _aa in enumerate(
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    _b = "TCAG"
    CODON_TABLE[_b[_i // 16] + _b[(_i // 4) % 4] + _b[_i % 4]] = _aa

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3].upper(), "X")
                   for i in range(0, len(seq) - 2, 3))


# ── 1. gene models ───────────────────────────────────────────────────────────
def read_gtf(path, genes):
    """gene -> (chrom, strand, [(start, end), ...]) for its CDS exons, 1-based
    inclusive, sorted by genomic coordinate."""
    want, out = set(genes), {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r'(?:gene_id|Parent|gene)[ =]"?([^";]+)"?', f[8])
            if not m:
                continue
            g = m.group(1).strip()
            # transcript ids are often <gene>.1 / <gene>-T1
            g = re.sub(r"[.\-]T?\d+$", "", g)
            if g not in want:
                continue
            rec = out.setdefault(g, {"chrom": f[0], "strand": f[6], "exons": []})
            rec["exons"].append((int(f[3]), int(f[4])))
    for g, rec in out.items():
        rec["exons"].sort()
    missing = want - set(out)
    if missing:
        sys.exit(f"these genes have no CDS feature in the GTF: {sorted(missing)}")
    return out


def read_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = "".join(buf)
    return seqs


COMP_ONE = {c: d for c, d in zip("ACGTNacgtn", "TGCANtgcan")}


def gene_cds(rec, chrom_seq):
    positions = []
    for start, end in rec["exons"]:
        positions.extend(range(start, end + 1))
    bases = [chrom_seq[p - 1].upper() for p in positions]
    if rec["strand"] == "-":
        positions = positions[::-1]
        bases = [COMP_ONE[b] for b in bases[::-1]]
    return "".join(bases), positions


# ── 2. pileup ────────────────────────────────────────────────────────────────
def pileup_positions(bam, chrom, positions, strand, min_mapq=0):
    """{genomic position: Counter of bases in the GENE's orientation}. No MAPQ
    filter by default -- the multi-mapping reads are the informative ones."""
    want = set(positions)
    counts = {p: collections.Counter() for p in positions}
    lo, hi = min(positions), max(positions)
    for col in bam.pileup(chrom, lo - 1, hi, truncate=True, min_base_quality=0,
                          min_mapping_quality=min_mapq, ignore_overlaps=False,
                          max_depth=100000, stepper="nofilter"):
        p = col.reference_pos + 1
        if p not in want:
            continue
        for read in col.pileups:
            if read.is_del or read.is_refskip:
                continue
            b = read.alignment.query_sequence[read.query_position].upper()
            if strand == "-":
                b = COMP_ONE.get(b, "N")
            counts[p][b] += 1
    return counts


def codon_counts(bam, rec, positions, nt_start, strand, min_mapq=0):
    """Counts for the three bases of the codon starting at CDS position
    `nt_start` (1-based), in the gene's orientation."""
    gpos = [positions[nt_start - 1 + k] for k in range(3)]
    cnt = pileup_positions(bam, rec["chrom"], gpos, strand, min_mapq)
    return [cnt[p] for p in gpos], gpos


# ── 3. linked-site test ──────────────────────────────────────────────────────
def linked_reads(bam, rec, positions, strand):
    """Reads that cover both codon 40 and codon 69: the pair of amino acids each
    one spells. A read with Gly-40 + Phe-69 can only come from an evasion
    paralog, whatever locus it mapped to."""
    g40 = [positions[NT_118 - 1 + k] for k in range(3)]
    g69 = [positions[NT_205 - 1 + k] for k in range(3)]
    need = set(g40) | set(g69)
    lo, hi = min(need), max(need)

    def codon_from(bases_at):
        """`bases_at` maps a genomic position to its base as read off the plus
        strand. g40/g69 are already in the gene's own 5'->3' order, so a
        minus-strand gene needs its bases COMPLEMENTED, not reverse-
        complemented: reversing them again would spell the codon backwards."""
        def one(gp):
            s = "".join(bases_at[p] for p in gp)
            if strand == "-":
                s = s.translate(COMP)
            return CODON_TABLE.get(s.upper(), "X")
        return one

    pairs = collections.Counter()
    for read in bam.fetch(rec["chrom"], lo - 1, hi):
        if read.is_unmapped or read.query_sequence is None:
            continue
        at = {r + 1: q for q, r in read.get_aligned_pairs(matches_only=True)}
        if not need <= set(at):
            continue
        codon = codon_from({p: read.query_sequence[at[p]] for p in need})
        pairs[(codon(g40), codon(g69))] += 1
    return pairs


def check_linked_decoding(gene, cds, strand):
    """The same decoding applied to the reference CDS must give the residues the
    reference carries. This is what catches a strand or ordering mistake before
    any read is classified."""
    def res(nt_start):
        return CODON_TABLE.get(cds[nt_start - 1:nt_start + 2].upper(), "X")
    got = (res(NT_118), res(NT_205))
    want = EXPECTED.get(gene)
    if want and got != want:
        sys.exit(f"{gene}: the linked-site decoding gives {got[0]}-40/{got[1]}-69 "
                 f"from the reference CDS but it carries {want[0]}/{want[1]}")
    return got


# ── 4. error floor ───────────────────────────────────────────────────────────
def error_floor(bam, rec, cds, positions, strand):
    """Per-position fraction of reads carrying a base other than the reference
    one, over the whole CDS. This is the combined sequencing-error and
    cross-paralog-mismatch background: the level below which a count of zero at
    one position carries no information."""
    cnt = pileup_positions(bam, rec["chrom"], positions, strand)
    rows = []
    for i, p in enumerate(positions):
        c = cnt[p]
        total = sum(c.values())
        if total == 0:
            continue
        ref = cds[i]
        rows.append({"cds_pos": i + 1, "genomic_pos": p, "ref": ref,
                     "depth": total, "n_mismatch": total - c[ref],
                     "frac_mismatch": (total - c[ref]) / total})
    return rows


# ── 5. self-test of the coordinate arithmetic ────────────────────────────────
def self_test():
    """Checks the exon-aware, strand-aware CDS mapping on a synthetic locus with
    an intron, on both strands. No data files needed."""
    chrom = "".join("ACGT"[i % 4] for i in range(400))
    plus = {"chrom": "c", "strand": "+", "exons": [(11, 40), (61, 90)]}
    seq_p, pos_p = gene_cds(plus, chrom)
    assert len(seq_p) == 60 and len(pos_p) == 60
    assert pos_p[0] == 11 and pos_p[29] == 40 and pos_p[30] == 61
    assert seq_p == "".join(chrom[p - 1] for p in pos_p).upper()

    minus = {"chrom": "c", "strand": "-", "exons": [(11, 40), (61, 90)]}
    seq_m, pos_m = gene_cds(minus, chrom)
    assert pos_m[0] == 90 and pos_m[29] == 61 and pos_m[30] == 40
    assert seq_m == "".join(COMP_ONE[chrom[p - 1].upper()] for p in pos_m)
    # the minus-strand CDS is the reverse complement of the plus-strand one
    assert seq_m == revcomp(seq_p)

    # a minus-strand gene whose codon must come out as the gene spells it, not
    # reversed -- the case a linked-site test gets wrong when it reverses it.
    # ATT (Ile) is used on purpose: its complement spelled backwards is TTA
    # (Leu), so a reverse-complement mistake changes the answer. TTT would not
    # have caught it, since complementing and reversing TTT gives TTT again.
    gene_cds_seq = "ATG" + "GGC" * 20 + "ATT"        # ... Gly ... Ile
    plus_strand_region = revcomp(gene_cds_seq)
    chrom2 = "N" * 10 + plus_strand_region + "N" * 10
    rec2 = {"chrom": "c", "strand": "-", "exons": [(11, 10 + len(gene_cds_seq))]}
    seq2, pos2 = gene_cds(rec2, chrom2)
    assert seq2 == gene_cds_seq, (seq2, gene_cds_seq)
    assert translate(seq2)[-1] == "I" and translate(seq2)[1] == "G"
    bases_at = {p: chrom2[p - 1] for p in pos2}
    last = [pos2[len(gene_cds_seq) - 3 + k] for k in range(3)]
    spelled = "".join(bases_at[p] for p in last).translate(COMP)
    assert CODON_TABLE[spelled] == "I", spelled
    assert CODON_TABLE[revcomp(spelled)] != "I"      # the mistake is detectable
    g_cod = [pos2[3 + k] for k in range(3)]
    assert CODON_TABLE["".join(bases_at[p] for p in g_cod).translate(COMP)] == "G"

    assert NT_205 == 205 and NT_118 == 118
    assert translate("TTT") == "F" and translate("GCA") == "A"
    assert translate("ATA") == "I" and translate("GGA") == "G"
    print("self-test passed: CDS coordinates (both strands, with intron), "
          "codon arithmetic and the two site positions are consistent")


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam")
    ap.add_argument("--fasta")
    ap.add_argument("--gtf")
    ap.add_argument("--out", default="avrblb2_snp205")
    ap.add_argument("--genes", default=",".join(CORE_PARALOGS),
                    help="comma-separated; defaults to the seven Avrblb2 paralogs")
    ap.add_argument("--controls", default=",".join(DEFAULT_CONTROLS),
                    help="single-copy loci used only for the error floor")
    ap.add_argument("--mapq-report", type=int, default=30,
                    help="MAPQ level reported alongside the unfiltered counts")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        self_test()
        return
    for req in ("bam", "fasta", "gtf"):
        if not getattr(a, req):
            sys.exit(f"--{req} is required (or use --self-test)")

    import pysam
    import pandas as pd

    genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    controls = [g.strip() for g in a.controls.split(",") if g.strip()]
    os.makedirs(a.out, exist_ok=True)

    models = read_gtf(a.gtf, genes + controls)
    chroms = read_fasta(a.fasta, {m["chrom"] for m in models.values()})
    bam = pysam.AlignmentFile(a.bam, "rb")

    cds, gpos = {}, {}
    for g, rec in models.items():
        if rec["chrom"] not in chroms:
            sys.exit(f"{rec['chrom']} (for {g}) is not in {a.fasta}")
        cds[g], gpos[g] = gene_cds(rec, chroms[rec["chrom"]])
        prot = translate(cds[g])
        if g in EXPECTED:
            want40, want69 = EXPECTED[g]
            got40, got69 = prot[RES_40 - 1], prot[RES_69 - 1]
            if (got40, got69) != (want40, want69):
                sys.exit(f"{g}: expected {want40}-40/{want69}-69 from the "
                         f"reference, the gene model gives {got40}/{got69}. The "
                         f"annotation does not match the one this script was "
                         f"built for; stopping rather than counting the wrong "
                         f"codon.")
            if len(prot.rstrip("*")) != 100:
                print(f"  warning: {g} translates to {len(prot.rstrip('*'))} aa, "
                      f"not 100")
    print(f"gene models verified against the reference residues for "
          f"{len([g for g in genes if g in EXPECTED])} paralogs")

    # ── codon 69, per paralog ────────────────────────────────────────────────
    rows, pooled, pooled_q = [], collections.Counter(), collections.Counter()
    for g in genes:
        rec = models[g]
        cnts, gp = codon_counts(bam, rec, gpos[g], NT_205, rec["strand"])
        cnts_q, _ = codon_counts(bam, rec, gpos[g], NT_205, rec["strand"],
                                 min_mapq=a.mapq_report)
        for k in range(3):
            c, cq = cnts[k], cnts_q[k]
            rows.append({
                "gene": g, "codon": 69, "codon_base": k + 1,
                "cds_pos": NT_205 + k, "genomic_pos": gp[k],
                "ref_base": cds[g][NT_205 - 1 + k],
                "depth": sum(c.values()), **{f"n_{b}": c[b] for b in "ACGT"},
                f"depth_mapq{a.mapq_report}": sum(cq.values()),
                **{f"n_{b}_mapq{a.mapq_report}": cq[b] for b in "ACGT"},
            })
        pooled.update(cnts[0])          # SNP 205 = first base of the codon
        pooled_q.update(cnts_q[0])
        obs = "".join(sorted((b for b in "ACGT" if cnts[0][b]),
                             key=lambda b: -cnts[0][b]))
        print(f"  {g}: codon-69 depth {sum(cnts[0].values()):>5}  "
              f"205 bases {dict(cnts[0])}  (reference {cds[g][NT_205 - 1]})")
    pd.DataFrame(rows).to_csv(f"{a.out}/codon69_pileup.tsv", sep="\t", index=False)

    # ── pooled SNP 205 ───────────────────────────────────────────────────────
    tot = sum(pooled[b] for b in "ACGT") or 1
    pool_rows = [{"base": b, "meaning": SNP205_MEANING[b], "n_reads": pooled[b],
                  "proportion": pooled[b] / tot,
                  f"n_reads_mapq{a.mapq_report}": pooled_q[b]} for b in "ACGT"]
    pd.DataFrame(pool_rows).to_csv(f"{a.out}/pooled_snp205.tsv", sep="\t",
                                   index=False)
    genotype = "205" + "".join(b for b in "GAT" if pooled[b] > 0)
    print(f"\npooled SNP 205 over {len(genes)} paralogs, {tot} reads:")
    for b in "GATC":
        print(f"  205{b}  {pooled[b]:>6} reads  {100 * pooled[b] / tot:5.1f}%   "
              f"{SNP205_MEANING[b]}")
    print(f"  genotype: {genotype}"
          + ("   <- no Phe-69 allele in the library" if pooled["T"] == 0 else ""))

    # ── error floor ──────────────────────────────────────────────────────────
    floor = []
    for g in genes + controls:
        for r in error_floor(bam, models[g], cds[g], gpos[g], models[g]["strand"]):
            r["gene"] = g
            r["set"] = "avrblb2" if g in genes else "control"
            floor.append(r)
    fdf = pd.DataFrame(floor)
    fdf.to_csv(f"{a.out}/error_floor.tsv", sep="\t", index=False)
    fam = fdf[fdf["set"] == "avrblb2"]
    n_bases = int(fam["depth"].sum())
    n_mm = int(fam["n_mismatch"].sum())
    rate = n_mm / n_bases if n_bases else float("nan")
    print(f"\nerror floor over {len(fam)} covered coding positions of the family:")
    print(f"  {n_mm} non-reference bases in {n_bases} read bases "
          f"= {rate:.2e} per base")
    print("  (the per-position fractions are nearly all exactly zero at this "
          "depth, which is why the limit below is a confidence bound, not a "
          "percentile of that distribution)")

    # The claim is negative, so what matters is the upper bound it supports.
    # With no T observed in `tot` pooled reads, the one-sided 95% upper bound on
    # the true frequency of the Phe-69 allele is 1 - 0.05 ** (1 / tot).
    obs_T = pooled["T"]
    upper95 = 1.0 - 0.05 ** (1.0 / tot) if tot else float("nan")
    # Power: a further paralog present at the depth the present ones reach would
    # contribute that many reads to the pool.
    depths = sorted(r["depth"] for r in rows
                    if r["codon_base"] == 1 and r["depth"] > 0)
    med_depth = depths[len(depths) // 2] if depths else 0
    f_expected = med_depth / (tot + med_depth) if (tot + med_depth) else 0.0
    p_zero = (1 - f_expected) ** (tot + med_depth) if f_expected else float("nan")
    print(f"\ndetection limit for the Phe-69 allele:")
    print(f"  observed 205T: {obs_T} reads of {tot} pooled")
    print(f"  one-sided 95% upper bound on its frequency: {100 * upper95:.1f}%")
    print(f"  a further Phe-69 paralog present at the depth of the paralogs that "
          f"are present ({med_depth} reads) would make up "
          f"{100 * f_expected:.0f}% of the pool;")
    print(f"  the probability of then observing zero T reads is {p_zero:.2e}")
    pd.DataFrame([{"pooled_reads_at_205": tot, "observed_T": obs_T,
                   "upper95_frequency": upper95,
                   "error_rate_per_base": rate,
                   "n_mismatch_bases": n_mm, "n_read_bases": n_bases,
                   "median_depth_present_paralog": med_depth,
                   "expected_frequency_if_present": f_expected,
                   "p_zero_if_present": p_zero}]).to_csv(
        f"{a.out}/detection_limit.tsv", sep="\t", index=False)

    # ── linked sites ─────────────────────────────────────────────────────────
    link_rows = []
    for g in genes:
        rec = models[g]
        check_linked_decoding(g, cds[g], rec["strand"])
        for (aa40, aa69), n in linked_reads(bam, rec, gpos[g], rec["strand"]).items():
            link_rows.append({"gene": g, "res40": aa40, "res69": aa69,
                              "n_reads": n,
                              "evasion_haplotype": (aa40 == "G" and aa69 == "F")})
    ldf = pd.DataFrame(link_rows)
    ldf.to_csv(f"{a.out}/linked_reads.tsv", sep="\t", index=False)
    n_ev = int(ldf.loc[ldf["evasion_haplotype"], "n_reads"].sum()) if len(ldf) else 0
    n_sp = int(ldf["n_reads"].sum()) if len(ldf) else 0
    print(f"\nlinked-site test: {n_sp} reads span both codon 40 and codon 69")
    print(f"  carrying Gly-40 + Phe-69 (evasion haplotype): {n_ev}")

    make_figure(a.out, pd.DataFrame(rows), pd.DataFrame(pool_rows), fdf, ldf,
                genes, a.mapq_report)
    print(f"\nwrote {a.out}/*.tsv and {a.out}/avrblb2_snp205.(pdf|png)")


def make_figure(out, pile, pool, floor, linked, genes, mapq):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    BASE_COLORS = {"A": "#4daf4a", "C": "#984ea3", "G": "#377eb8", "T": "#e41a1c"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2),
                             gridspec_kw={"width_ratios": [1.5, 0.7, 1.2]})

    # A: per-paralog composition at SNP 205
    ax = axes[0]
    snp = pile[pile["codon_base"] == 1].set_index("gene").reindex(genes)
    y = np.arange(len(genes))[::-1]
    left = np.zeros(len(genes))
    for b in "ACGT":
        v = snp[f"n_{b}"].fillna(0).to_numpy(float)
        ax.barh(y, v, left=left, color=BASE_COLORS[b], height=0.7,
                edgecolor="white", linewidth=0.4)
        left += v
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=8)
    for lbl in ax.get_yticklabels():
        if lbl.get_text() in EVASION_PARALOGS:
            lbl.set_color("#c0252b")
            lbl.set_fontweight("bold")
    ax.set_xlabel("reads at SNP 205 (all MAPQ)")
    ax.set_title("Base composition at the first base of codon 69", fontsize=9)
    ax.legend(handles=[Patch(facecolor=BASE_COLORS[b],
                             label=f"205{b} — {SNP205_MEANING[b]}") for b in "ACGT"],
              fontsize=6.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)

    # B: pooled genotype, the quantity Oliva et al. measure by pyrosequencing
    ax = axes[1]
    prop = pool.set_index("base")["proportion"]
    bottom = 0.0
    for b in "GATC":
        v = float(prop.get(b, 0.0))
        if v <= 0:
            continue
        ax.bar(0, v, bottom=bottom, color=BASE_COLORS[b], width=0.55,
               edgecolor="white")
        ax.text(0, bottom + v / 2, f"205{b}\n{100 * v:.0f}%", ha="center",
                va="center", fontsize=8, color="white", fontweight="bold")
        bottom += v
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel("proportion of the pooled genome signal")
    gt = "205" + "".join(b for b in "GAT" if float(prop.get(b, 0)) > 0)
    ax.set_title(f"Pooled genotype: {gt}", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    # C: what the negative result rests on -- the frequency the data exclude
    ax = axes[2]
    n_pool = int(pool["n_reads"].sum())
    obs_T = int(pool.set_index("base")["n_reads"].get("T", 0))
    f = np.linspace(0, 0.35, 400)
    ax.plot(100 * f, (1 - f) ** n_pool, color="0.3", lw=1.5)
    upper95 = 1.0 - 0.05 ** (1.0 / n_pool) if n_pool else np.nan
    ax.axvline(100 * upper95, color="#ff7f00", ls="--", lw=1.2,
               label=f"95% upper bound = {100 * upper95:.1f}%")
    dep = pile[pile["codon_base"] == 1]
    dep = dep.loc[dep["depth"] > 0, "depth"]
    med_depth = float(dep.median()) if len(dep) else 0.0
    f_exp = med_depth / (n_pool + med_depth) if (n_pool + med_depth) else 0.0
    if f_exp:
        ax.axvline(100 * f_exp, color="#c0252b", lw=1.6,
                   label=f"a present Phe-69 paralog\nwould give {100 * f_exp:.0f}%")
    ax.axhline(0.05, color="0.75", lw=0.8, ls=":")
    ax.set_xlabel("true frequency of the Phe-69 allele in the pool (%)")
    ax.set_ylabel(f"P(observing 0 reads | {n_pool} pooled)")
    ax.set_title(f"Excluded: observed 205T = {obs_T}/{n_pool} reads", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    n_ev = int(linked.loc[linked["evasion_haplotype"], "n_reads"].sum()) if len(linked) else 0
    n_sp = int(linked["n_reads"].sum()) if len(linked) else 0
    fig.text(0.5, -0.03, f"Linked-site test: {n_sp} reads span both codon 40 and "
             f"codon 69; {n_ev} carry the Gly-40 + Phe-69 evasion haplotype.",
             ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{out}/avrblb2_snp205.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{out}/avrblb2_snp205.png", bbox_inches="tight", dpi=200)


if __name__ == "__main__":
    main()
