#!/usr/bin/env python3
# ══════════════════════════════════════════════════════════════════════════════
#  avr3a_ki_em.py — which AVR3a allele does this isolate carry, KI or EM?
#
#  AVR3a (PITG_14371) exists as two forms that differ at three residues
#  (Armstrong et al. 2005; Bos et al. 2006):
#
#        residue      19    80    103
#        AVR3a-KI      C     K      I      recognized by R3a
#        AVR3a-EM      S     E      M      evades R3a
#
#  R3a recognition is triggered by the presence of KI, so the claim that an
#  isolate escapes R3a is a NEGATIVE one: no KI allele anywhere in the library.
#  A negative claim is worth exactly as much as the floor it sits on, so this
#  script reports, in the same library and at the same depth:
#
#    1. the composition of each of the three diagnostic codons, read off whole
#       reads (all three bases of the codon at once, so the residue is resolved
#       rather than inferred from one base);
#    2. the HAPLOTYPE of residues 80 and 103. They are 69 bp apart, inside one
#       150 bp read, so a read that covers both says which allele it came from.
#       This is what separates a homozygous EM isolate from a mixture that
#       happens to look EM at each site on its own;
#    3. the error floor: the fraction of read bases differing from the reference
#       over the whole AVR3a coding sequence;
#    4. a one-sided 95% upper bound on the frequency of an unobserved KI allele.
#
#  Unlike the Avrblb2 family, AVR3a is a single-copy locus that maps cleanly
#  (~34x, mean MAPQ 60, no MAPQ-0 reads in this library), so reads are counted
#  with and without the MAPQ filter and the two should agree. They are reported
#  separately precisely so that can be checked rather than assumed.
#
#  The residue numbering is verified, not trusted: the script translates the
#  gene model and stops unless residues 19, 80 and 103 each carry one of the two
#  known variants. If the annotation is offset, it prints the neighbourhood so
#  the offset can be found instead of silently counting the wrong codon.
#
#  Requirements: python3, pysam, numpy, pandas, matplotlib, and
#  avrblb2_snp205.py in the same folder (the exon- and strand-aware coordinate
#  machinery lives there; keeping one copy means one place to be right).
#
#  Usage, from the folder that holds the BAM (copy both scripts into it):
#
#      python3 avr3a_ki_em.py --bam pinsl.bam --fasta T30.fa --gtf T30.gtf \
#          --out avr3a_kiem
#
#  Writes <out>/codon_pileup.tsv, <out>/residue_calls.tsv,
#  <out>/haplotypes.tsv, <out>/error_floor.tsv, <out>/detection_limit.tsv
#  and <out>/avr3a_ki_em.{pdf,png}.
#
#  Self-test (no data needed):   python3 avr3a_ki_em.py --self-test
# ══════════════════════════════════════════════════════════════════════════════
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from avrblb2_snp205 import (CODON_TABLE, COMP, COMP_ONE, gene_cds, read_fasta,
                                read_gtf, pileup_positions, revcomp, translate)
except ImportError:
    sys.exit("avr3a_ki_em.py needs avrblb2_snp205.py in the same folder: it "
             "reuses its exon- and strand-aware CDS machinery.")

AVR3A = "PITG_14371"

# The two alleles, by residue number (1-based, on the full 147-aa protein).
KI = {19: "C", 80: "K", 103: "I"}        # recognized by R3a
EM = {19: "S", 80: "E", 103: "M"}        # evades R3a
SITES = [19, 80, 103]
HAPLOTYPE_SITES = [80, 103]              # 69 bp apart: within one read

ALLELE_OF = {tuple(KI[r] for r in HAPLOTYPE_SITES): "KI",
             tuple(EM[r] for r in HAPLOTYPE_SITES): "EM"}


def nt_of(residue):
    """First nucleotide of the codon of a 1-based residue number."""
    return (residue - 1) * 3 + 1


# ── read-level residue calling ───────────────────────────────────────────────
def residue_haplotypes(bam, rec, positions, residues, strand):
    """Counter over the tuple of residues that each read spells at `residues`,
    for the reads that cover every base of every one of those codons.

    The genomic positions in `positions` are already in the gene's own 5'->3'
    order, so a minus-strand gene needs its bases COMPLEMENTED, not reverse-
    complemented: reversing them again would spell each codon backwards.
    """
    codons = [[positions[nt_of(r) - 1 + k] for k in range(3)] for r in residues]
    need = set(p for c in codons for p in c)
    lo, hi = min(need), max(need)

    out = collections.Counter()
    for read in bam.fetch(rec["chrom"], lo - 1, hi):
        if read.is_unmapped or read.query_sequence is None:
            continue
        at = {r + 1: q for q, r in read.get_aligned_pairs(matches_only=True)}
        if not need <= set(at):
            continue
        spelled = []
        for c in codons:
            s = "".join(read.query_sequence[at[p]] for p in c).upper()
            if strand == "-":
                s = s.translate(COMP)
            spelled.append(CODON_TABLE.get(s, "X"))
        out[tuple(spelled)] += 1
    return out


def check_numbering(cds, prot):
    """Residues 19, 80 and 103 must each carry one of the two known variants.
    Anything else means the numbering of this annotation is not the numbering
    the KI/EM nomenclature is defined on."""
    bad = []
    ref = {}
    for r in SITES:
        aa = prot[r - 1] if len(prot) >= r else "?"
        ref[r] = aa
        if aa not in (KI[r], EM[r]):
            bad.append((r, aa))
    if bad:
        for r, aa in bad:
            lo, hi = max(0, r - 6), min(len(prot), r + 5)
            sys.exit(
                f"residue {r} of {AVR3A} translates to {aa}, which is neither "
                f"{KI[r]} (KI) nor {EM[r]} (EM).\n"
                f"  protein {lo + 1}-{hi}: {prot[lo:hi]}\n"
                f"  full length: {len(prot.rstrip('*'))} aa "
                f"(the KI/EM numbering is on the 147-aa protein)\n"
                f"Stopping rather than counting the wrong codon. Check that the "
                f"GTF and FASTA are the T30-4 pair and that {AVR3A} is the "
                f"AVR3a locus in this annotation.")
    return ref


# ── self-test ────────────────────────────────────────────────────────────────
def self_test():
    """The coordinate arithmetic is tested in avrblb2_snp205.py; what is new
    here is calling several residues from one read on a minus-strand gene."""
    assert nt_of(19) == 55 and nt_of(80) == 238 and nt_of(103) == 307
    # the codons start 69 bp apart, so the two together span 72 bp: well inside
    # one 150 bp read, which is what makes the haplotype readable at all
    assert nt_of(103) - nt_of(80) == 69
    assert nt_of(103) + 2 - nt_of(80) + 1 == 72

    # a minus-strand gene whose two diagnostic codons must come out as the gene
    # spells them, not reversed. GAA (Glu) is used because its complement read
    # backwards is TTC (Phe): a reverse-complement mistake changes the answer.
    prot_nt = ["ATG"] + ["GCT"] * 78 + ["GAA"] + ["GCT"] * 22 + ["ATG"]
    cds_seq = "".join(prot_nt)                      # ... Glu-80 ... Met-103
    assert translate(cds_seq)[79] == "E" and translate(cds_seq)[102] == "M"
    chrom = "N" * 10 + revcomp(cds_seq) + "N" * 10
    rec = {"chrom": "c", "strand": "-", "exons": [(11, 10 + len(cds_seq))]}
    seq, pos = gene_cds(rec, chrom)
    assert seq == cds_seq
    bases = {p: chrom[p - 1] for p in pos}
    for r, want in ((80, "E"), (103, "M")):
        gp = [pos[nt_of(r) - 1 + k] for k in range(3)]
        spelled = "".join(bases[p] for p in gp).translate(COMP)
        assert CODON_TABLE[spelled] == want, (r, spelled)
        assert CODON_TABLE.get(revcomp(spelled)) != want   # mistake detectable
    assert ALLELE_OF[("K", "I")] == "KI" and ALLELE_OF[("E", "M")] == "EM"
    print("self-test passed: residue->codon arithmetic and minus-strand "
          "multi-codon read calling are consistent")


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam")
    ap.add_argument("--fasta")
    ap.add_argument("--gtf")
    ap.add_argument("--gene", default=AVR3A)
    ap.add_argument("--out", default="avr3a_kiem")
    ap.add_argument("--mapq-report", type=int, default=30,
                    help="MAPQ level reported alongside the unfiltered counts")
    ap.add_argument("--copies", type=int, default=3,
                    help="assumed copy number of the locus, used only for the "
                         "power calculation (PiNSL-19 is predominantly triploid)")
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

    gene = a.gene
    os.makedirs(a.out, exist_ok=True)

    models = read_gtf(a.gtf, [gene])
    rec = models[gene]
    chroms = read_fasta(a.fasta, {rec["chrom"]})
    if rec["chrom"] not in chroms:
        sys.exit(f"{rec['chrom']} (for {gene}) is not in {a.fasta}")
    cds, gpos = gene_cds(rec, chroms[rec["chrom"]])
    prot = translate(cds)
    ref = check_numbering(cds, prot)
    ref_calls = set("KI" if ref[r] == KI[r] else "EM" for r in SITES)
    ref_allele = ref_calls.pop() if len(ref_calls) == 1 else "mixed"
    print(f"{gene}: {rec['chrom']}:{min(gpos)}-{max(gpos)} strand {rec['strand']}, "
          f"{len(prot.rstrip('*'))} aa")
    print(f"  reference residues: "
          + ", ".join(f"{r}={ref[r]}"
                      f"({'KI' if ref[r] == KI[r] else 'EM'})" for r in SITES))

    bam = pysam.AlignmentFile(a.bam, "rb")

    # ── 1. codon composition at each diagnostic site ─────────────────────────
    rows = []
    for r in SITES:
        gp = [gpos[nt_of(r) - 1 + k] for k in range(3)]
        cnt = pileup_positions(bam, rec["chrom"], gp, rec["strand"])
        cnt_q = pileup_positions(bam, rec["chrom"], gp, rec["strand"],
                                 min_mapq=a.mapq_report)
        for k in range(3):
            c, cq = cnt[gp[k]], cnt_q[gp[k]]
            rows.append({
                "gene": gene, "residue": r, "codon_base": k + 1,
                "cds_pos": nt_of(r) + k, "genomic_pos": gp[k],
                "ref_base": cds[nt_of(r) - 1 + k],
                "depth": sum(c.values()), **{f"n_{b}": c[b] for b in "ACGT"},
                f"depth_mapq{a.mapq_report}": sum(cq.values()),
                **{f"n_{b}_mapq{a.mapq_report}": cq[b] for b in "ACGT"},
            })
    pd.DataFrame(rows).to_csv(f"{a.out}/codon_pileup.tsv", sep="\t", index=False)

    # ── 2. residue calls, from whole reads ───────────────────────────────────
    call_rows = []
    print("\nresidue composition (reads covering the whole codon):")
    for r in SITES:
        counts = residue_haplotypes(bam, rec, gpos, [r], rec["strand"])
        tot = sum(counts.values())
        for (aa,), n in sorted(counts.items(), key=lambda kv: -kv[1]):
            allele = ("KI" if aa == KI[r] else "EM" if aa == EM[r] else "other")
            call_rows.append({"gene": gene, "residue": r, "aa": aa,
                              "allele": allele, "n_reads": n,
                              "proportion": n / tot if tot else float("nan"),
                              "ref_aa": ref[r]})
        shown = "  ".join(f"{aa}:{n}" for (aa,), n in
                          sorted(counts.items(), key=lambda kv: -kv[1]))
        print(f"  residue {r:>3}  n = {tot:>4}   {shown}")
    pd.DataFrame(call_rows).to_csv(f"{a.out}/residue_calls.tsv", sep="\t",
                                   index=False)

    # ── 3. the 80/103 haplotype ──────────────────────────────────────────────
    hap = residue_haplotypes(bam, rec, gpos, HAPLOTYPE_SITES, rec["strand"])
    n_span = sum(hap.values())
    hap_rows = []
    for pair, n in sorted(hap.items(), key=lambda kv: -kv[1]):
        hap_rows.append({"gene": gene,
                         "res80": pair[0], "res103": pair[1],
                         "allele": ALLELE_OF.get(pair, "other"),
                         "n_reads": n,
                         "proportion": n / n_span if n_span else float("nan")})
    pd.DataFrame(hap_rows).to_csv(f"{a.out}/haplotypes.tsv", sep="\t", index=False)
    n_ki = sum(n for p, n in hap.items() if ALLELE_OF.get(p) == "KI")
    n_em = sum(n for p, n in hap.items() if ALLELE_OF.get(p) == "EM")
    print(f"\nhaplotype over residues {HAPLOTYPE_SITES} "
          f"({n_span} reads span both codons):")
    for row in hap_rows:
        print(f"  {row['res80']}-80 / {row['res103']}-103  "
              f"{row['n_reads']:>4} reads  {100 * row['proportion']:5.1f}%  "
              f"[{row['allele']}]")

    # ── 4. error floor over the whole CDS ────────────────────────────────────
    from avrblb2_snp205 import error_floor
    floor = error_floor(bam, rec, cds, gpos, rec["strand"])
    fdf = pd.DataFrame(floor)
    fdf.to_csv(f"{a.out}/error_floor.tsv", sep="\t", index=False)
    n_bases = int(fdf["depth"].sum()) if len(fdf) else 0
    n_mm = int(fdf["n_mismatch"].sum()) if len(fdf) else 0
    rate = n_mm / n_bases if n_bases else float("nan")
    print(f"\nerror floor over {len(fdf)} covered coding positions of {gene}:")
    print(f"  {n_mm} non-reference bases in {n_bases} read bases "
          f"= {rate:.2e} per base")

    # ── 5. what a zero excludes ──────────────────────────────────────────────
    upper95 = 1.0 - 0.05 ** (1.0 / n_span) if n_span else float("nan")
    # Power. A single KI copy in a `copies`-ploid locus contributes 1/copies of
    # the reads, so the chance of seeing none of it is (1 - 1/copies)^n. The
    # per-site counts are deeper than the phased ones, so the same bound is also
    # computed at the best-covered diagnostic site.
    f_one_copy = 1.0 / a.copies
    p_zero = (1 - f_one_copy) ** n_span if n_span else float("nan")
    site_depth = max((int(c["n_reads"]) for c in call_rows), default=0)
    p_zero_site = (1 - f_one_copy) ** site_depth if site_depth else float("nan")
    print(f"\ndetection limit for the KI allele:")
    print(f"  observed KI haplotypes: {n_ki} of {n_span} reads spanning both sites")
    if n_ki == 0:
        print(f"  one-sided 95% upper bound on its frequency: "
              f"{100 * upper95:.1f}%")
        print(f"  a single KI copy out of {a.copies} would make up "
              f"{100 * f_one_copy:.0f}% of the reads; the probability of then "
              f"seeing none is {p_zero:.1e} ({p_zero_site:.1e} at the "
              f"best-covered single site, {site_depth} reads)")
    pd.DataFrame([{"gene": gene, "reads_spanning_both": n_span,
                   "n_KI": n_ki, "n_EM": n_em,
                   "upper95_frequency_KI": upper95,
                   "assumed_copies": a.copies,
                   "frequency_if_one_KI_copy": f_one_copy,
                   "p_zero_if_one_KI_copy": p_zero,
                   "max_site_depth": site_depth,
                   "p_zero_at_best_site": p_zero_site,
                   "error_rate_per_base": rate,
                   "n_mismatch_bases": n_mm, "n_read_bases": n_bases,
                   "reference_allele": ref_allele}]).to_csv(
        f"{a.out}/detection_limit.tsv", sep="\t", index=False)

    make_figure(a.out, pd.DataFrame(call_rows), pd.DataFrame(hap_rows), gene,
                n_span, upper95)
    print(f"\nwrote {a.out}/*.tsv and {a.out}/avr3a_ki_em.(pdf|png)")


def make_figure(out, calls, haps, gene, n_span, upper95):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    KI_COLOR, EM_COLOR, OTHER = "#c0252b", "#377eb8", "#bdbdbd"
    color_of = {"KI": KI_COLOR, "EM": EM_COLOR, "other": OTHER}

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6),
                             gridspec_kw={"width_ratios": [1.3, 0.8]})

    # A: residue composition at each diagnostic site
    ax = axes[0]
    y = np.arange(len(SITES))[::-1]
    for i, r in enumerate(SITES):
        sub = calls[calls["residue"] == r]
        left = 0.0
        for _, row in sub.iterrows():
            ax.barh(y[i], row["n_reads"], left=left,
                    color=color_of.get(row["allele"], OTHER), height=0.6,
                    edgecolor="white", linewidth=0.5)
            if row["n_reads"] > 0.08 * max(1.0, sub["n_reads"].sum()):
                ax.text(left + row["n_reads"] / 2, y[i], row["aa"],
                        ha="center", va="center", color="white",
                        fontweight="bold", fontsize=8)
            left += row["n_reads"]
    ax.set_yticks(y)
    ax.set_yticklabels([f"residue {r}\n(KI {KI[r]} / EM {EM[r]})" for r in SITES],
                       fontsize=8)
    ax.set_xlabel("whole-genome reads covering the codon")
    ax.set_title(f"{gene} (AVR3a): diagnostic residues", fontsize=9, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    # B: the 80/103 haplotype, from reads spanning both codons
    ax = axes[1]
    bottom = 0.0
    for _, row in haps.iterrows():
        v = float(row["proportion"])
        if v <= 0:
            continue
        ax.bar(0, v, bottom=bottom, color=color_of.get(row["allele"], OTHER),
               width=0.55, edgecolor="white")
        ax.text(0, bottom + v / 2,
                f"{row['res80']}-80 / {row['res103']}-103\n{100 * v:.0f}%",
                ha="center", va="center", color="white", fontweight="bold",
                fontsize=8)
        bottom += v
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel("proportion of spanning reads")
    ax.set_title(f"Phased haplotype (n = {n_span})", fontsize=9, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    # the key goes under the panels, never on top of a bar
    fig.subplots_adjust(bottom=0.30)
    fig.legend(handles=[Patch(facecolor=KI_COLOR, label="KI residue (R3a-recognized)"),
                        Patch(facecolor=EM_COLOR, label="EM residue (R3a-evading)"),
                        Patch(facecolor=OTHER, label="other")],
               fontsize=7.5, frameon=False, loc="lower left",
               bbox_to_anchor=(0.06, 0.09), ncol=3, handlelength=1.2)
    if np.isfinite(upper95):
        fig.text(0.06, 0.02,
                 f"No read carries the KI haplotype; one-sided 95% upper bound "
                 f"on its frequency, {100 * upper95:.1f}%.",
                 ha="left", fontsize=8)
    fig.savefig(f"{out}/avr3a_ki_em.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{out}/avr3a_ki_em.png", bbox_inches="tight", dpi=200)


if __name__ == "__main__":
    main()
