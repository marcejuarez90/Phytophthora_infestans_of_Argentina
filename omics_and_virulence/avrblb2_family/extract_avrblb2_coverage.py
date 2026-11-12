#!/usr/bin/env python3
# ══════════════════════════════════════════════════════════════════════════════
#  extract_avrblb2_coverage.py — per-position WGS depth for the AVR-blb2 paralogs
#
#  Cluster-side companion of the AVR coverage data behind Figure 3D: the same
#  per-position depth profile (CDS +/- a flank, total depth and MAPQ-filtered
#  depth), for the 19 AVR-blb2 paralogs annotated in T30-4 instead of the five
#  AVR loci. Writes a CSV with the same columns as
#  omics_and_virulence/inputs/avr_gene_coverage_profile_mapq30.csv, so it is a
#  drop-in input for plot_avr_coverage_panel():
#
#      Gene_ID, Gene_name, Category, chrom, strand, genomic_pos, gene_pos_bp,
#      in_cds, depth, depth_mapq30
#
#  genomic_pos is 1-based and ascending; gene_pos_bp is relative to the CDS start
#  *in the gene's own orientation* (0 = first coding base, negative = upstream),
#  so on a minus-strand gene it descends as genomic_pos ascends. in_cds marks the
#  coding positions (exon-aware: for a multi-exon gene the introns come out False).
#
#  Needs only python3 and samtools on PATH; the BAM must be indexed (.bai).
#
#  Usage (on the cluster, from wherever the BAM and the annotation are):
#
#      python3 extract_avrblb2_coverage.py \
#          --bam pinsl.bam \
#          --annotation P_infestans_T30-4.gtf \
#          --out avrblb2_gene_coverage_profile_mapq30.csv
#
#  Then copy the CSV into omics_and_virulence/inputs/ (or next to the figure
#  script that will read it). Options: --flank (300), --mapq (30),
#  --genes / --genes-file to use another gene set, --category, --samtools,
#  --include-secondary, --max-depth.
# ══════════════════════════════════════════════════════════════════════════════
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile

# The 19 characterized AVR-blb2-family paralogs, with the strand and the genomic
# position of the first base of their residue-69 codon (from
# inputs/avrblb2_pos69_readcounts.csv). The coordinates are not used to build the
# profile — the annotation is — but they are checked against it, which catches a
# mismatched annotation or assembly before a whole run is wasted.
AVRBLB2_PARALOGS = {
    "PITG_04081": ("Supercontig_1.5", "-", 570581),
    "PITG_04085": ("Supercontig_1.5", "-", 592838),
    "PITG_04086": ("Supercontig_1.5", "+", 594061),
    "PITG_04089": ("Supercontig_1.5", "-", 705120),
    "PITG_04090": ("Supercontig_1.5", "-", 762895),
    "PITG_04097": ("Supercontig_1.5", "+", 911523),
    "PITG_04194": ("Supercontig_1.5", "-", 2020422),
    "PITG_09632": ("Supercontig_1.15", "-", 1550516),
    "PITG_18670": ("Supercontig_1.79", "+", 367076),
    "PITG_18675": ("Supercontig_1.79", "+", 407221),
    "PITG_18683": ("Supercontig_1.79", "+", 515103),
    "PITG_18685": ("Supercontig_1.79", "+", 569803),
    "PITG_20300": ("Supercontig_1.125", "+", 117394),
    "PITG_20301": ("Supercontig_1.125", "-", 158059),
    "PITG_20303": ("Supercontig_1.125", "-", 195255),
    "PITG_20857": ("Supercontig_1.216", "-", 52620),
    "PITG_22604": ("Supercontig_1.4287", "+", 1074),
    "PITG_22727": ("Supercontig_1.5", "-", 934005),
    "PITG_23193": ("Supercontig_1.216", "+", 55177),
}

TOKEN = re.compile(r"[A-Za-z0-9_.:\-]+")
TRANSCRIPT_SUFFIX = re.compile(r"(?:T\d+|\.\d+|-T\d+|\.t\d+)$")


def gene_ids_in_attributes(attributes):
    """Gene ids mentioned in a GTF/GFF3 attribute field.

    Handles both `gene_id "PITG_04090"` (GTF) and `ID=CDS:PITG_04090T0;
    Parent=transcript:PITG_04090T0` (GFF3) by taking each token, dropping a
    `gene:` / `transcript:` / `CDS:` prefix and a transcript suffix, so the
    bare gene id is what comes out.
    """
    out = set()
    for token in TOKEN.findall(attributes):
        bare = token.split(":")[-1].strip('"')
        out.add(bare)
        stripped = TRANSCRIPT_SUFFIX.sub("", bare)
        if stripped:
            out.add(stripped)
    return out


def read_annotation(path, wanted, feature="CDS"):
    """{gene_id: {"chrom", "strand", "intervals": [(start, end), ...]}} (1-based,
    inclusive), from the `feature` lines of a GTF or GFF3 file."""
    found = {}
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            hits = gene_ids_in_attributes(f[8]) & wanted
            if not hits:
                continue
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            for gid in hits:
                g = found.setdefault(gid, {"chrom": chrom, "strand": strand, "intervals": []})
                if g["chrom"] != chrom or g["strand"] != strand:
                    sys.exit(f"ERROR: {gid} appears on {g['chrom']}{g['strand']} and on "
                             f"{chrom}{strand} in {path}; cannot build one profile for it.")
                g["intervals"].append((start, end))
    for g in found.values():
        g["intervals"].sort()
    return found


def check_against_codon69(genes):
    """Warn if a paralog's known residue-69 codon does not fall in its CDS."""
    problems = []
    for gid, g in genes.items():
        ref = AVRBLB2_PARALOGS.get(gid)
        if ref is None:
            continue
        chrom, strand, pos = ref
        if g["chrom"] != chrom or g["strand"] != strand:
            problems.append(f"  {gid}: annotation says {g['chrom']}({g['strand']}), "
                            f"expected {chrom}({strand})")
        elif not any(s <= pos <= e for s, e in g["intervals"]):
            problems.append(f"  {gid}: codon-69 position {pos:,} is outside the annotated CDS "
                            f"({g['intervals'][0][0]:,}-{g['intervals'][-1][1]:,})")
    if problems:
        print("WARNING: the annotation does not match the known AVR-blb2 coordinates.\n"
              "         Check that it is the T30-4 annotation these ids come from.\n"
              + "\n".join(problems), file=sys.stderr)


def samtools_depth(samtools, bam, bed, mapq=None, max_depth=0, include_secondary=False):
    """{(chrom, pos): depth} over the BED regions, zero-coverage positions included."""
    # NOTE: in `samtools depth` it is -Q that filters on MAPPING quality and -q on
    # base quality. `samtools mpileup` uses them the other way round, so copying a
    # mpileup command here would silently filter on the wrong thing.
    cmd = [samtools, "depth", "-a", "-d", str(max_depth), "-b", bed]
    if mapq is not None:
        cmd += ["-Q", str(mapq)]
    if include_secondary:
        cmd += ["-g", "SECONDARY"]
    cmd.append(bam)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"ERROR: {' '.join(cmd)}\n{proc.stderr.strip()}")
    depth = {}
    for line in proc.stdout.splitlines():
        chrom, pos, value = line.split("\t")[:3]
        depth[(chrom, int(pos))] = int(value)
    return depth


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="indexed WGS BAM (pinsl.bam)")
    p.add_argument("--annotation", required=True, help="T30-4 annotation, GTF or GFF3 (.gz ok)")
    p.add_argument("--out", default="avrblb2_gene_coverage_profile_mapq30.csv")
    p.add_argument("--flank", type=int, default=300, help="bp of flanking sequence (default 300)")
    p.add_argument("--mapq", type=int, default=30, help="MAPQ threshold of the filtered column (default 30)")
    p.add_argument("--genes", help="comma-separated gene ids (default: the 19 AVR-blb2 paralogs)")
    p.add_argument("--genes-file", help="file with one gene id per line")
    p.add_argument("--category", default="RXLR", help="value of the Category column (default RXLR)")
    p.add_argument("--feature", default="CDS", help="annotation feature to use (default CDS)")
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--include-secondary", action="store_true",
                   help="also count secondary alignments (off by default, as in samtools depth)")
    p.add_argument("--max-depth", type=int, default=0, help="samtools depth -d (0 = no limit)")
    args = p.parse_args()

    if args.genes:
        wanted = [g.strip() for g in args.genes.split(",") if g.strip()]
    elif args.genes_file:
        with open(args.genes_file) as fh:
            wanted = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    else:
        wanted = sorted(AVRBLB2_PARALOGS)

    if shutil.which(args.samtools) is None:
        sys.exit(f"ERROR: '{args.samtools}' not found on PATH.")
    for f in (args.bam, args.annotation):
        if not os.path.exists(f):
            sys.exit(f"ERROR: {f} does not exist.")
    if not any(os.path.exists(args.bam + ext) for ext in (".bai", ".csi")) and \
            not os.path.exists(os.path.splitext(args.bam)[0] + ".bai"):
        print(f"WARNING: no index found next to {args.bam}; samtools depth -b will be slow "
              f"or fail. Run `{args.samtools} index {args.bam}` first.", file=sys.stderr)

    genes = read_annotation(args.annotation, set(wanted), feature=args.feature)
    missing = [g for g in wanted if g not in genes]
    if missing:
        print(f"WARNING: no '{args.feature}' lines for {len(missing)} gene(s): "
              f"{', '.join(missing)}", file=sys.stderr)
    if not genes:
        sys.exit("ERROR: none of the requested genes was found in the annotation.")
    check_against_codon69(genes)

    # One BED with every region, so the BAM is read once per depth pass.
    regions = {}
    for gid in wanted:
        g = genes.get(gid)
        if g is None:
            continue
        cds_start, cds_end = g["intervals"][0][0], g["intervals"][-1][1]
        regions[gid] = (g["chrom"], max(1, cds_start - args.flank), cds_end + args.flank)

    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as bed:
        for gid, (chrom, start, end) in regions.items():
            bed.write(f"{chrom}\t{start - 1}\t{end}\n")       # BED is 0-based, half-open
        bed_path = bed.name
    try:
        total = samtools_depth(args.samtools, args.bam, bed_path,
                               max_depth=args.max_depth, include_secondary=args.include_secondary)
        filtered = samtools_depth(args.samtools, args.bam, bed_path, mapq=args.mapq,
                                  max_depth=args.max_depth, include_secondary=args.include_secondary)
    finally:
        os.unlink(bed_path)

    # A MAPQ-filtered depth above the unfiltered one would mean the two passes did
    # not filter on what we think they do (the -q / -Q trap above).
    bad = [k for k, v in filtered.items() if v > total.get(k, 0)]
    if bad:
        sys.exit(f"ERROR: MAPQ-filtered depth exceeds total depth at {len(bad)} position(s), "
                 f"e.g. {bad[0]}. Check this samtools build's depth -Q flag.")

    rows, summary = [], []
    for gid in wanted:
        g = genes.get(gid)
        if g is None:
            continue
        chrom, start, end = regions[gid]
        cds_start, cds_end = g["intervals"][0][0], g["intervals"][-1][1]
        cds_depths, cds_depths_f = [], []
        for pos in range(start, end + 1):
            in_cds = any(s <= pos <= e for s, e in g["intervals"])
            gene_pos = pos - cds_start if g["strand"] != "-" else cds_end - pos
            d, df = total.get((chrom, pos), 0), filtered.get((chrom, pos), 0)
            rows.append({"Gene_ID": gid, "Gene_name": gid, "Category": args.category,
                         "chrom": chrom, "strand": g["strand"], "genomic_pos": pos,
                         "gene_pos_bp": gene_pos, "in_cds": in_cds, "depth": d,
                         "depth_mapq30": df})
            if in_cds:
                cds_depths.append(d)
                cds_depths_f.append(df)
        n = len(cds_depths) or 1
        summary.append((gid, chrom, g["strand"], len(cds_depths),
                        sum(cds_depths) / n, sum(1 for d in cds_depths if d > 0) / n,
                        sum(cds_depths_f) / n, sum(1 for d in cds_depths_f if d > 0) / n))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["Gene_ID", "Gene_name", "Category", "chrom", "strand",
                                           "genomic_pos", "gene_pos_bp", "in_cds", "depth",
                                           "depth_mapq30"])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {args.out}: {len(rows):,} positions, {len(summary)} genes "
          f"(CDS +/- {args.flank} bp, MAPQ threshold {args.mapq})\n")
    print(f"{'Gene':<12} {'contig':<18} {'str':<4} {'CDS bp':>7} "
          f"{'mean dp':>8} {'cov>=1x':>8} {'mean dpQ':>9} {'covQ>=1x':>9}")
    for gid, chrom, strand, cds_bp, mean_d, breadth, mean_f, breadth_f in summary:
        print(f"{gid:<12} {chrom:<18} {strand:<4} {cds_bp:>7,} "
              f"{mean_d:>8.1f} {breadth:>8.2f} {mean_f:>9.1f} {breadth_f:>9.2f}")
    if missing:
        print(f"\nNot in the annotation: {', '.join(missing)}")


if __name__ == "__main__":
    main()
