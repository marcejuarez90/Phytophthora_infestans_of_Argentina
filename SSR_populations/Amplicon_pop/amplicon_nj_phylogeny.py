#!/usr/bin/env python3
"""
amplicon_nj_phylogeny.py
========================
Neighbor-joining phylogeny (Newick output) from the SNPs of the 98-amplicon
panel, straight from a multi-sample VCF (BWA-MEM + bcftools).

Two analyses are run by default, each on its own sample set:

    all     every sample in the VCF
    ar_eu   the samples listed in --subset-list (20 Argentine isolates + 35
            European reference isolates, see samples_ar_eu.txt)

For each analysis (output folder <outdir>/<analysis>/):

    <analysis>_nj.nwk                 NJ tree, branch lengths + bootstrap support
                                      (percent, on internal nodes)
    <analysis>_distance_matrix.tsv    pairwise Prevosti distances
    <analysis>_sample_qc.tsv          per-sample QC: SNPs called, missing and
                                      heterozygous fractions, nearest neighbor;
                                      with --genotypes, also a label check (the
                                      genotype group each sample is closest to,
                                      on average, versus the one it is labeled with)
    <analysis>_samples_used.txt       samples that made it into the tree
    <analysis>_report.txt             what was filtered and why

Method
------
* Only biallelic SNPs (single-base REF and ALT) are used; indels and
  multi-allelic sites are skipped and counted in the report. Sites whose FILTER
  column is neither PASS nor "." are skipped too (unless --keep-nonpass).
* The genotype of a sample at a site is its ALT-allele fraction, computed from
  the GT field whatever the ploidy (diploid 0/1 -> 1/2; triploid 0/0/1 -> 1/3,
  0/1/1 -> 2/3, ...). Note that `bcftools call` is diploid: a triploid isolate
  is then reported as 0/0, 0/1 or 1/1.
* Optional per-call quality filters: --min-dp (needs FORMAT/DP) and --min-gq
  (genotype quality computed from FORMAT/PL as the gap between the two lowest
  PL values); calls below the threshold become missing.
* Distance = Prevosti's: mean over the sites called in both samples of
  |alt fraction_i - alt fraction_j|. It uses the allele dosage, so
  heterozygous calls count as intermediate rather than as missing.
* Filters (applied within each analysis, on its own sample set): sites with
  more than --max-site-missing missing calls are dropped; then samples with
  more than --max-sample-missing missing calls over the remaining sites; then
  sites that are no longer polymorphic.
* Bootstrap: SNP columns are resampled with replacement (--nboot replicates),
  an NJ tree is built for each, and each internal edge of the reference tree
  gets the percentage of replicates containing the same bipartition.
* The tree is written unrooted (trifurcating root) unless --outgroup is given.

Requirements: Python >= 3.8 and numpy. No other packages; bcftools is not
needed to run this script.

Example
-------
    python amplicon_nj_phylogeny.py --vcf totalpop_final_filter.vcf \\
        --outdir results --nboot 1000
"""
import argparse
import gzip
import os
import re
import sys
from collections import defaultdict, Counter

import numpy as np

sys.setrecursionlimit(20000)

# ── Input ────────────────────────────────────────────────────────────────────
def clean_name(name):
    """Sample name as written in the VCF header -> plain sample name
    (drops any directory and the .bam / .sorted suffixes)."""
    base = os.path.basename(name.strip())
    return re.sub(r"(\.bam|\.sorted|_sorted)+$", "", base)


def read_sample_list(path):
    names = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                names.append(clean_name(line))
    return names


def read_vcf(path, min_dp=0, min_gq=0, keep_nonpass=False):
    """Returns names, chrom, pos, ref, alt, X (sites x samples, ALT fraction,
    NaN = missing) and a dict of counters with what was skipped."""
    opener = gzip.open if path.endswith(".gz") else open
    names, chrom, pos, ref, alt, rows = None, [], [], [], [], []
    stats = Counter()
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                raw = line.rstrip("\n").split("\t")[9:]
                names = [clean_name(n) for n in raw]
                dup = [n for n, c in Counter(names).items() if c > 1]
                if dup:
                    sys.exit("ERROR: duplicated sample names in the VCF header after "
                             "cleaning: " + ", ".join(dup))
                continue
            f = line.rstrip("\n").split("\t")
            stats["sites_total"] += 1
            if not keep_nonpass and f[6] not in ("PASS", "."):
                stats["sites_skipped_filter_not_pass"] += 1
                continue
            r, a = f[3], f[4]
            if len(r) != 1 or len(a) != 1 or r not in "ACGT" or a not in "ACGT":
                stats["sites_skipped_not_biallelic_snp"] += 1
                continue
            fmt = f[8].split(":")
            if "GT" not in fmt:
                stats["sites_skipped_no_GT"] += 1
                continue
            gi = fmt.index("GT")
            di = fmt.index("DP") if "DP" in fmt else None
            qi = fmt.index("PL") if "PL" in fmt else None
            row = np.full(len(names), np.nan, dtype=np.float32)
            for k, s in enumerate(f[9:]):
                p = s.split(":")
                if gi >= len(p):
                    continue
                if min_dp > 0 and di is not None and di < len(p):
                    try:
                        if int(float(p[di])) < min_dp:
                            continue
                    except ValueError:
                        continue
                if min_gq > 0 and qi is not None and qi < len(p):
                    try:
                        pl = sorted(int(x) for x in p[qi].split(","))
                    except ValueError:
                        continue
                    if len(pl) < 2 or pl[1] - pl[0] < min_gq:
                        continue
                al = re.split(r"[/|]", p[gi])
                if not al or any(x not in ("0", "1") for x in al):
                    continue            # missing (.) or an allele other than REF/ALT
                row[k] = al.count("1") / len(al)
            chrom.append(f[0]); pos.append(int(f[1])); ref.append(r); alt.append(a)
            rows.append(row)
    if names is None:
        sys.exit("ERROR: no #CHROM header line found in " + path)
    stats["sites_used_biallelic_snps"] = len(rows)
    X = np.vstack(rows).T if rows else np.empty((len(names), 0), dtype=np.float32)
    return names, chrom, pos, ref, alt, X, stats


def read_genotype_table(path):
    """sample<TAB>genotype table (lines starting with # and a header line
    starting with 'sample' are skipped). Returns {sample: genotype}."""
    g = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] == "sample":
                continue
            g[clean_name(p[0])] = p[1]
    return g


def sample_qc(sub_names, S, D, genotypes, outpath, P):
    """Per-sample QC table. Label check (needs a genotype table): for every
    labeled sample, the mean distance to the other members of each genotype
    group; a sample is flagged when it is closer, on average, to another group
    than to its own."""
    n = len(sub_names)
    called = ~np.isnan(S)
    n_called = called.sum(1)
    miss = 1 - n_called / max(S.shape[1], 1)
    het = ((S > 0) & (S < 1)).sum(1) / np.maximum(n_called, 1)
    Dm = D.copy()
    np.fill_diagonal(Dm, np.inf)
    nn = Dm.argmin(1)
    labels = [genotypes.get(x) if genotypes else None for x in sub_names]
    groups = defaultdict(list)
    for i, l in enumerate(labels):
        if l and l != "Not assigned":
            groups[l].append(i)
    flagged = []
    with open(outpath, "w") as fh:
        fh.write("\t".join(["sample", "n_snps_called", "missing_frac", "het_frac",
                            "mean_dist_all", "nearest_neighbor", "nn_distance", "label",
                            "closest_group", "mean_dist_own_group", "mean_dist_closest_group",
                            "flag_closer_to_other_group"]) + "\n")
        for i, nm in enumerate(sub_names):
            others = [j for j in range(n) if j != i]
            lab = labels[i] or ""
            means = {h: float(D[i, [j for j in m if j != i]].mean())
                     for h, m in groups.items() if any(j != i for j in m)}
            own = means.get(lab)
            best = min(means, key=means.get) if means else ""
            flag = int(own is not None and best != lab and means[best] < own)
            if flag:
                flagged.append((nm, lab, best, own, means[best]))
            fh.write("\t".join([safe(nm), str(int(n_called[i])), format(miss[i], ".4f"),
                                format(het[i], ".4f"), format(float(D[i, others].mean()), ".4f"),
                                safe(sub_names[nn[i]]), format(float(Dm[i, nn[i]]), ".4f"), lab,
                                best if lab else "", "" if own is None else format(own, ".4f"),
                                format(means[best], ".4f") if (lab and means) else "",
                                str(flag) if lab else ""]) + "\n")
    P("wrote " + os.path.basename(outpath) + " (per-sample QC)")
    for nm, lab, best, own, bd in flagged:
        P("LABEL CHECK: %s is labeled %s but is closer, on average, to %s "
          "(mean distance %.3f vs %.3f to its own group)" % (nm, lab, best, bd, own))


# ── Distances and neighbor-joining ───────────────────────────────────────────
def prevosti_distance(X):
    """X: samples x sites (ALT fraction, NaN missing). Returns n x n matrix."""
    n = X.shape[0]
    M = ~np.isnan(X)
    Z = np.nan_to_num(X)
    D = np.zeros((n, n))
    for i in range(n):
        both = M[i] & M
        s = (np.abs(Z[i] - Z) * both).sum(1)
        c = both.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            D[i] = np.where(c > 0, s / c, np.nan)
    D = (D + D.T) / 2
    np.fill_diagonal(D, 0.0)
    n_nan = int(np.isnan(D).sum())
    if n_nan:
        D[np.isnan(D)] = np.nanmean(D)
    return D, n_nan


def neighbor_joining(D):
    """Classic NJ (Saitou & Nei). Returns (edges, n_tips, root_id); tips are
    0..n-1, internal nodes n.. ; the last join is trifurcating (unrooted)."""
    n = D.shape[0]
    M = np.zeros((2 * n, 2 * n))
    M[:n, :n] = D
    active, edges, nxt = list(range(n)), [], n
    while len(active) > 3:
        a = np.array(active)
        m = len(a)
        sub = M[np.ix_(a, a)]
        r = sub.sum(1)
        Q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(np.argmin(Q)), m)
        ia, ja, dij = a[i], a[j], sub[i, j]
        di = 0.5 * dij + (r[i] - r[j]) / (2 * (m - 2))
        dj = dij - di
        u = nxt
        nxt += 1
        for t in active:
            if t != ia and t != ja:
                M[u, t] = M[t, u] = 0.5 * (M[ia, t] + M[ja, t] - dij)
        edges += [(u, int(ia), max(di, 0.0)), (u, int(ja), max(dj, 0.0))]
        active = [t for t in active if t not in (ia, ja)] + [u]
    i, j, k = active
    dij, dik, djk = M[i, j], M[i, k], M[j, k]
    u = nxt
    edges += [(u, i, max((dij + dik - djk) / 2, 0.0)),
              (u, j, max((dij + djk - dik) / 2, 0.0)),
              (u, k, max((dik + djk - dij) / 2, 0.0))]
    return edges, n, u


def build_adjacency(edges):
    adj = defaultdict(list)
    for a, b, l in edges:
        adj[a].append((b, l))
        adj[b].append((a, l))
    return adj


def clades(adj, n, root):
    """Tip sets below every node when the tree is rooted at `root`.
    Returns (below, parent)."""
    parent, order, stack = {root: None}, [root], [root]
    while stack:
        v = stack.pop()
        for w, _ in adj[v]:
            if w != parent[v]:
                parent[w] = v
                order.append(w)
                stack.append(w)
    below = {}
    for v in reversed(order):
        if v < n:
            below[v] = frozenset([v])
        else:
            s = set()
            for w, _ in adj[v]:
                if w != parent[v]:
                    s |= below[w]
            below[v] = frozenset(s)
    return below, parent


def bipartitions(edges, n, root):
    """Non-trivial bipartitions of a tree, as a dict {child_node: key} where key
    is the side NOT containing tip 0 (so it does not depend on the rooting)."""
    adj = build_adjacency(edges)
    below, parent = clades(adj, n, root)
    full = frozenset(range(n))
    out = {}
    for v, s in below.items():
        if v < n or v == root:
            continue
        key = full - s if 0 in s else s
        if 2 <= len(key) <= n - 2:
            out[v] = key
    return out


# ── Newick ───────────────────────────────────────────────────────────────────
def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def write_newick(edges, n, root, names, support, outgroup_idx=None):
    """Newick string. `support`: {frozenset({a, b}): percent} for internal edges.
    Unrooted (trifurcating) unless outgroup_idx is given, in which case the tree
    is rooted on the edge above the smallest clade containing the outgroup."""
    adj = build_adjacency(edges)
    below, parent = clades(adj, n, root)

    def sub(v, par, label=True):
        if v < n:
            return safe(names[v])
        kids = [(w, l) for w, l in adj[v] if w != par]
        inner = ",".join(sub(w, v) + ":" + format(l, ".6g") for w, l in kids)
        lab = support.get(frozenset((v, par)), "") if (label and par is not None) else ""
        return "(" + inner + ")" + str(lab)

    warn = None
    if outgroup_idx:
        og = frozenset(outgroup_idx)
        cands = [v for v in below if v != root and og <= below[v]]
        v = min(cands, key=lambda x: len(below[x]))
        if below[v] != og:
            warn = ("outgroup is not monophyletic in the tree: rooted on the smallest "
                    "clade that contains all of it (%d tips)" % len(below[v]))
        p = parent[v]
        L = dict(adj[v])[p]
        left = sub(v, p)                  # carries the support of the root edge
        right = sub(p, v, label=False)
        return "(" + left + ":" + format(L / 2, ".6g") + "," + right + ":" + \
               format(L / 2, ".6g") + ");", warn
    return sub(root, None) + ";", warn


# ── Analysis ─────────────────────────────────────────────────────────────────
def run_analysis(label, idx, names, chrom, pos, ref, alt, X, args, outdir, excluded=(),
                 genotypes=None):
    rep = []
    P = lambda s: (rep.append(s), print("  " + s))
    os.makedirs(outdir, exist_ok=True)
    print("\n== Analysis '%s' ==" % label)
    if excluded:
        P("excluded with --exclude: " + ", ".join(excluded))

    sub_names = [names[i] for i in idx]
    S = X[idx, :].copy()                      # samples x sites
    P("samples requested: %d; sites (biallelic SNPs) before filtering: %d"
      % (len(sub_names), S.shape[1]))

    site_miss = np.isnan(S).mean(0)
    keep_sites = site_miss <= args.max_site_missing
    P("sites dropped for > %.0f%% missing calls: %d"
      % (100 * args.max_site_missing, int((~keep_sites).sum())))
    S, sites_idx = S[:, keep_sites], np.where(keep_sites)[0]

    samp_miss = np.isnan(S).mean(1)
    keep_samp = samp_miss <= args.max_sample_missing
    for k in np.where(~keep_samp)[0]:
        P("sample dropped (%.0f%% missing calls): %s" % (100 * samp_miss[k], sub_names[k]))
    S = S[keep_samp, :]
    sub_names = [n for n, k in zip(sub_names, keep_samp) if k]

    with np.errstate(all="ignore"):
        poly = np.nanmax(S, 0) > np.nanmin(S, 0)
    poly = np.where(np.isnan(poly), False, poly).astype(bool)
    P("sites dropped for being monomorphic in this sample set: %d" % int((~poly).sum()))
    S, sites_idx = S[:, poly], sites_idx[poly]
    n, ns = S.shape
    P("final: %d samples x %d SNPs (%.1f%% missing cells)"
      % (n, ns, 100 * np.isnan(S).mean() if S.size else 0))
    if n < 4 or ns < 10:
        P("SKIPPED: need at least 4 samples and 10 SNPs")
        open(os.path.join(outdir, label + "_report.txt"), "w").write("\n".join(rep) + "\n")
        return

    D, n_nan = prevosti_distance(S)
    if n_nan:
        P("WARNING: %d sample pairs share no called site; their distance was set to the mean"
          % (n_nan // 2))
    with open(os.path.join(outdir, label + "_distance_matrix.tsv"), "w") as fh:
        fh.write("\t" + "\t".join(safe(x) for x in sub_names) + "\n")
        for i, nm in enumerate(sub_names):
            fh.write(safe(nm) + "\t" + "\t".join(format(v, ".6g") for v in D[i]) + "\n")
    with open(os.path.join(outdir, label + "_samples_used.txt"), "w") as fh:
        fh.write("\n".join(sub_names) + "\n")
    sample_qc(sub_names, S, D, genotypes, os.path.join(outdir, label + "_sample_qc.tsv"), P)

    edges, nt, root = neighbor_joining(D)
    ref_bip = bipartitions(edges, nt, root)
    counts = Counter()
    rng = np.random.default_rng(args.seed)
    for b in range(args.nboot):
        cols = rng.integers(0, ns, ns)
        Db, _ = prevosti_distance(S[:, cols])
        eb, _, rb = neighbor_joining(Db)
        counts.update(set(bipartitions(eb, nt, rb).values()))
    parent = _parent_of(edges, nt, root)
    support = {}
    if args.nboot:
        for v, key in ref_bip.items():
            support[frozenset((v, parent[v]))] = int(round(100.0 * counts[key] / args.nboot))
    og_idx = []
    if args.outgroup:
        wanted = [clean_name(x) for x in args.outgroup.split(",")]
        missing = [w for w in wanted if w not in sub_names]
        if missing:
            P("WARNING: outgroup sample(s) not in this analysis, ignored: " + ", ".join(missing))
        og_idx = [sub_names.index(w) for w in wanted if w in sub_names]
    nwk, warn = write_newick(edges, nt, root, sub_names, support, og_idx or None)
    if warn:
        P("WARNING: " + warn)
    with open(os.path.join(outdir, label + "_nj.nwk"), "w") as fh:
        fh.write(nwk + "\n")
    if args.nboot:
        vals = list(support.values())
        P("bootstrap (%d replicates): %d internal edges, median support %s%%, "
          ">= 70%%: %d" % (args.nboot, len(vals), int(np.median(vals)) if vals else "NA",
                           sum(v >= 70 for v in vals)))
    P("wrote " + os.path.join(outdir, label + "_nj.nwk"))

    with open(os.path.join(outdir, label + "_report.txt"), "w") as fh:
        fh.write("\n".join(rep) + "\n")


def _parent_of(edges, n, root):
    adj = build_adjacency(edges)
    _, parent = clades(adj, n, root)
    return parent


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="multi-sample VCF (.vcf or .vcf.gz)")
    ap.add_argument("--subset-list", default=os.path.join(here, "samples_ar_eu.txt"),
                    help="sample list for the 'ar_eu' analysis (default: samples_ar_eu.txt)")
    ap.add_argument("--outdir", default="amplicon_nj_results")
    ap.add_argument("--analyses", default="all,ar_eu",
                    help="comma-separated: all, ar_eu (default: both)")
    ap.add_argument("--nboot", type=int, default=100, help="bootstrap replicates (0 = none)")
    ap.add_argument("--max-site-missing", type=float, default=0.2)
    ap.add_argument("--max-sample-missing", type=float, default=0.5)
    ap.add_argument("--min-dp", type=int, default=0,
                    help="calls with FORMAT/DP below this become missing (default: off; "
                         "only works if the VCF has per-sample DP)")
    ap.add_argument("--min-gq", type=int, default=0,
                    help="calls whose PL-derived genotype quality is below this become "
                         "missing (default: off)")
    ap.add_argument("--keep-nonpass", action="store_true",
                    help="keep sites whose FILTER is not PASS or . (default: skipped)")
    ap.add_argument("--exclude", default="",
                    help="sample(s) to leave out of every analysis: comma-separated names "
                         "or a file with one name per line")
    ap.add_argument("--genotypes", default="",
                    help="sample<TAB>genotype table (e.g. sample_genotypes.tsv): adds the label "
                         "check to <analysis>_sample_qc.tsv")
    ap.add_argument("--outgroup", default="", help="comma-separated sample name(s) to root on")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    print("Reading " + args.vcf)
    names, chrom, pos, ref, alt, X, stats = read_vcf(args.vcf, args.min_dp, args.min_gq,
                                                         args.keep_nonpass)
    print("  %d samples; %s" % (len(names), ", ".join("%s=%d" % kv for kv in stats.items())))
    os.makedirs(args.outdir, exist_ok=True)

    excl = []
    if args.exclude:
        excl = (read_sample_list(args.exclude) if os.path.isfile(args.exclude)
                else [clean_name(x) for x in args.exclude.split(",") if x.strip()])
        not_found = [e for e in excl if e not in names]
        print("Excluding: " + ", ".join(excl) +
              ("   (NOT in the VCF: " + ", ".join(not_found) + ")" if not_found else ""))
    excl = [e for e in excl if e in names]
    genotypes = read_genotype_table(args.genotypes) if args.genotypes else None

    wanted = [a.strip() for a in args.analyses.split(",") if a.strip()]
    for a in wanted:
        if a not in ("all", "ar_eu"):
            sys.exit("ERROR: unknown analysis '%s' (use all and/or ar_eu)" % a)
    if "all" in wanted:
        run_analysis("all", [i for i, n in enumerate(names) if n not in excl], names, chrom,
                     pos, ref, alt, X, args, os.path.join(args.outdir, "all"), excl, genotypes)
    if "ar_eu" in wanted:
        req = read_sample_list(args.subset_list)
        in_vcf = [r for r in req if r in names]
        present = [r for r in in_vcf if r not in excl]
        absent = [r for r in req if r not in names]      # (excluded samples are not "absent")
        print("\n'ar_eu' list: %d requested, %d found in the VCF%s" % (
            len(req), len(in_vcf),
            (", %d left out by --exclude" % (len(in_vcf) - len(present))) if len(in_vcf) != len(present) else ""))
        if absent:
            print("  NOT in the VCF (check the names written in its header): " + ", ".join(absent))
        os.makedirs(os.path.join(args.outdir, "ar_eu"), exist_ok=True)
        with open(os.path.join(args.outdir, "ar_eu", "ar_eu_samples_missing_from_vcf.txt"), "w") as fh:
            fh.write("\n".join(absent) + ("\n" if absent else ""))
        run_analysis("ar_eu", [names.index(p) for p in present], names, chrom, pos, ref, alt, X,
                     args, os.path.join(args.outdir, "ar_eu"), [e for e in excl if e in req],
                     genotypes)


if __name__ == "__main__":
    main()
