#!/usr/bin/env python3
"""AVR-blb2 paralog family in isolate PiNSL-19.

Answers, for each of the 19 AVR-blb2 paralogs annotated in the T30-4 reference:

  1. which residue the reference carries at the position equivalent to 69, the
     polymorphism that controls Rpi-blb2 recognition (Oh et al. 2009: Ala/Ile/Val
     trigger the HR, Phe evades it);
  2. whether the paralog is present in the PiNSL-19 genome, from the per-base
     whole-genome read depth along its coding sequence;
  3. how much it is expressed during infection, as the gene-level mean
     log2(TPM + 1) per condition.

Point 1 is read off the reference proteins and the MAFFT alignment, never off
reads: in a family this similar, a read cannot be assigned to a paralog with any
confidence. That is also why point 2 reports the depth of uniquely mapping reads
(MAPQ >= 30) separately, and why the per-gene fraction of multi-mapped RNA-seq
reads is carried into the output table: a paralog with no whole-genome read of
its own can still collect transcripts that belong to a near-identical relative.

Run from this folder:   python avrblb2_family.py
Writes outputs/avrblb2_family.csv and outputs/avrblb2_family.{pdf,png}.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

IN = "../inputs"
COV_CSV = "inputs/avrblb2_gene_coverage_profile_mapq30.csv"
OUT = "outputs"

# position 69 is defined on PITG_20300, the paralog Oh et al. (2009) numbered
REF_GENE, REF_RESIDUE = "PITG_20300", 69
CONDITIONS = ["Myc", "24 hpi", "48 hpi", "72 hpi"]
COV_YMAX = 52.0            # common depth scale for every sparkline
COV_XLIM = (-300, 600)     # CDS +/- 300 bp, as extract_avrblb2_coverage.py writes it
CALL_COLORS = {"HR-triggering (Ala/Ile/Val)": "#4daf4a", "Evasion (Phe)": "#e41a1c",
               "Other": "#ff7f00", "No data": "#d9d9d9"}


def read_fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def classify(aa):
    if aa is None or aa == "-":
        return "No data"
    if aa in ("A", "I", "V"):
        return "HR-triggering (Ala/Ile/Val)"
    if aa == "F":
        return "Evasion (Phe)"
    return "Other"


# ── 1. residue 69 of every paralog, from the reference proteins ──────────────
prot = read_fasta(f"{IN}/paralogs_protein.fasta")
aln = read_fasta(f"{IN}/paralogs_protein_aligned.fasta")
GENES = sorted(aln)

col, seen = None, 0
for j, ch in enumerate(aln[REF_GENE]):
    if ch != "-":
        seen += 1
        if seen == REF_RESIDUE:
            col = j
            break
assert col is not None, f"residue {REF_RESIDUE} not found in {REF_GENE}"

residue, own_number = {}, {}
for g in GENES:
    ch = aln[g][col]
    residue[g] = None if ch == "-" else ch
    own_number[g] = None if ch == "-" else len(aln[g][:col + 1].replace("-", ""))
    if residue[g] is not None:
        assert prot[g][own_number[g] - 1] == residue[g], f"alignment/protein mismatch at {g}"


# ── 2. neighbor joining on the p-distance of the alignment ───────────────────
class Node:
    __slots__ = ("name", "children")

    def __init__(self, name):
        self.name, self.children = name, []


def pdistance(aln_dict, genes):
    m = np.array([list(aln_dict[g]) for g in genes])
    n = len(genes)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = np.mean(m[i] != m[j])
    return D


def neighbor_joining(labels, D):
    """ape::nj(), with the negative branch lengths of near-identical duplicates
    clamped at zero. The tree is returned as NJ builds it: an NJ tree without an
    outgroup is unrooted, so no direction should be read into its basal branches."""
    nodes = {l: Node(l) for l in labels}
    active = list(labels)
    Dm = {a: {b: D[i, j] for j, b in enumerate(labels)} for i, a in enumerate(labels)}
    nid = 0
    while len(active) > 2:
        m = len(active)
        r = {a: sum(Dm[a][b] for b in active if b != a) for a in active}
        best = None
        for ii in range(m):
            for jj in range(ii + 1, m):
                a, b = active[ii], active[jj]
                q = (m - 2) * Dm[a][b] - r[a] - r[b]
                if best is None or q < best[0]:
                    best = (q, a, b)
        _, a, b = best
        dab = Dm[a][b]
        da = max(0.5 * dab + (r[a] - r[b]) / (2 * (m - 2)), 0.0)
        u = Node(f"__internal_{nid}")
        nid += 1
        u.children = [(nodes[a], da), (nodes[b], max(dab - da, 0.0))]
        nodes[u.name] = u
        Dm[u.name] = {}
        for c in active:
            if c in (a, b):
                continue
            dc = 0.5 * (Dm[a][c] + Dm[b][c] - dab)
            Dm[u.name][c] = Dm[c][u.name] = dc
        Dm[u.name][u.name] = 0.0
        active = [c for c in active if c not in (a, b)] + [u.name]
    a, b = active
    root = Node("__root")
    root.children = [(nodes[a], max(Dm[a][b], 0.0)), (nodes[b], 0.0)]
    return root


tree = neighbor_joining(GENES, pdistance(aln, GENES))

leaf_order, depth = [], {}


def walk(node, d):
    depth[node.name] = d
    if not node.children:
        leaf_order.append(node.name)
        return
    for child, bl in node.children:
        walk(child, d + bl)


walk(tree, 0.0)


# ── 3. presence in the genome, from the whole-genome depth profile ───────────
cov_raw = pd.read_csv(COV_CSV)
COV, presence = {}, {}
for g, d in cov_raw.groupby("Gene_ID"):
    d = d.sort_values("gene_pos_bp")
    cds = d[d["in_cds"].astype(bool)]
    COV[g] = {"pos": d["gene_pos_bp"].to_numpy(),
              "depth": d["depth"].to_numpy(),
              "depth_q": d["depth_mapq30"].to_numpy(),
              "cds": (cds["gene_pos_bp"].min(), cds["gene_pos_bp"].max())}
    presence[g] = {
        "cds_len_bp": len(cds),
        "mean_depth": cds["depth"].mean(),
        "mean_depth_mapq30": cds["depth_mapq30"].mean(),
        "breadth_1x_mapq30": (cds["depth_mapq30"] >= 1).mean(),
        "present": bool(cds["depth"].max() > 0),
    }


# ── 4. gene-level expression, from the genome-wide TPM matrix ────────────────
tmat = pd.read_csv(f"{IN}/TPM_matrix.tsv", sep="\t").set_index("gene_id")
time_of = {"3": "24 hpi", "4": "48 hpi", "5": "72 hpi"}
col_cond = {c: ("Myc" if c.startswith("Ctrl") else time_of[c[1]]) for c in tmat.columns}
TPM, TPM_LOG2 = {}, {}
for g in GENES:
    row = tmat.loc[g] if g in tmat.index else None
    TPM[g] = {c: (float(np.mean([row[s] for s, cc in col_cond.items() if cc == c]))
                  if row is not None else np.nan) for c in CONDITIONS}
    TPM_LOG2[g] = {c: (float(np.log2(TPM[g][c] + 1)) if row is not None else 0.0)
                   for c in CONDITIONS}

qc = pd.read_csv(f"{IN}/effector_mapping_qc_per_gene_ext.csv").set_index("Gene_ID")


# ── 5. the table ─────────────────────────────────────────────────────────────
os.makedirs(OUT, exist_ok=True)
rows = []
for g in GENES:
    p, q = presence.get(g, {}), (qc.loc[g] if g in qc.index else {})
    rows.append({
        "Gene_ID": g,
        "residue_69_equivalent": residue[g],
        "residue_number_in_this_paralog": own_number[g],
        "Rpi_blb2_HR": {"HR-triggering (Ala/Ile/Val)": "yes", "Evasion (Phe)": "no"}
                       .get(classify(residue[g]), "unknown"),
        "wgs_mean_depth": round(p.get("mean_depth", np.nan), 2),
        "wgs_mean_depth_mapq30": round(p.get("mean_depth_mapq30", np.nan), 2),
        "wgs_breadth_1x_mapq30": round(p.get("breadth_1x_mapq30", np.nan), 3),
        "present_in_PiNSL19": p.get("present"),
        "rnaseq_reads": (q.get("n_reads_overlapping") if len(q) else None),
        "rnaseq_frac_multimapped": (q.get("frac_NHgt1") if len(q) else None),
        **{f"TPM_{c.replace(' ', '')}": round(TPM[g][c], 2) for c in CONDITIONS},
    })
tab = pd.DataFrame(rows)
tab.to_csv(f"{OUT}/avrblb2_family.csv", index=False)

absent = tab.loc[~tab["present_in_PiNSL19"].astype(bool), "Gene_ID"].tolist()
evasion = tab.loc[tab["Rpi_blb2_HR"] == "no", "Gene_ID"].tolist()
print(f"{len(GENES)} paralogs; reference residue {REF_RESIDUE} of {REF_GENE} "
      f"= alignment column {col}")
print(f"  evasion allele (Phe) in the reference: {evasion}")
print(f"  no whole-genome read over the CDS:     {absent}")
print(f"  both:                                  {sorted(set(absent) & set(evasion))}")


# ── 6. the panel ─────────────────────────────────────────────────────────────
n = len(leaf_order)
y_of = {g: n - 1 - i for i, g in enumerate(leaf_order)}
fig, ax = plt.subplots(figsize=(9.5, 6.2))

TREE_W, CW, GAP, COV_W = 1.0, 0.13, 0.06, 0.85
sx = TREE_W / max(depth[g] for g in leaf_order)


def draw(node, acc=0.0):
    if not node.children:
        return y_of[node.name]
    ys = []
    for child, bl in node.children:
        cy = draw(child, acc + bl)
        ax.plot([(acc + bl) * sx] * 2, [cy, cy], color="black", lw=0.8)
        ax.plot([acc * sx, (acc + bl) * sx], [cy, cy], color="black", lw=0.8)
        ys.append(cy)
    ax.plot([acc * sx] * 2, [min(ys), max(ys)], color="black", lw=0.8)
    return 0.5 * (min(ys) + max(ys))


draw(tree)

x_lab = TREE_W + 0.04
for g in leaf_order:
    ax.plot([depth[g] * sx + 0.01, x_lab - 0.01], [y_of[g]] * 2,
            color="0.6", lw=0.5, ls=(0, (1, 2)))
    ax.text(x_lab, y_of[g], g, ha="left", va="center", fontsize=7)

x1 = x_lab + 0.80
xc = x1 + CW + GAP
xd = xc + COV_W + GAP
h = 0.9
vmax = max(max(d.values()) for d in TPM_LOG2.values()) or 1.0
cov_x = lambda bp: xc + COV_W * (bp - COV_XLIM[0]) / (COV_XLIM[1] - COV_XLIM[0])

for g in leaf_order:
    y = y_of[g]
    ax.add_patch(Rectangle((x1, y - h / 2), CW, h, lw=0,
                           facecolor=CALL_COLORS[classify(residue[g])]))
    cov = COV.get(g)
    if cov is not None:
        base, top = y - h / 2 + 0.04, y + h / 2 - 0.02
        ax.add_patch(Rectangle((cov_x(cov["cds"][0]), base - 0.02),
                               cov_x(cov["cds"][1]) - cov_x(cov["cds"][0]),
                               top - base + 0.02, facecolor="#f2f2f2", lw=0, zorder=1))
        xs = cov_x(cov["pos"])
        ax.plot([xc, xc + COV_W], [base] * 2, color="0.75", lw=0.3, zorder=1)
        ax.fill_between(xs, base, base + (top - base) * cov["depth"] / COV_YMAX,
                        color="#9ecae1", lw=0, zorder=2)
        ax.plot(xs, base + (top - base) * cov["depth_q"] / COV_YMAX,
                color="#08519c", lw=0.5, zorder=3)
    for k, c in enumerate(CONDITIONS):
        ax.add_patch(Rectangle((xd + k * CW, y - h / 2), CW, h, lw=0.3,
                               edgecolor="white",
                               facecolor=plt.cm.viridis(TPM_LOG2[g][c] / vmax)))

top_y = n - 0.5 + 0.15
for x, lab in ([(x1 + CW / 2, "Rpi-blb2 HR"), (xc + COV_W / 2, "WGS coverage")] +
               [(xd + (k + 0.5) * CW, c) for k, c in enumerate(CONDITIONS)]):
    ax.text(x, top_y, lab, rotation=45, rotation_mode="anchor", ha="left",
            va="bottom", fontsize=7, fontweight="bold")

sb = 0.05
ax.plot([0, sb * sx], [-1.0, -1.0], color="black", lw=0.8)
ax.text(sb * sx / 2, -1.25, f"{sb}", ha="center", va="top", fontsize=7)
ax.set_xlim(-0.02, xd + 4 * CW + 0.02)
ax.set_ylim(-1.8, n + 2.4)
ax.axis("off")

# keys in a band of their own under the panel, so nothing sits on the tree
fig.subplots_adjust(bottom=0.20)
keys = [Patch(facecolor=CALL_COLORS[k], label=l) for k, l in
        [("HR-triggering (Ala/Ile/Val)", "Yes"), ("Other", "Unknown"), ("Evasion (Phe)", "No")]]
leg1 = fig.legend(handles=keys, title="Rpi-blb2 HR", loc="upper left",
                  bbox_to_anchor=(0.02, 0.17), fontsize=7, title_fontsize=7.5,
                  frameon=False, handlelength=1.0, labelspacing=0.4)
fig.add_artist(leg1)
fig.legend(handles=[Patch(facecolor="#9ecae1", label="All reads"),
                    Line2D([0], [0], color="#08519c", lw=1.0, label=r"MAPQ $\geq$ 30")],
           title=f"WGS depth (0-{COV_YMAX:.0f}x)", loc="upper left",
           bbox_to_anchor=(0.30, 0.17), fontsize=7, title_fontsize=7.5,
           frameon=False, handlelength=1.0, labelspacing=0.4)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0, vmax=vmax))
cax = fig.add_axes([0.62, 0.10, 0.24, 0.022])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label("RNA-seq expression, mean log$_2$(TPM + 1)", fontsize=7, labelpad=2)
cb.ax.tick_params(labelsize=6.5, pad=1.5)

fig.savefig(f"{OUT}/avrblb2_family.pdf", bbox_inches="tight", dpi=300)
fig.savefig(f"{OUT}/avrblb2_family.png", bbox_inches="tight", dpi=200)
print(f"wrote {OUT}/avrblb2_family.csv and {OUT}/avrblb2_family.(pdf|png)")
