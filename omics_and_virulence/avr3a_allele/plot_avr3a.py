#!/usr/bin/env python3
"""Supplementary figure — PiNSL-19 carries AVR3a-EM and no AVR3a-KI allele.

Three panels: the whole-genome read depth along the AVR3a coding sequence with
the three polymorphic codons marked, the residue composition of each of those
codons, and the phased 80/103 haplotype of the reads that cover both.

Redraws the figure of avr3a_ki_em.py from the tables it wrote, so the panel is
reproducible from this repository without the BAM, exactly as plot_snp205.py is
for the Avrblb2 family.

Run from this folder:   python plot_avr3a.py
Reads *.tsv, writes outputs/figS5_avr3a.{pdf,png}.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

IN, OUT = ".", "outputs"

KI = {19: "C", 80: "K", 103: "I"}        # recognized by R3a
EM = {19: "S", 80: "E", 103: "M"}        # evades R3a
SITES = [19, 80, 103]
# Residue 19 (C/S) belongs to the same haplotype but sits in the signal peptide,
# and is not one of the two residues of the mature protein that determine R3a
# recognition -- which is why the alleles are named KI and EM. It is drawn, but
# in grey rather than in the colour of the allele, so the figure does not imply
# that it carries the same weight as residues 80 and 103.
GREY_SITES = {19}
GREY = "#9e9e9e"
KI_COLOR, EM_COLOR, OTHER = "#c0252b", "#377eb8", "#bdbdbd"
COLOR_OF = {"KI": KI_COLOR, "EM": EM_COLOR, "other": OTHER}

calls = pd.read_csv(f"{IN}/residue_calls.tsv", sep="\t")
haps = pd.read_csv(f"{IN}/haplotypes.tsv", sep="\t")
floor = pd.read_csv(f"{IN}/error_floor.tsv", sep="\t")
lim = pd.read_csv(f"{IN}/detection_limit.tsv", sep="\t").iloc[0]

n_span = int(lim["reads_spanning_both"])
n_ki = int(lim["n_KI"])
upper95 = float(lim["upper95_frequency_KI"])
rate = float(lim["error_rate_per_base"])
n_mm, n_read_bases = int(lim["n_mismatch_bases"]), int(lim["n_read_bases"])

os.makedirs(OUT, exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3),
                         gridspec_kw={"width_ratios": [1.35, 1.0, 0.5]})

# ── A: depth along the coding sequence, with the diagnostic codons marked ────
ax = axes[0]
ax.fill_between(floor["cds_pos"], 0, floor["depth"], color="#9ecae1", lw=0)
ax.plot(floor["cds_pos"], floor["depth"], color="#08519c", lw=0.6)
ymax = float(floor["depth"].max()) * 1.35
for r in SITES:
    nt = (r - 1) * 3 + 1
    c = GREY if r in GREY_SITES else EM_COLOR
    ax.axvspan(nt - 0.5, nt + 2.5, color=c, alpha=0.25, lw=0)
    ax.annotate(f"{r}", xy=(nt + 1, float(floor["depth"].max())),
                xytext=(0, 9), textcoords="offset points", ha="center",
                va="bottom", fontsize=7.5, color=c, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=c, lw=0.8))
ax.set_xlim(0, float(floor["cds_pos"].max()) + 1)
ax.set_ylim(0, ymax)
ax.set_xlabel("position in the AVR3a coding sequence (bp)")
ax.set_ylabel("whole-genome read depth")
ax.set_title("Coverage and the three polymorphic codons", fontsize=9,
             loc="left")
ax.spines[["top", "right"]].set_visible(False)

# ── B: residue composition at each diagnostic codon ─────────────────────────
ax = axes[1]
y = np.arange(len(SITES))[::-1]
for i, r in enumerate(SITES):
    sub = calls[calls["residue"] == r]
    total = float(sub["n_reads"].sum())
    left = 0.0
    for _, row in sub.iterrows():
        bar_color = (GREY if r in GREY_SITES
                     else COLOR_OF.get(row["allele"], OTHER))
        ax.barh(y[i], row["n_reads"], left=left, color=bar_color, height=0.6,
                edgecolor="white", linewidth=0.5)
        if row["n_reads"] > 0.12 * max(total, 1.0):
            ax.text(left + row["n_reads"] / 2, y[i], row["aa"], ha="center",
                    va="center", color="white", fontweight="bold", fontsize=8)
        left += row["n_reads"]
    ax.text(total + 0.6, y[i], f"n = {int(total)}", va="center", fontsize=7,
            color="0.35")
ax.set_yticks(y)
ax.set_yticklabels([f"residue {r}\nKI {KI[r]} / EM {EM[r]}" for r in SITES],
                   fontsize=8)
ax.set_xlim(0, float(calls.groupby("residue")["n_reads"].sum().max()) * 1.28)
ax.set_xlabel("reads covering the whole codon")
ax.set_title("Residue at each polymorphic site", fontsize=9, loc="left")
ax.spines[["top", "right"]].set_visible(False)

# ── C: the phased 80/103 haplotype ──────────────────────────────────────────
ax = axes[2]
bottom = 0.0
for _, row in haps.iterrows():
    v = float(row["proportion"])
    if v <= 0:
        continue
    ax.bar(0, v, bottom=bottom, color=COLOR_OF.get(row["allele"], OTHER),
           width=0.55, edgecolor="white")
    ax.text(0, bottom + v / 2,
            f"{row['res80']}-80\n{row['res103']}-103\n{100 * v:.0f}%",
            ha="center", va="center", color="white", fontweight="bold",
            fontsize=7.5)
    bottom += v
ax.set_xlim(-0.6, 0.6)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_ylabel("proportion of reads")
ax.set_title(f"Phased\n(n = {n_span})", fontsize=9, loc="left")
ax.spines[["top", "right"]].set_visible(False)

for letter, ax in zip("ABC", axes):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-34, 14),
                textcoords="offset points", fontsize=11, fontweight="bold",
                va="bottom", ha="left", annotation_clip=False)

fig.tight_layout()
fig.subplots_adjust(bottom=0.30, top=0.84)
# the numbers this negative result rests on go in the figure caption, not on the
# canvas: the key is all the figure itself needs
fig.legend(handles=[Patch(facecolor=KI_COLOR, label="KI residue (R3a-recognized)"),
                    Patch(facecolor=EM_COLOR, label="EM residue (R3a-evading)"),
                    Patch(facecolor=GREY, label="residue 19 (signal peptide)")],
           fontsize=7.5, frameon=False, loc="lower left",
           bbox_to_anchor=(0.055, 0.01), ncol=3, handlelength=1.2)
fig.savefig(f"{OUT}/figS5_avr3a.pdf", bbox_inches="tight", dpi=300)
fig.savefig(f"{OUT}/figS5_avr3a.png", bbox_inches="tight", dpi=200)
print(f"KI reads {n_ki}/{n_span}; 95% upper bound {100 * upper95:.1f}%; "
      f"error floor {rate:.2e} per base")
print(f"wrote {OUT}/figS5_avr3a.(pdf|png)")
