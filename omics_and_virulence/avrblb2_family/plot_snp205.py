#!/usr/bin/env python3
"""Figure 5B — the Rpi-blb2 evasion allele is absent from PiNSL-19.

Two panels: the base composition at SNP 205 in each Avrblb2 paralog, and the
same counts pooled. The detection limit that the negative result rests on (the
95% upper bound and the error floor) is in detection_limit.tsv and is reported
in the text rather than drawn.

Redraws the figure of avrblb2_snp205.py from the tables it wrote, so the panel
is reproducible from this repository without the BAM, exactly as
avrblb2_family.py is reproducible from the committed coverage profile.

Run from this folder:   python plot_snp205.py
Reads avrblb2_snp205/*.tsv, writes outputs/fig5B_snp205.{pdf,png}.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

IN, OUT = "avrblb2_snp205", "outputs"

CORE_PARALOGS = ["PITG_04085", "PITG_04086", "PITG_04090", "PITG_18683",
                 "PITG_20300", "PITG_20301", "PITG_20303"]
EVASION_PARALOGS = ["PITG_20301", "PITG_20303"]
BASE_COLORS = {"A": "#4daf4a", "C": "#984ea3", "G": "#377eb8", "T": "#e41a1c"}
SNP205_MEANING = {"A": "Ile-69", "G": "Ala-69 / Val-69", "T": "Phe-69 (evasion)",
                  "C": "not observed in natural populations"}
EVASION_COLOR = "#c0252b"

pile = pd.read_csv(f"{IN}/codon69_pileup.tsv", sep="\t")
pool = pd.read_csv(f"{IN}/pooled_snp205.tsv", sep="\t")
lim = pd.read_csv(f"{IN}/detection_limit.tsv", sep="\t").iloc[0]
linked = pd.read_csv(f"{IN}/linked_reads.tsv", sep="\t")

n_pool = int(pool["n_reads"].sum())
obs_T = int(pool.set_index("base")["n_reads"].get("T", 0))
upper95 = float(lim["upper95_frequency"])
f_exp = float(lim["expected_frequency_if_present"])
med_depth = int(lim["median_depth_present_paralog"])

os.makedirs(OUT, exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(9, 4.0),
                         gridspec_kw={"width_ratios": [1.5, 0.65]})

# ── A: per-paralog base composition at SNP 205 ───────────────────────────────
ax = axes[0]
snp = pile[pile["codon_base"] == 1].set_index("gene").reindex(CORE_PARALOGS)
y = np.arange(len(CORE_PARALOGS))[::-1]
left = np.zeros(len(CORE_PARALOGS))
for b in "ACGT":
    v = snp[f"n_{b}"].fillna(0).to_numpy(float)
    ax.barh(y, v, left=left, color=BASE_COLORS[b], height=0.7, edgecolor="white",
            linewidth=0.4)
    left += v
for i, g in enumerate(CORE_PARALOGS):
    if float(snp.loc[g, "depth"]) == 0:
        ax.text(0.3, y[i], "no read", va="center", fontsize=7, style="italic",
                color="0.45")
ax.set_yticks(y)
ax.set_yticklabels(CORE_PARALOGS, fontsize=8)
for lbl in ax.get_yticklabels():
    if lbl.get_text() in EVASION_PARALOGS:
        lbl.set_color(EVASION_COLOR)
        lbl.set_fontweight("bold")
ax.set_xlabel("whole-genome reads at SNP 205 (no MAPQ filter)")
ax.set_title("Base at the first position of codon 69", fontsize=9, loc="left")
ax.legend(handles=[Patch(facecolor=BASE_COLORS[b],
                         label=f"205{b} — {SNP205_MEANING[b]}") for b in "ACGT"],
          fontsize=6.5, frameon=False, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

# ── B: the pooled genotype, as Oliva et al. (2015) measure it ────────────────
ax = axes[1]
prop = pool.set_index("base")["proportion"]
bottom = 0.0
for b in "GATC":
    v = float(prop.get(b, 0.0))
    if v <= 0:
        continue
    ax.bar(0, v, bottom=bottom, color=BASE_COLORS[b], width=0.55, edgecolor="white")
    ax.text(0, bottom + v / 2, f"205{b}\n{100 * v:.0f}%", ha="center", va="center",
            fontsize=8, color="white", fontweight="bold")
    bottom += v
genotype = "205" + "".join(b for b in "GAT" if float(prop.get(b, 0)) > 0)
ax.set_xlim(-0.6, 0.6)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_ylabel("proportion of the pooled signal")
ax.set_title(f"Genotype {genotype}", fontsize=9, loc="left")
ax.spines[["top", "right"]].set_visible(False)

n_ev = int(linked.loc[linked["evasion_haplotype"], "n_reads"].sum()) if len(linked) else 0
n_sp = int(linked["n_reads"].sum()) if len(linked) else 0
fig.text(0.5, -0.04,
         f"Linked-site test: of {n_sp} reads spanning both codon 40 and codon 69, "
         f"{n_ev} carry the Gly-40 / Phe-69 haplotype diagnostic of the evasion "
         f"paralogs.", ha="center", fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/fig5B_snp205.pdf", bbox_inches="tight", dpi=300)
fig.savefig(f"{OUT}/fig5B_snp205.png", bbox_inches="tight", dpi=200)
print(f"genotype {genotype}; 205T = {obs_T}/{n_pool}; "
      f"95% upper bound {100 * upper95:.1f}%; linked evasion reads {n_ev}/{n_sp}")
print(f"wrote {OUT}/fig5B_snp205.(pdf|png)")
