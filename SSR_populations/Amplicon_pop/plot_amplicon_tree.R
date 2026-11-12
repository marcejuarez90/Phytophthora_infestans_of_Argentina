#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
#  plot_amplicon_tree.R
#
#  Draws a Newick tree from amplicon_nj_phylogeny.py with every leaf colored by
#  the genotype assigned to the sample: a colored dot on the leaf and the sample
#  name in the same (darkened) color. Argentine isolates are written in bold.
#  Support values (bootstrap, %) are shown on the internal nodes at or above a
#  threshold.
#
#  Usage (from the Amplicon_pop folder):
#
#    Rscript plot_amplicon_tree.R <tree.nwk> [genotypes.tsv] [out_prefix] [layout] [min_support] [width_in] [height_in] [legend]
#
#    tree.nwk      Newick file, e.g. results/ar_eu/ar_eu_nj.nwk
#    genotypes.tsv sample / genotype table (default: sample_genotypes.tsv next to
#                  this script). Needs the columns sample, genotype and group (others are ignored).
#                  Tree samples missing from the table are drawn as "Not assigned".
#    out_prefix    output files <out_prefix>.pdf / .png / .svg
#                  (default: <tree>_genotype next to the tree file)
#    layout        midpoint (default: rectangular tree rooted at its midpoint),
#                  unrooted (radial, no root needed) or asis (as written in the file)
#    min_support   show support >= this value (default 70)
#    width_in, height_in  figure size in inches (default: 9.5 x 0.17 per leaf + 1.8)
#    legend        right (default: outside, to the right of the tree) or bottomleft
#                  (inside the plot, in the empty area under the root; for small panels)
#
#  Requires: R with the ape package; phangorn for the midpoint root (if it is not
#  installed the tree is drawn unrooted instead).
# ══════════════════════════════════════════════════════════════════════════════
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript plot_amplicon_tree.R <tree.nwk> [genotypes.tsv] [out_prefix] [layout] [min_support]")
script_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                       error = function(e) ".")
tree_file   <- args[1]
meta_file   <- if (length(args) >= 2 && nzchar(args[2])) args[2] else file.path(script_dir, "sample_genotypes.tsv")
out_prefix  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else sub("\\.[A-Za-z]+$", "_genotype", tree_file)
layout      <- if (length(args) >= 4) args[4] else "midpoint"
min_support <- if (length(args) >= 5) as.numeric(args[5]) else 70
fig_w       <- if (length(args) >= 6 && nzchar(args[6])) as.numeric(args[6]) else NA
fig_h       <- if (length(args) >= 7 && nzchar(args[7])) as.numeric(args[7]) else NA
legend_pos  <- if (length(args) >= 8 && nzchar(args[8])) args[8] else "right"

# Genotype colors: EU_2_A1, EU_13_A2, EU_8_A1, EC1, US_1 and PE_3 use the same colors as the
# SSR figures (Figure 1); the other genotypes get their own.
geno_colors <- c(
  "EU_2_A1" = "#E63946", "EU_13_A2" = "#1f78b4", "EU_8_A1" = "#FFD700",
  "EC1" = "#FFA500", "US_1" = "#9400D3", "PE_3" = "#FF69B4",
  "PE_7" = "#00A6A6", "PE_7_Puno" = "#00695C",
  "EU_6_A1" = "#8c564b", "EU_22_A2" = "#a6cee3", "EU_5_A1" = "#33a02c",
  "EU_3_A2" = "#bcbd22", "EU_10_A2" = "#7f7f7f", "EU_17_A2" = "#cab2d6",
  "EU_36_A2" = "#b15928", "EU_37_A2" = "#393b79", "EU_12_A1" = "#999999",
  "US_6" = "#222222", "Not assigned" = "#cccccc")

tr   <- read.tree(tree_file)
meta <- read.delim(meta_file, stringsAsFactors = FALSE, comment.char = "#")
i    <- match(tr$tip.label, meta$sample)
geno  <- ifelse(is.na(i), "Not assigned", meta$genotype[i])
group <- ifelse(is.na(i), NA_character_, meta$group[i])
if (any(is.na(i)))
  message("Not in the genotype table (drawn as 'Not assigned'): ", sum(is.na(i)), " samples")

extra <- setdiff(unique(geno), names(geno_colors))     # genotypes without a preset color
if (length(extra)) geno_colors[extra] <- hcl.colors(length(extra), "Dark 3")
cols <- unname(geno_colors[geno])
darken <- function(x, f = 0.75) { m <- col2rgb(x) / 255 * f; rgb(m[1, ], m[2, ], m[3, ]) }

# Root and order. Support values are edge labels, so they have to move with the root.
if (layout == "midpoint") {
  if (requireNamespace("phangorn", quietly = TRUE)) {
    tr <- phangorn::midpoint(tr, node.labels = "support")
  } else { message("phangorn not installed: drawing unrooted"); layout <- "unrooted" }
}
if (layout != "unrooted") tr <- ladderize(tr)
# midpoint()/ladderize() may reorder tips: recompute the per-tip vectors
i     <- match(tr$tip.label, meta$sample)
geno  <- ifelse(is.na(i), "Not assigned", meta$genotype[i])
group <- ifelse(is.na(i), NA_character_, meta$group[i])
cols  <- unname(geno_colors[geno])
fonts <- ifelse(!is.na(group) & group == "Argentina", 2, 1)
sup   <- suppressWarnings(as.numeric(tr$node.label))

n <- Ntip(tr)
# Branch lengths of the NJ tree are in Prevosti distance.
scale_label <- "Prevosti distance"
draw <- function() {
  if (layout == "unrooted") {
    par(mar = c(1, 1, 1, 9))
    plot.phylo(tr, type = "unrooted", lab4ut = "axial", show.tip.label = TRUE,
               tip.color = darken(cols), font = fonts, cex = 0.7, edge.width = 1.1,
               edge.color = "grey25", no.margin = FALSE, underscore = TRUE,
               label.offset = 0.025 * max(cophenetic.phylo(tr)))   # label starts after the dot
    tiplabels(pch = 21, bg = cols, col = "grey15", cex = 1.1, adj = 0.5)
  } else {
    par(mar = c(3, 1, 1, if (legend_pos == "right") 10 else 1))
    xmax <- max(node.depth.edgelength(tr))
    plot.phylo(tr, type = "phylogram", show.tip.label = FALSE, edge.width = 1.1,
               edge.color = "grey25", x.lim = c(0, xmax * 1.5), no.margin = FALSE)
    tiplabels(pch = 21, bg = cols, col = "grey15", cex = 1.1, adj = 0.5)
    tiplabels(text = tr$tip.label, col = darken(cols), font = fonts, cex = 0.75,
              frame = "none", adj = c(-0.25, 0.5))
    show <- which(!is.na(sup) & sup >= min_support)
    if (length(show))
      nodelabels(text = sup[show], node = n + show, frame = "none", cex = 0.55,
                 col = "grey35", adj = c(1.15, -0.35))
    add.scale.bar(cex = 0.7, lwd = 1.2)
    mtext(scale_label, side = 1, line = 2, cex = 0.7, adj = 0)
  }
  tab <- table(factor(geno, levels = names(geno_colors)))
  tab <- tab[tab > 0]
  usr <- par("usr")
  leg_labels <- paste0(names(tab), "  (n = ", tab, ")")
  if (legend_pos == "bottomleft") {
    legend("bottomleft", bty = "n", cex = 0.72, title = "Genotype", title.adj = 0,
           legend = leg_labels, pch = 21, pt.cex = 1.1, pt.bg = unname(geno_colors[names(tab)]),
           col = "grey15", inset = c(0, 0.09))   # above the scale bar
    if (any(fonts == 2))
      legend("bottomright", bty = "n", cex = 0.7, legend = "bold: Argentine isolates", text.font = 2)
  } else {
    legend(x = usr[2], y = usr[4], xpd = NA, bty = "n", cex = 0.8, title = "Genotype", title.adj = 0,
           legend = leg_labels, pch = 21, pt.cex = 1.3,
           pt.bg = unname(geno_colors[names(tab)]), col = "grey15")
    if (any(fonts == 2))
      legend(x = usr[2], y = usr[3] + 0.0 * (usr[4] - usr[3]), xpd = NA, bty = "n", cex = 0.75,
             legend = "bold: Argentine isolates", text.font = 2, yjust = 0)
  }
}

w <- 9.5; h <- max(5, 0.17 * n + 1.8)
if (layout == "unrooted") { w <- 10; h <- 8 }
if (!is.na(fig_w)) w <- fig_w
if (!is.na(fig_h)) h <- fig_h
pdf(paste0(out_prefix, ".pdf"), width = w, height = h);               draw(); invisible(dev.off())
png(paste0(out_prefix, ".png"), width = w, height = h, units = "in", res = 300); draw(); invisible(dev.off())
svg(paste0(out_prefix, ".svg"), width = w, height = h);               draw(); invisible(dev.off())
cat("Wrote ", out_prefix, ".pdf / .png / .svg (", n, " leaves; layout: ", layout, ")\n", sep = "")
print(table(Genotype = geno))
