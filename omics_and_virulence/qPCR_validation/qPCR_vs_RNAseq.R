#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# RT-qPCR validation of the PiNSL-19 RNA-seq.
#
# Compares log2 fold change (infection timepoint vs. in vitro mycelium)
# measured by RT-qPCR (2^-DDCt method, log2FC = -DDCt) with the value from the
# RNA-seq of the same experiment, for four genes at 24, 48 and 72 hpi.
#
# Input : qPCR_vs_RNAseq.csv   (see README.md for the provenance of each column)
# Output: qPCR_vs_RNAseq.pdf, qPCR_vs_RNAseq.png, qPCR_vs_RNAseq_stats.txt
#
# Rscript qPCR_vs_RNAseq.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

here <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
setwd(here)

d <- read.csv("qPCR_vs_RNAseq.csv", stringsAsFactors = FALSE)

GENE_ORDER <- c("Avr2 (PITG_08943)", "EPIC2B (PITG_09173)",
                "PITG_13567", "PITG_07059")
d$gene_label <- factor(d$gene_label, levels = GENE_ORDER)
d$timepoint  <- factor(d$timepoint_hpi, levels = c(24, 48, 72),
                       labels = c("24 hpi", "48 hpi", "72 hpi"))

PAL_GENE   <- c("Avr2 (PITG_08943)"   = "#1b6ca8",
                "EPIC2B (PITG_09173)" = "#d1495b",
                "PITG_13567"          = "#2e8b57",
                "PITG_07059"          = "#e08a1e")
PAL_METHOD <- c("RNA-seq" = "#1f6f78", "RT-qPCR" = "#d96c2c")

# ------------------------------- statistics --------------------------------
x <- d$rnaseq_log2FC          # RNA-seq, the measurement being validated
y <- d$qpcr_log2FC            # RT-qPCR, the independent measurement

fit  <- lm(y ~ x)
ct   <- cor.test(x, y, method = "pearson")
sp   <- suppressWarnings(cor.test(x, y, method = "spearman"))
n    <- length(x)
r    <- unname(ct$estimate)
r2   <- summary(fit)$r.squared
b    <- unname(coef(fit)[2]); a0 <- unname(coef(fit)[1])
pval <- ct$p.value
same <- sum(sign(x) == sign(y))

# Lin's concordance correlation coefficient
ccc <- 2 * cov(x, y) * (n - 1) / n /
       (var(x) * (n - 1) / n + var(y) * (n - 1) / n + (mean(x) - mean(y))^2)

fmt_p <- function(p) if (p < 1e-4) sprintf("%.1e", p) else sprintf("%.4f", p)

out <- c(
  "RT-qPCR vs RNA-seq log2 fold change (infection vs. in vitro mycelium)",
  "----------------------------------------------------------------------",
  sprintf("genes          : %s", paste(levels(d$gene_label), collapse = ", ")),
  sprintf("timepoints     : 24, 48, 72 hpi"),
  sprintf("observations   : n = %d", n),
  "",
  sprintf("Pearson r      : %.3f  (95%% CI %.3f to %.3f, p = %s)",
          r, ct$conf.int[1], ct$conf.int[2], fmt_p(pval)),
  sprintf("R-squared      : %.3f", r2),
  sprintf("Spearman rho   : %.3f  (p = %s)", unname(sp$estimate), fmt_p(sp$p.value)),
  sprintf("Lin's CCC      : %.3f", ccc),
  sprintf("OLS qPCR ~ seq : slope = %.3f (SE %.3f), intercept = %.3f (SE %.3f)",
          b, summary(fit)$coefficients[2, 2], a0, summary(fit)$coefficients[1, 2]),
  sprintf("same direction : %d of %d observations", same, n),
  "",
  "Per-observation values",
  "----------------------"
)
tab <- d[order(d$gene_label, d$timepoint_hpi),
         c("gene", "timepoint_hpi", "qpcr_log2FC", "rnaseq_log2FC", "rnaseq_p")]
writeLines(out)
print(tab, row.names = FALSE)
con <- file("qPCR_vs_RNAseq_stats.txt", "w")
writeLines(out, con)
capture.output(print(tab, row.names = FALSE), file = con)
close(con)

# ------------------------------- panel A -----------------------------------
# Per-gene temporal profile, the two methods side by side.
long <- rbind(
  data.frame(gene_label = d$gene_label, timepoint = d$timepoint,
             method = "RNA-seq", log2FC = d$rnaseq_log2FC),
  data.frame(gene_label = d$gene_label, timepoint = d$timepoint,
             method = "RT-qPCR", log2FC = d$qpcr_log2FC)
)
long$method <- factor(long$method, levels = names(PAL_METHOD))

pA <- ggplot(long, aes(timepoint, log2FC, colour = method, group = method)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.1) +
  facet_wrap(~ gene_label, nrow = 1) +
  scale_colour_manual(values = PAL_METHOD, name = NULL) +
  scale_y_continuous(breaks = seq(-8, 12, 4)) +
  labs(x = NULL, y = expression(log[2]~"fold change vs. mycelium")) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey94", colour = "grey70"),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "top",
    legend.margin = margin(b = -4),
    legend.key.height = unit(9, "pt"),
    axis.text.x = element_text(size = 7.5)
  )

# ------------------------------- panel B -----------------------------------
lim <- range(c(x, y)); pad <- diff(lim) * 0.10
lim <- c(lim[1] - pad, lim[2] + pad)

lab <- sprintf("n = %d\nr = %.2f (p = %s)\nR² = %.2f\nslope = %.2f",
               n, r, fmt_p(pval), r2, b)

pB <- ggplot(d, aes(rnaseq_log2FC, qpcr_log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              linewidth = 0.45, colour = "grey45") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
              colour = "black", fill = "grey80", alpha = 0.45,
              linewidth = 0.6) +
  geom_point(aes(colour = gene_label, shape = timepoint),
             size = 2.6, stroke = 0.7) +
  annotate("text", x = lim[2] - diff(lim) * 0.02, y = lim[1] + diff(lim) * 0.02,
           hjust = 1, vjust = 0, label = lab, size = 2.9, lineheight = 1.15) +
  scale_colour_manual(values = PAL_GENE, name = "Gene") +
  scale_shape_manual(values = c(16, 17, 15), name = "Timepoint") +
  coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
  labs(x = expression("RNA-seq"~log[2]~"fold change"),
       y = expression("RT-qPCR"~log[2]~"fold change")) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape  = guide_legend(order = 2)) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.key.height = unit(11, "pt"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.5),
    plot.margin = margin(4, 4, 4, 4)
  )

# ------------------------------- compose -----------------------------------
fig <- plot_grid(pA, pB, ncol = 1, rel_heights = c(1, 1.55),
                 labels = c("A", "B"), label_size = 11,
                 label_fontface = "bold", align = "none")

W <- 6.7; H <- 7.1
ggsave("qPCR_vs_RNAseq.pdf", fig, width = W, height = H, device = cairo_pdf)
ggsave("qPCR_vs_RNAseq.png", fig, width = W, height = H, dpi = 400, bg = "white")
cat("\nwrote qPCR_vs_RNAseq.pdf, qPCR_vs_RNAseq.png, qPCR_vs_RNAseq_stats.txt\n")
