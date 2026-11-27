# ══════════════════════════════════════════════════════════════════════════════
#  P. infestans SSR — networks, PCA, PCoA, UPGMA, DAPC, AMOVA, diversity, maps
#  v4: agrega PCoA (Bruvo/cmdscale), UPGMA (ape), AMOVA (poppr.amova,
#       clone-corrected, ~lineage/location, 999 perms), ia() con permutaciones
# ══════════════════════════════════════════════════════════════════════════════

library(poppr); library(adegenet); library(ade4); library(ape)
library(ggplot2); library(dplyr); library(tidyr); library(plotly); library(knitr)
library(maps); library(scatterpie); library(ggforce); library(scales)

# ── Multi-format figure saving ─────────────────────────────────────────────────
# Every figure in this script is saved as .png (needed by make_figura1.R, which
# reads PNGs directly to compose Figure 1), .svg and .pdf (publication/vector
# formats). basepath has NO extension — each call appends the three it needs.
save_ggplot_multi <- function(plot, basepath, width, height, dpi = 600, bg = "white") {
  ggsave(paste0(basepath, ".png"), plot = plot, width = width, height = height,
         dpi = dpi, bg = bg)
  ggsave(paste0(basepath, ".svg"), plot = plot, width = width, height = height, bg = bg)
  ggsave(paste0(basepath, ".pdf"), plot = plot, width = width, height = height, bg = bg)
}

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 1 — GLOBAL PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
loci <- c("D13","PinfSSR8","PinfSSR4","Pi04","Pi70","PinfSSR6",
          "Pi63","PiG11","Pi02","PinfSSR11","PinfSSR2","Pi4B")

# Repeat length of each locus, for Bruvo's distance. Ten of the twelve are
# dinucleotide; Pi63 and Pi70 are trinucleotide, and the alleles observed here
# step accordingly (Pi70: 189/192/195; Pi63: 270/273/279).
# RL() matches the vector to the loci an object actually carries, by name:
# missingno() and poppr.amova() drop loci with more than 5% missing data, so a
# vector indexed by position would give the wrong repeat length to the loci
# that remain.
replen <- setNames(rep.int(2, length(loci)), loci)
replen[c("Pi63", "Pi70")] <- 3
RL <- function(g) replen[locNames(g)]

geno_colors <- c(
  "BR_1"      = "#a65628",
  "EC1"       = "orange",
  "EU_2_A1"   = "red",
  "EU_8_A1"   = "#FFD700",
  "Unknown"   = "#cccccc"
)

# Linetypes for genotype-group ellipses in the Fig. 1b PCoA (country-colored
# version). BR_1 and EU_8_A1 are merged into a single "old lineage" analysis
# group there (see Section 7) — three groups map cleanly onto three linetypes.
geno_linetypes <- c(
  "EU_2_A1"        = "solid",
  "EC1"            = "dashed",
  "BR_1 / EU_8_A1" = "dotted"
)

# Colors for the same merged 3-group scheme, used in the genotype-colored
# version of Fig. 1b — there both points and ellipses are colored (not
# linetyped) by this grouping, so ellipse lines match point fill exactly.
geno_colors_pcoa <- c(
  "EU_2_A1"        = "red",
  "EC1"            = "orange",
  "BR_1 / EU_8_A1" = "#FFD700"
)

# ── Unified genotype categorization + palette ─────────────────────────────────
# Single source of truth for "which named genotype does this raw Genotype
# label belong to, and what color does that named genotype get" — used
# EVERYWHERE a figure colors points/tips by genotype across the *full*
# dataset (s1): the S1-by-genotype PCoA (Section 7), the UPGMA trees
# (Section 3/run_analysis), and the map-genotype PCoA (Section 9). Defined
# here (Section 1) rather than down in Section 9 so prep_subset() (Section 2)
# and run_analysis() (Section 3) can use it too.
#
# Raw Genotype values have label variants that need folding into one category
# before coloring (otherwise e.g. "EC_1", the actual majority-label spelling
# in the CSV, silently falls outside "EC1" 's color and renders grey): EC1/
# EC_1, US7/US_7/US_7_A2, EU_2_A1/PE_7.1 (PE_7.1 shares MLGs with EU_2_A1 —
# see Results). BR_1 and EU_8_A1 (the two older reference lineages used for
# context in Section 7) are merged into a single "BR_1 / EU_8_A1" category,
# as they were in geno_colors_pcoa above. Everything else (CO_4, EU_36_A2,
# CeA, PE_7, US_8, and other low-n labels) falls into "Other" (all category
# names/labels here are the literal strings that show up in plot legends —
# kept in English throughout, no Spanish category values anywhere).
map_genotype_category <- function(genotype_vec) {
  case_when(
    genotype_vec %in% c("EU_2_A1", "PE_7.1")      ~ "EU_2_A1",
    genotype_vec %in% c("EC1", "EC_1")             ~ "EC1",
    genotype_vec == "EU_13_A2"                     ~ "EU_13_A2",
    genotype_vec %in% c("US7", "US_7", "US_7_A2")  ~ "US7",
    genotype_vec == "US_1"                         ~ "US_1",
    genotype_vec == "PE_3"                         ~ "PE_3",
    genotype_vec %in% c("BR_1", "EU_8_A1")         ~ "BR_1 / EU_8_A1",
    TRUE                                            ~ "Other"
  )
}

# Provenance of an isolate: "This work" (Source = "this work", or the legacy
# label "Flor_Lucca_check" used in some tables for the same isolates) vs.
# "Context samples" (every other source: public datasets, Euroblight reference
# panel, other published studies). Case-insensitive. Used by the provenance map
# (Section 6) and by the Source column of the Figure 1 tree (Section 7).
provenance_category <- function(src) {
  ifelse(tolower(trimws(src)) %in% c("this work", "flor_lucca_check"),
         "This work", "Context samples")
}

# Source column of the Figure 1 tree (blue / red).
provenance_colors <- c(
  "This work"       = "#2166AC",
  "Context samples" = "#D6322B"
)

map_geno_colors <- c(
  "EU_2_A1"        = "#E63946",
  "EC1"            = "orange",
  "EU_13_A2"       = "#1f78b4",
  "US7"            = "#33a02c",
  "US_1"           = "darkviolet",
  "PE_3"           = "hotpink",
  "BR_1 / EU_8_A1" = "#FFD700",
  "Other"          = "grey80"
)

map_geno_linetypes <- c(
  "EU_2_A1"        = "solid",
  "EC1"            = "dashed",
  "EU_13_A2"       = "dotted",
  "US7"            = "dotdash",
  "US_1"           = "longdash",
  "PE_3"           = "twodash",
  "BR_1 / EU_8_A1" = "solid"
)

country_colors <- c(
  "Argentina"             = "skyblue",
  "Bolivia"               = "#009E73",
  "Brazil"                = "#FFD700",
  "Chile"                 = "red",
  "Colombia"              = "purple",
  "Costa Rica"            = "#33a02c",
  "Ecuador"               = "orange",
  "EU_controls"           = "#888888",
  "Guatemala"             = "#1f78b4",
  "Honduras"              = "darkblue",
  "Panama"                = "#e7298a",
  "Peru"                  = "#8B4513",
  "Uruguay"               = "blue",
  "Multiple"              = "black"
)

period_colors <- c(
  "pre-2000"  = "darkblue",
  "2001-2010" = "#74add1",
  "2011-2020" = "#fee090",
  "2021+"     = "#d73027",
  "Unknown"   = "#cccccc"
)

map_todos_paises <- c(
  "Argentina","Bolivia","Brazil","Chile","Colombia","Costa Rica","Ecuador",
  "Guatemala","Honduras","Panama","Peru","Uruguay","Venezuela","Paraguay",
  "Guyana","Suriname","French Guiana","Nicaragua","El Salvador","Belize",
  "Trinidad and Tobago","Mexico"
)

map_coords <- tribble(
  ~Country,        ~lon,    ~lat,
  "Argentina",    -64.0,  -34.0,
  "Bolivia",      -64.9,  -16.5,
  "Brazil",       -52.0,  -14.0,
  "Chile",        -71.0,  -35.5,
  "Colombia",     -74.0,    4.0,
  "Costa Rica",   -84.0,    9.9,
  "Ecuador",      -78.0,   -1.8,
  "Guatemala",    -90.4,   15.8,
  "Honduras",     -86.6,   14.8,
  "Mexico",       -99.0,   23.0,
  "Panama",       -80.1,    8.5,
  "Peru",         -75.0,   -9.5,
  "Uruguay",      -56.2,  -32.5,
  "Venezuela",    -66.0,    8.0,
  "Paraguay",     -58.0,  -23.0
)

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — DATA LOADING AND PREPARATION
# ══════════════════════════════════════════════════════════════════════════════
# Reads the already-filtered isolate set (see 00_prepare_filtered_dataset.R):
# genotype labels normalized, Mexican isolates dropped, EC1 capped to its
# 100 most recent isolates, restricted to the lineages analyzed in the
# paper (EU_2_A1, EU_8_A1, BR_1, PE_7.1, EC1), and Euroblight reference-panel
# country codes already pooled into "EU_controls". None of that filtering
# happens here anymore — this script starts directly from the analysis-ready
# table.
df1 <- read.csv("inputs/full_dataset_SSR.csv", stringsAsFactors = FALSE)
df1$Year <- as.integer(df1$Year)

df2 <- df1[df1$Genotype %in% c("EU_2_A1","PE_7.1"), ]

prep_subset <- function(df_in) {
  # Pop_geno = true genotype/lineage identity. PE_7.1 is folded into
  # EU_2_A1 (shares MLG with European EU_2_A1 reference isolates — see
  # Results). Argentina is NOT special-cased here: its isolates keep
  # whatever genotype they actually carry (mostly EU_2_A1, but a few
  # pre-2005 controls are EU_8_A1). Country identity lives separately in
  # Pop_country, so the two dimensions stay orthogonal.
  df_in$Pop_geno <- df_in$Genotype
  df_in$Pop_geno[df_in$Pop_geno == "PE_7.1"] <- "EU_2_A1"

  # Pop_geno_map = the same unified genotype categorization used by the map
  # figures (Section 9) and now shared with the S1 PCoA-by-genotype and the
  # UPGMA trees (Section 3/7) — see map_genotype_category() in Section 1.
  # Folds label variants (EC_1->EC1, US_7/US_7_A2->US7, PE_7.1->EU_2_A1) so
  # every plot that colors by genotype uses the same categories/palette.
  df_in$Pop_geno_map <- map_genotype_category(df_in$Genotype)

  # Pop_country = true country identity. Already has European/UK reference
  # isolates pooled as "EU_controls" upstream (00_prepare_filtered_dataset.R).
  df_in$Pop_country <- df_in$Country

  df_in$Period <- cut(as.numeric(df_in$Year),
                      breaks = c(1990, 2000, 2010, 2020, 2026),
                      labels = c("pre-2000","2001-2010","2011-2020","2021+"),
                      right  = TRUE)
  df_in$Period <- as.character(df_in$Period)
  df_in$Period[is.na(df_in$Period)] <- "Unknown"

  gm <- df_in[, loci]
  sn <- make.unique(df_in$Sample)
  rownames(gm)    <- sn
  rownames(df_in) <- sn
  list(df = df_in, gm = gm, sn = sn)
}

s1 <- prep_subset(df1)
s2 <- prep_subset(df2)

cat("Subset 1:\n"); print(table(s1$df$Pop_geno))
cat("\nSubset 2:\n"); print(table(s2$df$Pop_geno))

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 3 — HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── Minimum spanning network ──────────────────────────────────────────────────
msn_plot <- function(gc, bruvo, colors_vec, seed = 42, nodelab = 5) {
  pal <- function(n) unname(colors_vec[levels(pop(gc))])
  msn <- poppr.msn(gc, bruvo, showplot = FALSE,
                   vertex.label = NA, palette = pal, nodelab = nodelab)
  set.seed(seed)
  igraph::V(msn$graph)$label <- NA
  igraph::V(msn$graph)$name  <- NA
  plot_poppr_msn(gc, msn, palette = pal, nodelab = nodelab,
                 quantiles = FALSE, vertex.label.color = "transparent")
  invisible(msn)
}

# ── PCA: compute scores and collapse by MLG ───────────────────────────────────
pca_scores <- function(gen, gc, color_vec, year_vec) {
  x   <- tab(gen, freq = TRUE, NA.method = "mean")
  pca <- dudi.pca(x, scannf = FALSE, nf = 10, center = TRUE, scale = FALSE)
  var <- round(pca$eig / sum(pca$eig) * 100, 1)
  sc  <- as.data.frame(pca$li); colnames(sc)[1:3] <- c("PC1","PC2","PC3")
  sc$Sample  <- indNames(gen)
  sc$Country <- color_vec
  sc$Year    <- year_vec
  sc$MLG     <- as.character(mll(gc))
  sc_col <- do.call(rbind, lapply(split(sc, sc$MLG), function(g) {
    lines  <- paste(sort(paste0(g$Sample," (",g$Year,")")), collapse = "<br>")
    ctries <- sort(unique(g$Country))
    data.frame(PC1 = g$PC1[1], PC2 = g$PC2[1], PC3 = g$PC3[1],
               Country = if (length(ctries) > 1) "Multiple" else ctries[1],
               n = nrow(g),
               hover = paste0("<b>",nrow(g)," isolate(s) — ",
                              paste(ctries, collapse = "/"),"</b><br>",lines),
               stringsAsFactors = FALSE)
  }))
  list(var = var, sc = sc, sc_col = sc_col)
}

# ── PCA 2D ggplot ─────────────────────────────────────────────────────────────
pca_2d <- function(sc_col, colors_vec, var, xpc = "PC1", ypc = "PC2") {
  xi <- as.integer(sub("PC","",xpc)); yi <- as.integer(sub("PC","",ypc))
  ggplot(sc_col, aes(x = .data[[xpc]], y = .data[[ypc]], fill = Country, size = n)) +
    geom_point(shape = 21, alpha = 0.85, color = "grey30", stroke = 0.3) +
    scale_fill_manual(values = colors_vec, name = "Group") +
    scale_size_continuous(range = c(2,9), name = "N isolates") +
    labs(x = paste0(xpc," (",var[xi],"%)"), y = paste0(ypc," (",var[yi],"%)")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.key = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 4)),
           size = guide_legend(override.aes = list(shape = 21, fill = "grey50")))
}

# ── PCA 3D plotly ─────────────────────────────────────────────────────────────
pca_3d <- function(sc_col, colors_vec, var) {
  sc_col$sz <- log10(sc_col$n + 1) * 6
  groups <- unique(sc_col$Country)
  p <- plot_ly(type = "scatter3d", mode = "markers")
  for (grp in groups) {
    sub <- sc_col[sc_col$Country == grp, ]
    clr <- unname(colors_vec[grp])
    if (is.na(clr)) clr <- "#cccccc"
    p <- add_trace(p, data = sub,
                   x = ~PC1, y = ~PC2, z = ~PC3,
                   name = grp,
                   text = ~hover, hoverinfo = "text",
                   marker = list(color   = clr,
                                 size    = ~sz,
                                 symbol  = "circle",
                                 opacity = 0.85))
  }
  p %>% layout(
    scene = list(
      xaxis = list(title = paste0("PC1 (",var[1],"%)")),
      yaxis = list(title = paste0("PC2 (",var[2],"%)")),
      zaxis = list(title = paste0("PC3 (",var[3],"%)"))
    ),
    legend = list(title = list(text = "<b>Group</b>"))
  )
}

# ── PCoA sobre distancias de Bruvo ────────────────────────────────────────────
# Equivalente al PCoA del texto de referencia (polysat::meandistance.matrix +
# cmdscale), pero usando bruvo.dist de poppr. El resultado es una ordenación
# en el espacio de distancias de Bruvo (no en frecuencia alélica como dudi.pca).
# group_col: columna de df_in usada para colorear ("Pop_geno" o "Pop_country")
pcoa_bruvo <- function(bruvo_dist, df_in, colors_vec,
                        group_col = "Pop_country", label = NULL) {
  mat   <- as.matrix(bruvo_dist)
  k_max <- min(5, nrow(mat) - 1)
  pcoa  <- cmdscale(mat, k = k_max, eig = TRUE)

  eig_pos <- pcoa$eig[pcoa$eig > 0]
  var_exp <- round(eig_pos / sum(eig_pos) * 100, 1)

  sc        <- as.data.frame(pcoa$points[, 1:min(2, k_max)])
  colnames(sc) <- c("PCo1","PCo2")[seq_len(ncol(sc))]
  if (ncol(sc) < 2) sc$PCo2 <- 0
  sc$Sample <- rownames(sc)
  sc$Group  <- df_in[sc$Sample, group_col]

  p <- ggplot(sc, aes(x = PCo1, y = PCo2, fill = Group)) +
    geom_point(shape = 21, size = 3, color = "grey30", stroke = 0.3, alpha = 0.85) +
    scale_fill_manual(values = colors_vec, name = "Group", na.value = "#cccccc") +
    labs(x     = sprintf("PCo1 (%.1f%%)", var_exp[1]),
         y     = sprintf("PCo2 (%.1f%%)", var_exp[2]),
         title = paste0("PCoA — Bruvo distance",
                        if (!is.null(label)) paste0("  [", label, "]"))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.key = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 4)))

  print(p)
  invisible(list(pcoa = pcoa, sc = sc, var = var_exp))
}

# ── PCoA publication figure (Fig. 1b) ─────────────────────────────────────────
# Standalone, publication-ready version of pcoa_bruvo() for Figure 1b:
# PCoA on Bruvo distances for dataset 1 (EU_2_A1 + EC1 + EU_8_A1 + BR_1 +
# PE_7.1). Points are collapsed by MLG (gc required): isolates sharing an
# identical multilocus genotype coincide at (near-)identical PCoA
# coordinates, so plotting one point per individual would just stack
# overlapping circles. Instead each point represents one MLG, sized by how
# many isolates share it. Two point styles:
#   - point_style = "point" : plain colored point per MLG (ggplot size
#     aesthetic). MLGs whose members don't all share the same
#     point_group_col value are rendered in mixed_color ("Multiple").
#   - point_style = "pie"   : geom_circle for single-group MLGs, geom_scatterpie
#     for mixed MLGs showing the real proportion of each point_group_col
#     category within that MLG (radius still scales with total n). Only
#     makes sense when point_group_col has few enough categories to read
#     as pie slices (e.g. country).
# Ellipses are computed from the full individual-level data (collapsing
# would distort the covariance/95% CI estimate), and can be styled two ways:
#   - ellipse_colors_vec supplied  → ellipses colored to match ellipse_group_col,
#     all drawn with the same fixed linetype (ellipse_linetype, dotted by
#     default) — use when points share the SAME grouping variable, e.g. both
#     by genotype, so ellipse lines visually match the fill of the points
#     they enclose)
#   - ellipse_colors_vec = NULL    → ellipses use a fixed color
#     (ellipse_color) and are distinguished only by linetype
#     (ellipse_linetypes_vec) — use when points show a DIFFERENT grouping
#     (e.g. points by country, ellipses by genotype) and color can't be
#     reused to encode both without colliding.
# coord_fixed(ratio = 1) is applied in both point styles (not just "pie") so
# a PCo1 unit and a PCo2 unit render at the same visual scale, and so two
# panels built from the same PCoA (e.g. Fig1b's country- and genotype-colored
# versions, both computed from the same bruvo_dist) come out the same size
# and proportions when placed side by side: without it the point style
# stretches to fill its canvas while the pie style keeps its real aspect
# ratio, so two panels saved at the same width and height still look
# mismatched next to each other.
# Saves directly to file since this is meant to be assembled as a standalone
# panel alongside the map and the amplicon-sequencing tree (both generated
# separately) into the final Figure 1.
pcoa_publication_plot <- function(bruvo_dist, df_in, gc,
                                   point_group_col       = "Pop_country",
                                   point_colors_vec      = country_colors,
                                   point_legend_name     = "Country",
                                   point_style            = "point",
                                   point_alpha            = 0.85,
                                   mixed_color            = "black",
                                   pie_radius_range       = NULL,
                                   ellipse_group_col      = "Pop_geno",
                                   ellipse_linetypes_vec  = geno_linetypes,
                                   ellipse_colors_vec     = NULL,
                                   ellipse_color          = "grey25",
                                   ellipse_linetype       = "dotted",
                                   min_grp_size = 3,
                                   ellipse_level = 0.95,
                                   ellipse_exclude = "Other",
                                   ellipse_unit = c("isolate", "mlg"),
                                   ellipse_linewidth = NULL,
                                   ellipse_border = NULL,
                                   axis_limits = NULL,
                                   panel_size = NULL,
                                   show_legend = TRUE,
                                   point_size_range = c(2, 9),
                                   filename     = "Fig1b_PCoA_dataset1_country.png",
                                   width = 7, height = 6, dpi = 400) {
  ellipse_unit <- match.arg(ellipse_unit)
  mat   <- as.matrix(bruvo_dist)
  k_max <- min(5, nrow(mat) - 1)
  pcoa  <- cmdscale(mat, k = k_max, eig = TRUE)

  eig_pos <- pcoa$eig[pcoa$eig > 0]
  var_exp <- round(eig_pos / sum(eig_pos) * 100, 1)

  sc           <- as.data.frame(pcoa$points[, 1:min(2, k_max)])
  colnames(sc) <- c("PCo1","PCo2")[seq_len(ncol(sc))]
  if (ncol(sc) < 2) sc$PCo2 <- 0
  sc$Sample       <- rownames(sc)
  sc$PointGroup   <- df_in[sc$Sample, point_group_col]
  # ellipse_group_col = NULL -> no ellipses at all (points only).
  use_ellipses    <- !is.null(ellipse_group_col)
  sc$EllipseGroup <- if (use_ellipses) df_in[sc$Sample, ellipse_group_col] else NA_character_
  sc$MLG          <- as.character(mll(gc))[match(sc$Sample, indNames(gc))]

  # Ellipses: computed on the full individual-level data (not collapsed).
  grp_sizes  <- table(sc$EllipseGroup)
  # Groups listed in ellipse_exclude (default "Other", a catch-all with no
  # biological meaning as a population) never get an ellipse, however large.
  valid_grps <- setdiff(names(grp_sizes[grp_sizes >= min_grp_size]), ellipse_exclude)
  sc_ellipse <- sc[sc$EllipseGroup %in% valid_grps, ]
  # ellipse_unit = "isolate": every isolate is a data point (weights each MLG
  # by its frequency). "mlg": one point per distinct MLG within each group
  # (clone-corrected — removes the weight of heavily sampled clones). The
  # min_grp_size filter above always counts isolates, so the same groups get
  # an ellipse under either unit.
  if (ellipse_unit == "mlg")
    sc_ellipse <- sc_ellipse[!duplicated(sc_ellipse[, c("MLG", "EllipseGroup")]), ]

  use_ellipse_color <- !is.null(ellipse_colors_vec)
  draw_ellipses     <- use_ellipses && nrow(sc_ellipse) > 0

  # axis_limits = list(x = c(lo, hi), y = c(lo, hi)) pins the axis ranges and
  # panel_size = c(width_in, height_in) the plotting area itself (the square
  # bounded by the axes), so different figures can share an identical panel
  # regardless of legend width — needed to place two panels side by side.
  coord_pcoa <- if (is.null(axis_limits)) coord_fixed(ratio = 1) else
    coord_fixed(ratio = 1, xlim = axis_limits$x, ylim = axis_limits$y, expand = FALSE)

  p <- ggplot()

  if (point_style == "pie") {
    # One row per MLG × PointGroup, wide format for geom_scatterpie.
    mlg_coords <- sc[!duplicated(sc$MLG), c("MLG","PCo1","PCo2")]
    mlg_counts <- as.data.frame.matrix(table(sc$MLG, sc$PointGroup))
    mlg_counts$MLG <- rownames(mlg_counts)
    pie_df <- merge(mlg_coords, mlg_counts, by = "MLG")

    country_cols <- setdiff(names(mlg_counts), "MLG")
    pie_df$total <- rowSums(pie_df[country_cols])

    span <- max(diff(range(sc$PCo1)), diff(range(sc$PCo2)))
    if (is.null(pie_radius_range)) pie_radius_range <- c(span * 0.015, span * 0.05)
    pie_df$r <- scales::rescale(sqrt(pie_df$total), to = pie_radius_range)

    is_pure   <- rowSums(pie_df[country_cols] > 0) == 1
    pie_pure  <- pie_df[is_pure, ]
    pie_mixed <- pie_df[!is_pure, ]
    pie_pure$PointGroup <- apply(pie_pure[country_cols], 1,
                                  function(row) names(row)[row > 0][1])

    if (nrow(pie_mixed) > 0)
      message(nrow(pie_mixed), " MLG(s) span more than one '", point_group_col,
              "' value — drawn as pie charts.")

    p <- p +
      geom_circle(data = pie_pure,
                  aes(x0 = PCo1, y0 = PCo2, r = r, fill = PointGroup),
                  color = "grey20", linewidth = 0.25, alpha = point_alpha)
    if (nrow(pie_mixed) > 0)
      p <- p +
        geom_scatterpie(data = pie_mixed, aes(x = PCo1, y = PCo2, r = r),
                        cols = country_cols, color = "grey20",
                        linewidth = 0.25, alpha = point_alpha)
    p <- p +
      scale_fill_manual(values = point_colors_vec, name = point_legend_name,
                         na.value = "#cccccc",
                         # legend in palette order (otherwise: order of first
                         # appearance — single-group MLGs first, pie-only ones last)
                         breaks = c(intersect(names(point_colors_vec), country_cols),
                                    setdiff(country_cols, names(point_colors_vec)))) +
      coord_pcoa

    sc_mlg <- pie_df   # returned for inspection, in place of the point-mode sc_mlg

  } else {
    # Plain point per MLG. PointGroup = shared value if the whole MLG
    # agrees, "Multiple" otherwise. n = isolate count (point size).
    mixed_mlgs <- 0L
    sc_mlg <- do.call(rbind, lapply(split(sc, sc$MLG), function(g) {
      grp_unique <- unique(g$PointGroup)
      if (length(grp_unique) > 1) mixed_mlgs <<- mixed_mlgs + 1L
      data.frame(PCo1 = g$PCo1[1], PCo2 = g$PCo2[1],
                 PointGroup = if (length(grp_unique) > 1) "Multiple" else grp_unique,
                 n = nrow(g), stringsAsFactors = FALSE)
    }))
    if (mixed_mlgs > 0)
      message(mixed_mlgs, " MLG(s) span more than one '", point_group_col,
              "' value — rendered in ", mixed_color, ".")

    point_colors_vec2 <- c(point_colors_vec, "Multiple" = mixed_color)

    p <- p +
      geom_point(data = sc_mlg,
                 aes(x = PCo1, y = PCo2, fill = PointGroup, size = n),
                 shape = 21, color = "grey20", stroke = 0.25, alpha = point_alpha) +
      scale_fill_manual(values  = point_colors_vec2, name = point_legend_name,
                         na.value = "#cccccc",
                         breaks = c(intersect(names(point_colors_vec2), unique(sc_mlg$PointGroup)),
                                    setdiff(unique(sc_mlg$PointGroup), names(point_colors_vec2)))) +
      scale_size_continuous(range = point_size_range, name = "N isolates") +
      guides(size = guide_legend(override.aes = list(shape = 21, fill = "grey50"))) +
      coord_pcoa
  }

  # Ellipses are added AFTER the points so they are drawn on top of them.
  if (draw_ellipses) {
    if (use_ellipse_color) {
      lw <- if (is.null(ellipse_linewidth)) 0.9 else ellipse_linewidth
      # Optional outline: the same ellipse drawn first in ellipse_border
      # (e.g. "black") and wider, so the colored line sits on top of it with
      # a border showing on both sides.
      if (!is.null(ellipse_border))
        p <- p + stat_ellipse(data = sc_ellipse,
                               aes(x = PCo1, y = PCo2, group = EllipseGroup),
                               geom = "path", level = ellipse_level,
                               linewidth = lw + 0.9, color = ellipse_border,
                               linetype = ellipse_linetype)
      p <- p + stat_ellipse(data = sc_ellipse,
                             aes(x = PCo1, y = PCo2, group = EllipseGroup, color = EllipseGroup),
                             geom = "path", level = ellipse_level, linewidth = lw,
                             linetype = ellipse_linetype)
    } else {
      p <- p + stat_ellipse(data = sc_ellipse,
                             aes(x = PCo1, y = PCo2, group = EllipseGroup, linetype = EllipseGroup),
                             geom = "path", level = ellipse_level, linewidth = 0.8,
                             color = ellipse_color)
    }
  }

  p <- p +
    labs(x = sprintf("PCo1 (%.1f%%)", var_exp[1]),
         y = sprintf("PCo2 (%.1f%%)", var_exp[2])) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = if (show_legend) "right" else "none",
          legend.key       = element_blank(),
          axis.title       = element_text(face = "bold")) +
    guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))

  if (!draw_ellipses) {
    # points only: nothing to add to the legend
  } else if (use_ellipse_color) {
    p <- p +
      scale_color_manual(values = ellipse_colors_vec, name = "Genotype group",
                          na.value = "#cccccc") +
      guides(color = guide_legend(override.aes = list(linewidth = 1.1)))
  } else {
    p <- p +
      scale_linetype_manual(values = ellipse_linetypes_vec, name = "Genotype group",
                             na.value = "dotdash") +
      guides(linetype = guide_legend(override.aes = list(linewidth = 0.9, color = ellipse_color)))
  }

  p_out <- p
  if (!is.null(panel_size)) {
    # legend.justification = "top": a legend taller than the panel (e.g. the
    # country + genotype-group legends) would otherwise grow the layout
    # vertically and push the panel down relative to a figure with a shorter
    # legend.
    p <- p + theme(panel.widths  = unit(panel_size[1], "in"),
                   panel.heights = unit(panel_size[2], "in"),
                   legend.justification = "top")
    # With an absolute panel size ggplot centers the whole plot on the canvas,
    # so the panel would shift sideways whenever the legend width differs
    # between figures. A flexible ("null") column/row on the right/bottom
    # absorbs the leftover space and pins the panel to the top-left corner.
    p_out <- ggplotGrob(p)
    p_out <- gtable::gtable_add_cols(p_out, unit(1, "null"))
    p_out <- gtable::gtable_add_rows(p_out, unit(1, "null"))
    # coord_fixed() sets respect = TRUE, which would make grid keep the null
    # row/column proportional and re-center the leftover space; the panel is
    # already absolute, so turn it off.
    p_out$respect <- FALSE
  }

  print(p)
  save_ggplot_multi(p_out, tools::file_path_sans_ext(filename),
                     width = width, height = height, dpi = dpi, bg = "white")

  invisible(list(pcoa = pcoa, sc = sc, sc_mlg = sc_mlg, var = var_exp, plot = p))
}

# ── UPGMA dendrogram ──────────────────────────────────────────────────────────
# Genera el dendrograma UPGMA a partir de la matriz de distancias de Bruvo.
# group_col : columna de df_in para colorear los tip labels.
# cex_tip   : tamaño del texto de las etiquetas (reducir si hay muchas muestras).
upgma_plot <- function(bruvo_dist, df_in, colors_vec,
                        group_col = "Pop_country",
                        cex_tip   = 0.45,
                        label     = NULL,
                        filename  = NULL,
                        width     = 8, height = 10) {
  tree <- upgma(bruvo_dist)

  tip_groups <- df_in[tree$tip.label, group_col]
  tip_cols   <- unname(colors_vec[tip_groups])
  tip_cols[is.na(tip_cols)] <- "#cccccc"

  draw_tree <- function() {
    plot(tree,
         type         = "phylogram",
         tip.color    = tip_cols,
         cex          = cex_tip,
         label.offset = 0.001,
         no.margin    = TRUE,
         main         = paste0("UPGMA — Bruvo distance",
                               if (!is.null(label)) paste0("  [", label, "]")))

    present      <- unique(na.omit(tip_groups))
    present_cols <- colors_vec[present[present %in% names(colors_vec)]]
    if (length(present_cols) > 0)
      legend("topright", legend = names(present_cols),
             fill = present_cols, cex = 0.6, bty = "n", title = "Group")
  }

  draw_tree()  # on-screen, same behavior as before

  # base-R phylo plots aren't ggplot objects, so they can't go through
  # save_ggplot_multi() — redraw the same tree into svg/pdf devices instead.
  # No .png here (this function was never producing one, and nothing else in
  # the pipeline reads a UPGMA tree back in as an image).
  if (!is.null(filename)) {
    basepath <- tools::file_path_sans_ext(filename)
    svg(paste0(basepath, ".svg"), width = width, height = height); draw_tree(); dev.off()
    pdf(paste0(basepath, ".pdf"), width = width, height = height); draw_tree(); dev.off()
  }

  invisible(tree)
}

# ── AMOVA ─────────────────────────────────────────────────────────────────────
# Particiona la varianza genética según la fórmula indicada.
# Se aplica clone correction antes del cálculo, reteniendo un representante
# por MLG por combinación de estratos (equivalente al paper de referencia).
#
# formula_str : fórmula como string, e.g. "~Pop_geno/Pop_country"
# strata_cols : columnas de df_in incluidas en la fórmula,
#               en el mismo orden que aparecen en formula_str
#               e.g. c("Pop_geno","Pop_country")
# n_perm      : permutaciones Monte Carlo para significancia
#
# Interpretación de los resultados:
#   Between X            : varianza entre niveles del factor externo (e.g. linajes)
#   Between Y within X   : varianza entre localidades dentro de cada linaje
#   Within samples       : varianza dentro de los individuos
#   p-value < 0.05 → diferenciación significativa a ese nivel
amova_analysis <- function(gm, df_in, sn, formula_str, strata_cols,
                            prefix, n_perm = 999, seed = 42, outdir = ".") {
  tryCatch({
    cat(paste0("\n── AMOVA: ", formula_str, "  [", prefix, "] ──\n"))

    gc2      <- as.genclone(df2genind(gm, sep = "/", ploidy = 3, ind.names = sn))
    strata_df <- df_in[sn, strata_cols, drop = FALSE]
    strata(gc2) <- strata_df

    gc_cc    <- clonecorrect(gc2, strata = as.formula(formula_str))
    bruvo_cc <- bruvo.dist(gc_cc, replen = RL(gc_cc))

    amova_res <- poppr.amova(gc_cc, as.formula(formula_str),
                              dist = bruvo_cc, method = "ade4")
    cat("\nResultados AMOVA:\n");      print(amova_res$results)
    cat("\nComponentes de varianza:\n"); print(amova_res$componentsofcovariance)

    set.seed(seed)
    sig <- randtest(amova_res, nrepet = n_perm)
    cat(paste0("\nSignificancia (", n_perm, " permutaciones):\n"))
    print(sig)

    write.csv(amova_res$results,
              file.path(outdir, paste0("amova_", prefix, "_",
                                        gsub("[~/ ]","",formula_str), ".csv")),
              row.names = TRUE)

    invisible(list(amova = amova_res, signif = sig))
  }, error = function(e) {
    message("AMOVA error (", formula_str, "): ", e$message)
    NULL
  })
}

# ── Índice de asociación con permutaciones ────────────────────────────────────
# Llama a ia() explícitamente con n_perm permutaciones para evaluar la
# significancia del Ia y el rD (standardized index of association).
# Ia ≈ 0 y p > 0.05 → compatible con recombinación libre (panmixis)
# Ia > 0 y p < 0.05 → asociación alélica significativa, consistente con clonalidad
ia_test <- function(gc, prefix, n_perm = 999, seed = 42, plot = FALSE, outdir = ".") {
  cat(paste0("\n── Índice de asociación [", prefix, "] — ",
             n_perm, " permutaciones (n = ", nInd(gc), " aislamientos) ──\n"))
  cat("  arrancando... (esto puede tardar varios minutos con n_perm alto)\n")
  set.seed(seed)
  t0 <- Sys.time()
  out <- tryCatch({
    res <- ia(gc, sample = n_perm, plot = plot, quiet = FALSE)
    cat(sprintf("Ia  = %.4f  p = %.4f\n", res["Ia"],   res["p.Ia"]))
    cat(sprintf("rD  = %.4f  p = %.4f\n", res["rbarD"], res["p.rD"]))
    write.csv(t(as.data.frame(res)),
              file.path(outdir, paste0("ia_", prefix, ".csv")), row.names = TRUE)
    invisible(res)
  }, error = function(e) {
    message("ia() error: ", e$message); NULL
  })
  cat(sprintf("  terminado en %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
  out
}

# ── DAPC ──────────────────────────────────────────────────────────────────────
dapc_plot <- function(gen, pop_vec, colors_vec, n_pca = 50, seed = 42,
                      min_grp_size = 3) {
  grp_sizes  <- table(pop_vec)
  valid_grps <- names(grp_sizes[grp_sizes >= min_grp_size])
  tiny_grps  <- names(grp_sizes[grp_sizes < min_grp_size])
  keep       <- pop_vec %in% valid_grps

  if (length(tiny_grps) > 0)
    message("DAPC: grupos con <", min_grp_size, " muestras excluidos: ",
            paste(tiny_grps, collapse = ", "))

  gen     <- gen[keep, ]
  pop_vec <- pop_vec[keep]
  pop(gen) <- factor(pop_vec)
  grps     <- levels(factor(pop_vec))

  if (length(grps) < 2) {
    message("DAPC: solo hay 1 grupo ('", grps, "') tras filtrar — no hay ",
            "nada que discriminar, se omite este DAPC.")
    return(invisible(NULL))
  }

  n_da     <- length(grps) - 1

  # Loci sin datos (o casi) en este subconcreto de individuos hacen que
  # tab(x, NA.method="mean") interna de dapc() imputen NaN (media de puros
  # NA) -> svd() revienta con "infinite or missing values in 'x'". Se
  # descartan acá con el mismo criterio (5% de missingness) que ya aplica
  # poppr.amova() internamente.
  gen <- missingno(gen, type = "loci", cutoff = 0.05, quiet = TRUE)

  min_grp <- min(table(factor(pop_vec)))
  n_pca   <- min(n_pca, nInd(gen) - length(grps), min_grp - 1L)
  n_pca   <- max(n_pca, 1L)

  set.seed(seed)
  dp <- dapc(gen, n.pca = n_pca, n.da = n_da)

  sc <- as.data.frame(dp$ind.coord)
  if (ncol(sc) >= 2) {
    colnames(sc)[1:2] <- c("DF1","DF2")
  } else {
    colnames(sc)[1] <- "DF1"; sc$DF2 <- 0
  }
  sc$Group <- dp$grp

  var_exp <- round(dp$eig / sum(dp$eig) * 100, 1)

  p <- ggplot(sc, aes(x = DF1, y = DF2, fill = Group, color = Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95,
                 aes(fill = Group), show.legend = FALSE) +
    stat_ellipse(level = 0.95, linewidth = 0.7) +
    geom_point(shape = 21, size = 2,
               color = "white", stroke = 0.3, alpha = 0.85) +
    scale_fill_manual(values  = colors_vec, name = "Group") +
    scale_color_manual(values = colors_vec, name = "Group") +
    labs(x = paste0("DF1 (", var_exp[1], "%)"),
         y = paste0("DF2 (", var_exp[2], "%)")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "right",
          legend.key       = element_blank()) +
    guides(fill  = guide_legend(override.aes = list(size = 4, alpha = 1)),
           color = "none")

  print(p)
  invisible(dp)
}

# ── Barplot: N isolates per year, stacked by country ─────────────────────────
barplot_samples_year_country <- function(df_in, title = NULL, top_n = 12) {
  dat <- df_in %>%
    mutate(Year = as.integer(Year)) %>%
    filter(!is.na(Year))

  top_countries <- names(sort(table(dat$Country), decreasing = TRUE))[1:min(top_n, length(unique(dat$Country)))]
  dat$Country_plot <- ifelse(dat$Country %in% top_countries, dat$Country, "Other")

  dat_count <- dat %>%
    count(Year, Country_plot) %>%
    mutate(Country_plot = factor(Country_plot, levels = c(top_countries, "Other")))

  pal <- c(country_colors[top_countries[top_countries %in% names(country_colors)]],
           setNames(rep("#cccccc", sum(!top_countries %in% names(country_colors))),
                    top_countries[!top_countries %in% names(country_colors)]),
           "Other" = "#aaaaaa")

  ggplot(dat_count, aes(x = Year, y = n, fill = Country_plot)) +
    geom_bar(stat = "identity", color = "grey30", linewidth = 0.15) +
    scale_x_continuous(breaks = seq(min(dat_count$Year, na.rm = TRUE),
                                    max(dat_count$Year, na.rm = TRUE), by = 2)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_fill_manual(values = pal, name = "Country") +
    labs(x = "Year", y = "N isolates", title = title) +
    theme_bw(base_size = 12) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          plot.title      = element_text(hjust = 0.5, face = "bold"))
}

# ── Barplot top N MLGs by genotype/country ────────────────────────────────────
barplot_mlg_country <- function(sc, colors_vec, top_n = 20) {
  top_mlg <- names(sort(table(sc$MLG), decreasing = TRUE))[1:min(top_n, length(unique(sc$MLG)))]
  dat <- sc %>% filter(MLG %in% top_mlg) %>%
    count(MLG, Country) %>% mutate(MLG = factor(MLG, levels = top_mlg))
  ggplot(dat, aes(x = MLG, y = n, fill = Country)) +
    geom_bar(stat = "identity", color = "grey30", linewidth = 0.2) +
    geom_hline(yintercept = 10, linetype = "dotted", linewidth = 0.7) +
    scale_fill_manual(values = colors_vec, name = "Group") +
    scale_y_continuous(breaks = seq(0, 500, by = 10)) +
    labs(x = "MLG", y = "N isolates") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
}

# ── Barplot top N MLGs by year ────────────────────────────────────────────────
barplot_mlg_year <- function(sc, top_n = 20) {
  top_mlg <- names(sort(table(sc$MLG), decreasing = TRUE))[1:min(top_n, length(unique(sc$MLG)))]
  dat <- sc %>% filter(MLG %in% top_mlg) %>%
    count(MLG, Year) %>%
    mutate(MLG = factor(MLG, levels = top_mlg), Year = as.numeric(Year))
  ggplot(dat, aes(x = MLG, y = n, fill = Year)) +
    geom_bar(stat = "identity", color = "grey30", linewidth = 0.2) +
    geom_hline(yintercept = 10, linetype = "dotted", linewidth = 0.7) +
    scale_fill_gradientn(colours = c("blue","yellow","red"),
                         limits = range(as.numeric(sc$Year), na.rm = TRUE),
                         name = "Year") +
    scale_y_continuous(breaks = seq(0, 500, by = 10)) +
    labs(x = "MLG", y = "N isolates") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
}

# ── Barplot N isolates per year, stacked by MLG ───────────────────────────────
barplot_mlg_year_by_country <- function(sc, groups, top_n = 20) {
  dat <- sc[sc$Country %in% groups, ]
  top_mlg <- names(sort(table(dat$MLG), decreasing = TRUE))[1:min(top_n, length(unique(dat$MLG)))]
  dat$MLG_plot <- ifelse(dat$MLG %in% top_mlg, as.character(dat$MLG), "Other")
  dat_count <- dat %>%
    count(Year, MLG_plot) %>%
    mutate(Year     = as.integer(Year),
           MLG_plot = factor(MLG_plot, levels = c(top_mlg, "Other")))
  ggplot(dat_count, aes(x = Year, y = n, fill = MLG_plot)) +
    geom_bar(stat = "identity", color = "grey30", linewidth = 0.2) +
    scale_x_continuous(breaks = seq(min(dat_count$Year, na.rm = TRUE),
                                    max(dat_count$Year, na.rm = TRUE), by = 2)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    labs(x = "Year", y = "N isolates", fill = "MLG",
         title = paste(groups, collapse = " + ")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
}

# ── find.clusters + xvalDapc ─────────────────────────────────────────────────
dapc_cluster_validation <- function(gen, pop_country,
                                     hyp_groups,
                                     n_pca   = 50,
                                     max_k   = 8,
                                     n_rep   = 100,
                                     seed    = 42) {
  set.seed(seed)
  hyp_label <- paste(hyp_groups, collapse = " + ")

  grp_sizes <- table(pop_country)
  keep      <- pop_country %in% names(grp_sizes[grp_sizes >= 3])
  gen_fc    <- gen[keep, ]
  n_pca_fc  <- min(n_pca, nInd(gen_fc) - 1L)

  cat("\n── find.clusters (unsupervised) ──\n")
  fc <- find.clusters(gen_fc,
                      max.n.clust    = max_k,
                      n.pca          = n_pca_fc,
                      choose.n.clust = FALSE,
                      criterion      = "diffNgroup",
                      n.iter         = 1e5,
                      n.start        = 10)

  k_opt <- nlevels(fc$grp)
  cat("Optimal k:", k_opt, "\n")

  bic_df <- data.frame(k = seq_along(fc$Kstat), BIC = fc$Kstat)
  p_bic <- ggplot(bic_df, aes(x = k, y = BIC)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(size = 3, color = "steelblue") +
    geom_vline(xintercept = k_opt, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text", x = k_opt + 0.15, y = max(bic_df$BIC),
             label = paste0("k = ", k_opt), hjust = 0, color = "red", size = 4) +
    labs(title = "Unsupervised clustering — BIC by k",
         x = "Number of clusters (k)", y = "BIC") +
    theme_bw(base_size = 12)
  print(p_bic)

  ref <- ifelse(pop_country[keep] %in% hyp_groups, hyp_label, "Other")
  cat("\nCross-tabulation vs '", hyp_label, "':\n", sep = "")
  print(table(Cluster = fc$grp, Group = ref))

  cat("\n── xvalDapc — cross-validation of '", hyp_label, "' vs rest ──\n", sep = "")
  pop_binary <- factor(ifelse(pop_country %in% hyp_groups, hyp_label, "Other"))
  pop(gen)   <- pop_binary
  x          <- tab(gen, freq = TRUE, NA.method = "mean")
  n_pca_xv   <- min(n_pca, nInd(gen) - 2L)

  set.seed(seed)
  xval <- xvalDapc(x, pop(gen),
                   n.pca.max    = n_pca_xv,
                   training.set = 0.9,
                   result       = "groupMean",
                   center       = TRUE,
                   scale        = FALSE,
                   n.rep        = n_rep,
                   xval.plot    = TRUE)

  best_npca <- xval$`Number of PCs Achieving Highest Mean Success`
  best_acc  <- round(max(xval$`Mean Successful Assignment by Number of PCs of PCA`) * 100, 1)
  cat("Optimal n.pca:", best_npca, "— % correct assignment:", best_acc, "%\n")

  invisible(list(find_clusters = fc, xval = xval))
}

# ── Export diversity tables ───────────────────────────────────────────────────
# Every step is individually wrapped in tryCatch (same defensive pattern as
# amova_analysis/ia_test/safe_poppr elsewhere in this script). A failure in
# any one piece (locus table, a genclone construction, a poppr() call, a
# write.csv) prints a message and yields NULL/skips, but does NOT abort the
# whole function — which matters because export_diversity() runs inside
# run_analysis(), where an uncaught error kills the entire run_analysis()
# call and silently leaves res1/res2 unassigned (visible only as accumulated
# warnings with no clear error text, then "object 'resX' not found" on the
# next line that tries to use it).
export_diversity <- function(gc, gm, df_in, sn, prefix, outdir = ".") {
  safe_run <- function(expr, label) {
    tryCatch(expr, error = function(e) {
      message("export_diversity(): '", label, "' failed — ", e$message)
      NULL
    })
  }

  div_locus <- safe_run(locus_table(gc, info = FALSE), "locus_table")

  # df_in ya está alineado 1:1 con gm/sn (rownames(df_in) <- sn en
  # prep_subset()), así que las columnas se usan directo. OJO: df_in$col[sn]
  # (indexar por nombre un vector plano extraído con $) NO funciona — ese
  # vector no hereda los rownames de df_in como names(), así que devuelve
  # NA para todos los individuos. Eso era lo que producía genclones con
  # población NA y rompía el formateo interno de poppr() más abajo.
  gc_country <- safe_run(as.genclone(df2genind(gm, sep = "/", ploidy = 3,
                             pop = df_in$Pop_country, ind.names = sn)), "gc_country")
  gc_year    <- safe_run(as.genclone(df2genind(gm, sep = "/", ploidy = 3,
                             pop = as.character(df_in$Year), ind.names = sn)), "gc_year")
  gc_geno    <- safe_run(as.genclone(df2genind(gm, sep = "/", ploidy = 3,
                             pop = df_in$Pop_geno, ind.names = sn)), "gc_geno")

  safe_poppr <- function(x) {
    if (is.null(x)) return(NULL)
    tryCatch(poppr(x),
             error = function(e) { message("poppr() error: ", e$message); NULL })
  }

  div_country <- safe_poppr(gc_country)
  div_year    <- safe_poppr(gc_year)
  div_geno    <- safe_poppr(gc_geno)

  if (!is.null(div_locus)) {
    cat("\n── Diversity per locus ──\n")
    safe_run(print(round(as.data.frame(div_locus), 3)), "print div_locus")
  }
  if (!is.null(div_geno))    { cat("\n── Diversity per genotype ──\n"); print(div_geno) }
  if (!is.null(div_country)) { cat("\n── Diversity per country ──\n");  print(div_country) }
  if (!is.null(div_year))    { cat("\n── Diversity per year ──\n");     print(div_year) }

  if (!is.null(div_locus))
    safe_run(write.csv(as.data.frame(div_locus),
              file.path(outdir, paste0("diversity_",prefix,"_per_locus.csv"))),
             "write div_locus")
  if (!is.null(div_geno))
    safe_run(write.csv(div_geno[, !names(div_geno) == "File"],
              file.path(outdir, paste0("diversity_",prefix,"_per_genotype.csv")), row.names = FALSE),
              "write div_geno")
  if (!is.null(div_country))
    safe_run(write.csv(div_country[, !names(div_country) == "File"],
              file.path(outdir, paste0("diversity_",prefix,"_per_country.csv")), row.names = FALSE),
              "write div_country")
  if (!is.null(div_year))
    safe_run(write.csv(div_year[, !names(div_year) == "File"],
              file.path(outdir, paste0("diversity_",prefix,"_per_year.csv")), row.names = FALSE),
              "write div_year")

  invisible(list(div_locus = div_locus, div_country = div_country,
                 div_year = div_year, div_geno = div_geno))
}

# ── Map helpers: shared prep + drawing, reused by both category maps below ────
# prep_category_map_data(): counts isolates per Country x category (cat_fn
# decides the category per isolate), plus the oldest/newest sampling year per
# country (for the map label). cat_fn receives the already country-filtered
# df_in and must return a character vector the same length as nrow(df_in).
prep_category_map_data <- function(df_in, cat_fn, cat_levels) {
  df_country <- df_in %>% filter(Country %in% map_todos_paises)

  df_counts <- df_country %>%
    mutate(cat = cat_fn(df_country)) %>%
    count(Country, cat) %>%
    pivot_wider(names_from = cat, values_from = n, values_fill = 0)

  for (cl in cat_levels) if (!cl %in% names(df_counts)) df_counts[[cl]] <- 0L

  # Oldest/newest sampled isolate per country -> map label, e.g.
  # "Argentina\n(1999-2026)."
  year_range <- df_country %>%
    filter(!is.na(Year)) %>%
    group_by(Country) %>%
    summarise(year_min = min(Year), year_max = max(Year), .groups = "drop")

  df_counts %>%
    mutate(total = rowSums(across(all_of(cat_levels)))) %>%
    left_join(map_coords, by = "Country") %>%
    left_join(year_range, by = "Country") %>%
    filter(!is.na(lon), total > 0) %>%
    mutate(r = rescale(sqrt(total), to = c(1.5, 5)),
           year_label = ifelse(
             is.na(year_min), paste0(Country, "\n(Unknown)."),
             paste0(Country, "\n(", year_min, "-", year_max, ").")
           ))
}

# draw_category_map(): the actual ggplot (base map + circles/pies + count
# label + optional country/year-range label + legend/theme). Identical for
# every category map — only df_pie's category columns, cat_colors and
# show_year_label change between calls, so this is the single place that
# needs editing for any shared styling change (label size, offset, map
# extent, etc.).
draw_category_map <- function(df_pie, cat_levels, cat_colors, legend_name,
                               title = NULL, label_offset = 1,
                               show_year_label = TRUE, show_legend = TRUE) {
  n_nonzero <- rowSums(df_pie[cat_levels] > 0)
  df_pure   <- df_pie[n_nonzero == 1, ]
  df_mixed  <- df_pie[n_nonzero > 1, ]

  if (nrow(df_pure) > 0)
    df_pure$cat_pure <- apply(df_pure[cat_levels], 1,
                               function(row) cat_levels[row > 0][1])

  paises_con_datos <- df_pie$Country
  mapa_base        <- map_data("world") %>% filter(region %in% map_todos_paises)

  p <- ggplot() +
    geom_polygon(data = mapa_base %>% filter(!region %in% paises_con_datos),
                 aes(x = long, y = lat, group = group),
                 fill = "grey92", color = "grey75", linewidth = 0.3) +
    geom_polygon(data = mapa_base %>% filter(region %in% paises_con_datos),
                 aes(x = long, y = lat, group = group),
                 fill = "grey82", color = "grey60", linewidth = 0.3)

  if (nrow(df_pure) > 0)
    p <- p + geom_circle(data = df_pure,
                          aes(x0 = lon, y0 = lat, r = r, fill = cat_pure),
                          color = "white", linewidth = 0.5, alpha = 0.92)

  if (nrow(df_mixed) > 0)
    p <- p + geom_scatterpie(data = df_mixed,
                             aes(x = lon, y = lat, r = r),
                             cols = cat_levels,
                             color = "white", linewidth = 0.5, alpha = 0.92)

  p <- p +
    geom_text(data = df_pie,
              aes(x = lon, y = lat, label = total),
              size = 2.8, fontface = "bold", color = "black")

  # Country name + oldest-newest sampling year, just below each pie/circle
  # (offset scales with the pie's own radius so bigger pies don't run the
  # label into the pie below it) — nudge label_offset by eye if two nearby
  # countries' labels collide. Optional (show_year_label): the genotype map
  # skips it to keep that figure less cluttered; the provenance map keeps it.
  if (show_year_label)
    p <- p +
      geom_text(data = df_pie,
                aes(x = lon, y = lat - r - label_offset, label = year_label),
                size = 2.3, fontface = "plain", color = "black", lineheight = 0.85)

  p +
    scale_fill_manual(values = cat_colors, name = legend_name) +
    coord_fixed(ratio = 1, xlim = c(-92, -33), ylim = c(-57, 22)) +
    labs(title = title) +
    theme_void() +
    theme(plot.background  = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.title       = element_text(hjust = 0.5, size = 10, face = "bold"),
          legend.position  = if (show_legend) "bottom" else "none",
          legend.title     = element_text(size = 9),
          legend.text      = element_text(size = 8))
}

# ── Map: genotype proportion per country (post-2005) ──────────────────────────
# Six categories: EU_2_A1 (incl. PE_7.1) / EC1 (incl. EC_1 label variant) /
# EU_13_A2 / US7 (incl. US_7 / US_7_A2 label variants) / US_1 / Other.
# Countries with a single category present render as a solid circle;
# countries with more than one category present render as a pie
# (geom_scatterpie). cat_colors is a parameter so the palette (or the set of
# categories, together with the case_when below) can be adjusted later
# without touching the rest of the function. Deliberately its own smaller
# category set (no PE_3 / BR_1+EU_8_A1 slice) rather than the full
# map_genotype_category() used elsewhere — this map is meant to stay
# readable with only the handful of genotypes that matter for the
# post-2005 country-composition story.
make_eu2a1_map <- function(df_in, title = NULL,
                            # same colors as every other genotype-colored figure
                            cat_colors = map_geno_colors[c("EU_2_A1", "EC1", "EU_13_A2",
                                                           "US7", "US_1", "Other")],
                            show_legend = TRUE) {
  cat_levels <- names(cat_colors)

  cat_fn <- function(d) case_when(
    d$Genotype %in% c("EU_2_A1","PE_7.1")             ~ "EU_2_A1",
    d$Genotype %in% c("EC1","EC_1")                   ~ "EC1",
    d$Genotype == "EU_13_A2"                          ~ "EU_13_A2",
    d$Genotype %in% c("US7","US_7","US_7_A2")         ~ "US7",
    d$Genotype == "US_1"                              ~ "US_1",
    TRUE                                                ~ "Other"
  )

  df_pie <- prep_category_map_data(df_in, cat_fn, cat_levels)
  draw_category_map(df_pie, cat_levels, cat_colors, legend_name = "Genotype",
                     title = title, show_year_label = FALSE, show_legend = show_legend)
}

# ── Map: sample provenance per country (post-2005) ────────────────────────────
# Two categories: "This work" (isolates sourced from this study — the Source
# column carries "This work"/"this work", and the legacy label
# "Flor_Lucca_check" used in some earlier tables for the same isolates) vs.
# "Context samples" (every other isolate: public databases, Euroblight
# reference panel, other published studies). Matching is on the Source column
# if present, falling back to lowercase "source" (used in some older tables),
# and is case-insensitive so future label variants of the same two sources
# don't need a script edit. Same drawing logic as make_eu2a1_map() via
# draw_category_map(), including the country/year-range label (kept here,
# unlike the genotype map).
make_provenance_map <- function(df_in, title = NULL,
                                 cat_colors = c("This work"       = "#6a3d9a",
                                                "Context samples" = "#FFD700")) {
  cat_levels <- names(cat_colors)

  cat_fn <- function(d) {
    src <- if ("Source" %in% names(d)) d$Source else d$source
    provenance_category(src)
  }

  df_pie <- prep_category_map_data(df_in, cat_fn, cat_levels)
  draw_category_map(df_pie, cat_levels, cat_colors, legend_name = "Provenance", title = title)
}

# ── Master analysis wrapper ───────────────────────────────────────────────────
# amova_formula: fórmula AMOVA como string. Default ~Pop_geno/Pop_country
#                (lineaje anidado en localidad). Para subsets de un solo
#                linaje usar "~Pop_country".
# do_dapc/do_amova: para S1 (dataset disperso, muchos genotipos/países) se
#                decidió no reportar DAPC ni AMOVA jerárquico — el PCoA
#                descriptivo alcanza para la pregunta de esa sección (ver
#                Resultados) y evita el sobre-ajuste que DAPC/AMOVA
#                introducirían sobre un dataset tan heterogéneo. S2 (EU_2_A1)
#                sí necesita ambos para testear la hipótesis Argentina+Uruguay
#                vs. resto, así que quedan prendidos por default.
run_analysis <- function(s, label, amova_formula = "~Pop_geno/Pop_country",
                          do_dapc = TRUE, do_amova = TRUE,
                          pcoa_ellipses = TRUE) {
  d <- s$df; gm <- s$gm; sn <- s$sn
  # All files this function writes go to outputs/s1 or outputs/s2, matching
  # the folder split by dataset — label is always "S1" or "S2".
  outdir <- file.path("outputs", tolower(label))

  gc_base <- as.genclone(df2genind(gm, sep = "/", ploidy = 3, ind.names = sn))
  bruvo   <- bruvo.dist(gc_base, replen = RL(gc_base))

  # gc_geno/gen_geno are grouped by Pop_geno_map (the unified genotype
  # palette/categories, Section 1), not raw Pop_geno — otherwise the MSN and
  # PCA below would show one node/color per raw label (including
  # never-colored variants like "EC_1", "US_7_A2") instead of the same
  # EU_2_A1/EC1/EU_13_A2/US7/US_1/... categories used by the PCoA, UPGMA and
  # map figures. DAPC is unaffected: dapc_plot() always overrides pop() with
  # whatever pop_vec is explicitly passed to it (still raw d$Pop_geno there).
  gc_geno    <- as.genclone(df2genind(gm, sep = "/", ploidy = 3, pop = d$Pop_geno_map, ind.names = sn))
  gc_country <- as.genclone(df2genind(gm, sep = "/", ploidy = 3, pop = d$Pop_country,  ind.names = sn))
  gc_period  <- as.genclone(df2genind(gm, sep = "/", ploidy = 3, pop = d$Period,       ind.names = sn))

  gen_geno    <- df2genind(gm, sep = "/", ploidy = 3, pop = d$Pop_geno_map, ind.names = sn)
  gen_country <- df2genind(gm, sep = "/", ploidy = 3, pop = d$Pop_country,  ind.names = sn)

  # ── MSN ────────────────────────────────────────────────────────────────────
  cat(paste0("\n══ MSN — ", label, " ══\n"))
  msn_plot(gc_geno,    bruvo, map_geno_colors, seed = 23)
  msn_plot(gc_country, bruvo, country_colors,  seed = 23)
  msn_plot(gc_period,  bruvo, period_colors,   seed = 23)

  # ── PCA (frecuencias alélicas, ade4::dudi.pca) ─────────────────────────────
  cat(paste0("\n══ PCA — ", label, " ══\n"))
  pca_g <- pca_scores(gen_geno,    gc_geno,    color_vec = d$Pop_geno_map, year_vec = d$Year)
  pca_c <- pca_scores(gen_country, gc_country, color_vec = d$Pop_country, year_vec = d$Year)

  cat(paste0("Explained variance (%):\n")); print(pca_g$var[1:5])

  print(pca_2d(pca_g$sc_col, map_geno_colors, pca_g$var))
  print(pca_2d(pca_g$sc_col, map_geno_colors, pca_g$var, "PC1","PC3"))
  pca_3d(pca_g$sc_col, map_geno_colors, pca_g$var)
  print(pca_2d(pca_c$sc_col, country_colors, pca_c$var))
  pca_3d(pca_c$sc_col, country_colors, pca_c$var)

  # ── PCoA (distancias de Bruvo, cmdscale) ───────────────────────────────────
  # Uses pcoa_publication_plot() (not the older pcoa_bruvo()) so these
  # exploratory PCoAs get the same treatment as the publication panels: one
  # point per MLG (not one per raw individual — avoids stacking identical
  # clones on top of each other), sized by how many isolates share that MLG,
  # with a 95% ellipse per group. Genotype-colored version uses Pop_geno_map
  # (the unified genotype palette, Section 1) instead of raw Pop_geno, so
  # EU_13_A2/EC1/EU_2_A1/US7/US_1 render in their real colors instead of
  # falling back to grey (raw Pop_geno keeps label variants like "EC_1"
  # that geno_colors never had an entry for).
  cat(paste0("\n══ PCoA — ", label, " ══\n"))
  pcoa_g <- pcoa_publication_plot(
    bruvo_dist          = bruvo, df_in = d, gc = gc_base,
    point_group_col     = "Pop_geno_map",
    point_colors_vec    = map_geno_colors,
    point_legend_name   = "Genotype",
    ellipse_group_col   = if (pcoa_ellipses) "Pop_geno_map" else NULL,
    ellipse_colors_vec  = map_geno_colors,
    filename            = file.path(outdir, paste0(label, "_PCoA_by_geno.png"))
  )
  pcoa_c <- pcoa_publication_plot(
    bruvo_dist            = bruvo, df_in = d, gc = gc_base,
    point_group_col       = "Pop_country",
    point_colors_vec      = country_colors,
    point_legend_name     = "Country",
    ellipse_group_col     = if (pcoa_ellipses) "Pop_geno_map" else NULL,
    ellipse_linetypes_vec = map_geno_linetypes,
    filename              = file.path(outdir, paste0(label, "_PCoA_by_country.png"))
  )

  # ── UPGMA ──────────────────────────────────────────────────────────────────
  # Both groupings are written out, and the filename says which one colors
  # the tree, so the two files are never ambiguous.
  # By-genotype tree colors tips via Pop_geno_map/map_geno_colors (same
  # unified palette as the PCoA above and the map figures), not raw
  # Pop_geno/geno_colors — with the raw labels most tips ("EC_1", "EU_13_A2",
  # "US_7_A2", …) have no matching color entry and render grey.
  cat(paste0("\n══ UPGMA — ", label, " ══\n"))
  tree_geno    <- upgma_plot(bruvo, d, map_geno_colors, "Pop_geno_map", label = label,
                              filename = file.path(outdir, paste0("UPGMA_", label, "_by_geno")))
  tree_country <- upgma_plot(bruvo, d, country_colors, "Pop_country", label = label,
                              filename = file.path(outdir, paste0("UPGMA_", label, "_by_country")))

  # ── DAPC ───────────────────────────────────────────────────────────────────
  if (do_dapc) {
    cat(paste0("\n══ DAPC — ", label, " ══\n"))
    dapc_plot(gen_geno,    d$Pop_geno,    geno_colors)
    dapc_plot(gen_country, d$Pop_country, country_colors)
  }

  # ── Índice de asociación (Ia, rD) con 999 permutaciones ───────────────────
  # Clone-corrected within each country before testing: gc_country still has
  # one row per raw isolate, and S1 in particular is dominated by long runs
  # of identical EU_2_A1 MLGs within the same country (clonal reproduction).
  # Testing on the raw (non-clone-corrected) data both (a) inflates apparent
  # linkage disequilibrium — repeated copies of the same genotype trivially
  # look like "non-random association" — and (b) is why this permutation
  # step is so much slower than everything else: it's recomputing the
  # statistic 999 times over every duplicate clone in n = nrow(d), instead of
  # over the much smaller set of distinct genotypes per country.
  cat(paste0("\n══ Index of association — ", label, " ══\n"))
  gc_country_cc <- gc_country[!duplicated(data.frame(Pop = pop(gc_country), MLG = mll(gc_country))), ]
  cat("Clone-corrected (per country) from", nInd(gc_country), "to", nInd(gc_country_cc), "individuals.\n")
  ia_test(gc_country_cc, prefix = label, n_perm = 999, outdir = outdir)

  # ── Diversity ──────────────────────────────────────────────────────────────
  cat(paste0("\n══ Diversity — ", label, " ══\n"))
  export_diversity(gc_base, gm, d, sn, label, outdir = outdir)

  # ── AMOVA ──────────────────────────────────────────────────────────────────
  amova_res <- NULL
  if (do_amova) {
    cat(paste0("\n══ AMOVA — ", label, " ══\n"))
    # Determinar columnas de estratos a partir de la fórmula
    strata_cols <- trimws(unlist(strsplit(gsub("~","", amova_formula), "/")))
    amova_res <- amova_analysis(gm, d, sn,
                                 formula_str = amova_formula,
                                 strata_cols = strata_cols,
                                 prefix      = label,
                                 n_perm      = 999,
                                 outdir      = outdir)
  }

  # ── Barplots ───────────────────────────────────────────────────────────────
  cat(paste0("\n══ Barplots MLG — ", label, " ══\n"))
  tryCatch(print(barplot_mlg_country(pca_g$sc, map_geno_colors)),
           error = function(e) message("barplot_mlg_country() failed — ", e$message))
  tryCatch(print(barplot_mlg_year(pca_g$sc)),
           error = function(e) message("barplot_mlg_year() failed — ", e$message))

  invisible(list(pca_g      = pca_g,
                 pca_c      = pca_c,
                 pcoa_g     = pcoa_g,
                 pcoa_c     = pcoa_c,
                 tree_geno  = tree_geno,
                 tree_country = tree_country,
                 gc_base    = gc_base,
                 bruvo      = bruvo,
                 gc_geno    = gc_geno,
                 gc_country = gc_country,
                 amova      = amova_res))
}

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 4 — ANALYSIS 1: EU_2_A1 + EC1(100) + EU_8_A1 + BR_1 + PE_7.1
# ══════════════════════════════════════════════════════════════════════════════

# S1 (global South American dataset): PCoA without ellipses (see Section 7).
res1 <- run_analysis(s1, "S1", do_dapc = FALSE, do_amova = FALSE,
                     pcoa_ellipses = FALSE)

print(barplot_mlg_year_by_country(res1$pca_c$sc, "Argentina"))
print(barplot_mlg_year_by_country(res1$pca_c$sc, c("Argentina","Uruguay")))
print(barplot_mlg_year_by_country(res1$pca_g$sc, "EU_2_A1"))

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 5 — ANALYSIS 2: EU_2_A1 + PE_7.1
# ══════════════════════════════════════════════════════════════════════════════
# Para S2 se mantiene ~Pop_geno/Pop_country: permite evaluar si PE_7.1 y EU_2_A1
# difieren genotípicamente (hipótesis nula: son el mismo linaje).
res2 <- run_analysis(s2, "S2", amova_formula = "~Pop_geno/Pop_country")

# DAPC S2: Argentina + Uruguay (celeste) vs resto (rojo)
s2_countries      <- unique(s2$df$Pop_country)
s2_arg_uru_colors <- setNames(
  ifelse(s2_countries %in% c("Argentina","Uruguay"), "skyblue", "red"),
  s2_countries
)
gen_s2_country <- df2genind(s2$gm, sep = "/", ploidy = 3,
                             pop = s2$df$Pop_country, ind.names = s2$sn)
dapc_plot(gen_s2_country, s2$df$Pop_country, s2_arg_uru_colors)

print(barplot_mlg_year_by_country(res2$pca_c$sc, "Argentina"))
print(barplot_mlg_year_by_country(res2$pca_c$sc, c("Argentina","Uruguay")))

sc_arg_uru <- res2$pca_c$sc %>%
  filter(Country %in% c("Argentina","Uruguay","Chile","Brazil","Bolivia","Peru","Panama"))
print(barplot_mlg_country(sc_arg_uru, country_colors))

# Validación de la hipótesis Argentina + Uruguay vs resto
gen_s2_base <- df2genind(s2$gm, sep = "/", ploidy = 3, ind.names = s2$sn)
val_s2 <- dapc_cluster_validation(
  gen         = gen_s2_base,
  pop_country = s2$df$Pop_country,
  hyp_groups  = c("Argentina","Uruguay"),
  n_pca       = 50,
  max_k       = 8,
  n_rep       = 100
)

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 6 — MAPS: genotype AND provenance distribution (all isolates)
# ══════════════════════════════════════════════════════════════════════════════
# Map A: EU_2_A1(+PE_7.1) / EC1 / EU_13_A2 / US7 / US_1 / Other, over every
# isolate of the table (no filter on sampling year; isolates without a year
# are included). No legend: in Figure 1 the genotype colors are explained by
# the legend of the tree (same palette, map_geno_colors).
p_map <- make_eu2a1_map(df1, title = "", show_legend = FALSE)
print(p_map)
save_ggplot_multi(p_map, "outputs/s1/mapa_EU2A1_EC1_post2005",
                   width = 6, height = 8, dpi = 600, bg = "transparent")

# Map B: same layout/labels, colored by sample provenance instead of
# genotype — "This work" (purple) vs. "Context samples" (grey).
p_map_provenance <- make_provenance_map(df1, title = "")
print(p_map_provenance)
save_ggplot_multi(p_map_provenance, "outputs/s1/mapa_provenance_post2005",
                   width = 6, height = 8, dpi = 600, bg = "transparent")

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 7 — FIGURE 1b: PCoA of dataset 1, points by country (final figure)
#  Standalone block — does NOT depend on Section 4 / run_analysis(s1, ...)
#  having been run first, so this panel can be regenerated on its own.
#
#  PCoA (Bruvo distance, cmdscale) on dataset 1
#  (EU_2_A1 + EC1(100) + EU_8_A1 + BR_1 + PE_7.1). One point per MLG
#  (isolates sharing an identical multilocus genotype collapse to the same
#  coordinates), sized by how many isolates share that MLG, colored by
#  country — MLGs whose members come from more than one country are drawn
#  in black ("Multiple").
#
#  No ellipses: the three dominant lineages (EU_2_A1, EC1, EU_13_A2) already
#  separate visibly along PCo1/PCo2, and the genotype-colored panel next to
#  this one (same axes, same panel size) shows the lineage of every point, so
#  ellipses would add no information and could suggest a statistical
#  separation that this descriptive figure does not test.
#
#  Both panels use the same circle rendering (radius = f(n isolates), pies for
#  MLGs spanning >1 group), the same axis limits and the same panel size, so
#  they can be compared point by point.
# ══════════════════════════════════════════════════════════════════════════════
gc_base_s1 <- as.genclone(df2genind(s1$gm, sep = "/", ploidy = 3, ind.names = s1$sn))
bruvo_s1   <- bruvo.dist(gc_base_s1, replen = RL(gc_base_s1))

df_fig1b <- s1$df

# Both Fig1b panels share axis limits, panel size and canvas size so they line
# up exactly when placed side by side (make_figura1.R scales each PNG to its
# box). Limits come from the same PCoA the function computes, plus padding.
pc_s1 <- cmdscale(as.matrix(bruvo_s1), k = 2)
fig1b_axis_limits <- list(x = range(pc_s1[, 1]) + c(-0.03, 0.03),
                          y = range(pc_s1[, 2]) + c(-0.03, 0.03))
fig1b_panel_w  <- 4.4   # inches; height follows from the 1:1 axis ratio
fig1b_panel_size <- c(fig1b_panel_w,
                      fig1b_panel_w * diff(fig1b_axis_limits$y) / diff(fig1b_axis_limits$x))
fig1b_canvas <- c(5.4, 5)  # inches (width, height), identical for both panels; trimmed to the content (no legend)

fig1b_country <- pcoa_publication_plot(
  bruvo_dist            = bruvo_s1,
  df_in                 = df_fig1b,
  gc                    = gc_base_s1,
  point_group_col       = "Pop_country",
  point_colors_vec      = country_colors,
  point_legend_name     = "Country",
  point_style           = "pie",   # MLGs mixing >1 country render as pie charts
  pie_radius_range      = c(0.003, 0.018),   # shared with fig1b_genotype below
  point_alpha           = 0.75,
  ellipse_group_col     = NULL,    # no ellipses
  show_legend           = FALSE,   # colors are explained by the tree legends (same palettes)
  axis_limits           = fig1b_axis_limits,
  panel_size            = fig1b_panel_size,
  width = fig1b_canvas[1], height = fig1b_canvas[2],
  filename              = "outputs/s1/Fig1b_by_country.png"
)

# Genotype-colored version (Pop_geno_map, palette map_geno_colors — the same
# colors as the map figures (Section 9) and the UPGMA trees).
fig1b_genotype <- pcoa_publication_plot(
  bruvo_dist          = bruvo_s1,
  df_in               = df_fig1b,
  gc                  = gc_base_s1,
  point_group_col     = "Pop_geno_map",
  point_colors_vec    = map_geno_colors,
  point_legend_name   = "Genotype",
  point_style         = "pie",   # same circle/pie rendering and radius as fig1b_country,
  pie_radius_range    = c(0.003, 0.018),   # so point sizes match between the two panels
  point_alpha         = 0.75,
  ellipse_group_col   = NULL,    # no ellipses
  show_legend         = FALSE,   # colors are explained by the tree legends (same palettes)
  axis_limits         = fig1b_axis_limits,
  panel_size          = fig1b_panel_size,
  width = fig1b_canvas[1], height = fig1b_canvas[2],
  filename            = "outputs/s1/Fig1b_by_genotype.png"
)

# ── FIGURE 1 (left third): UPGMA tree with a color "table" at the leaves ─────
# UPGMA tree (Bruvo distance) of every isolate in dataset 1. No tip labels:
# each leaf carries three color columns instead — Country, Genotype and Source
# (This work vs. Context samples) — using the same palettes as every other
# figure (country_colors, map_geno_colors) and blue/red for Source, with the
# three legends side by side below the tree. Drawn directly with
# ggplot2 (segments + rectangles + in-plot legends) rather than with
# ggtree::gheatmap(): the installed ggtree 3.8.2 is not compatible with
# ggplot2 4.x. Saved at 4 x 10 in = one third of the 12 x 10 in Figure 1 page,
# so make_figura1.R places it at its native size.
upgma_color_table_plot <- function(tree, meta, col_palettes,
                                    col_w = 0.08, col_gap = 0.015,
                                    legend_text_size = 3.2, fig_width_in = 4,
                                    legends_below = FALSE) {
  tree <- ape::reorder.phylo(ape::ladderize(tree), "cladewise")
  n    <- length(tree$tip.label)
  N    <- n + tree$Nnode

  # Tip y = order of appearance in the cladewise edge matrix (what plot.phylo
  # does); internal node y = midpoint of its children's y, filled in postorder
  # (every edge below a node is visited before the edge above it).
  y <- numeric(N)
  tip_nodes <- tree$edge[tree$edge[, 2] <= n, 2]
  y[tip_nodes] <- seq_along(tip_nodes)
  cmin <- rep(Inf, N); cmax <- rep(-Inf, N)
  po <- ape::reorder.phylo(tree, "postorder")$edge
  for (i in seq_len(nrow(po))) {
    par <- po[i, 1]; chi <- po[i, 2]
    if (chi > n) y[chi] <- (cmin[chi] + cmax[chi]) / 2
    cmin[par] <- min(cmin[par], y[chi]); cmax[par] <- max(cmax[par], y[chi])
  }
  root <- n + 1
  y[root] <- (cmin[root] + cmax[root]) / 2

  x <- ape::node.depth.edgelength(tree)
  h <- max(x)
  flip <- function(v) n + 1 - v          # first tip at the top

  seg_h <- data.frame(x = x[tree$edge[, 1]], xend = x[tree$edge[, 2]],
                      y = flip(y[tree$edge[, 2]]))
  internal <- (n + 1):N
  seg_v <- data.frame(x = x[internal], y = flip(cmin[internal]), yend = flip(cmax[internal]))

  p <- ggplot() +
    geom_segment(data = seg_h, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.15) +
    geom_segment(data = seg_v, aes(x = x, xend = x, y = y, yend = yend), linewidth = 0.15)

  # Color columns: one rectangle per leaf per column.
  ytip <- flip(y[seq_len(n)])
  col_x0 <- h * (1 + col_gap) + (seq_along(col_palettes) - 1) * h * (col_w + col_gap)
  W      <- tail(col_x0, 1) + h * col_w      # total drawn width (x from 0 to W)
  for (k in seq_along(col_palettes)) {
    cname <- names(col_palettes)[k]
    vals  <- as.character(meta[tree$tip.label, cname])
    cols  <- unname(col_palettes[[k]][vals]); cols[is.na(cols)] <- "#cccccc"
    p <- p + geom_rect(data = data.frame(xmin = col_x0[k], xmax = col_x0[k] + h * col_w,
                                         ymin = ytip - 0.5, ymax = ytip + 0.5, fill = cols),
                       aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                       color = NA)
    p <- p + annotate("text", x = col_x0[k] + h * col_w / 2, y = n + 0.012 * n,
                      label = cname, angle = 90, hjust = 0, vjust = 0.5,
                      fontface = "bold", size = legend_text_size + 0.4)
  }

  # Optional legends below the tree (off by default: in Figure 1 the legends are
  # a separate panel, see upgma_color_legend_plot()).
  y_lo <- 0.5 - 0.004 * n
  if (legends_below) {
    # Legends below the tree, one per palette side by side (palette order, only
    # the levels present), justified to the width of the tree: the first legend
    # starts at the tree's left edge and the last one ends at the right edge of
    # the last color column, with the spare space split evenly between legends.
    # Text width is estimated from the number of characters (~0.55 em each).
    step    <- 0.024 * n
    box_h   <- 0.016 * n
    box_w   <- 0.045 * W
    n_leg   <- length(col_palettes)
    y_start <- -0.8 * step
    usable_in <- fig_width_in - 2 * 8 / 72                  # plot.margin = 8 pt
    char_x    <- (legend_text_size * 72.27 / 25.4 * 0.55 / 72) / usable_in * W
    leg_levels <- lapply(seq_len(n_leg), function(k)
      intersect(names(col_palettes[[k]]), unique(as.character(meta[, names(col_palettes)[k]]))))
    leg_w <- vapply(seq_len(n_leg), function(k)
      box_w * 1.45 + max(nchar(c(leg_levels[[k]], names(col_palettes)[k]))) * char_x, numeric(1))
    leg_gap <- if (n_leg > 1) max(0, (W - sum(leg_w)) / (n_leg - 1)) else 0
    leg_x0  <- c(0, cumsum(leg_w + leg_gap))[seq_len(n_leg)]
    max_rows <- 0
    for (k in seq_len(n_leg)) {
      pal    <- col_palettes[[k]]
      levels <- leg_levels[[k]]
      x0     <- leg_x0[k]
      p <- p + annotate("text", x = x0, y = y_start, label = names(col_palettes)[k],
                        hjust = 0, vjust = 0.5, fontface = "bold", size = legend_text_size + 0.4)
      for (j in seq_along(levels)) {
        yy <- y_start - j * step
        p <- p + annotate("rect", xmin = x0, xmax = x0 + box_w,
                          ymin = yy - box_h / 2, ymax = yy + box_h / 2,
                          fill = pal[[levels[j]]], color = "grey20", linewidth = 0.15) +
                 annotate("text", x = x0 + box_w * 1.45, y = yy, label = levels[j],
                          hjust = 0, vjust = 0.5, size = legend_text_size)
      }
      max_rows <- max(max_rows, length(levels))
    }
    y_lo <- y_start - (max_rows + 0.6) * step
  }

  p +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_lo, n * 1.10), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(8, 8, 8, 8), plot.background = element_rect(fill = "white", color = NA))
}

# Legends of the tree (Country / Genotype / Source) as a stand-alone panel, one
# column of color keys, drawn in inch units so it can be placed anywhere in
# Figure 1 at its native size. They are the single color reference for the
# whole figure (map, PCoAs and tree share the palettes).
upgma_color_legend_plot <- function(meta, col_palettes, width_in = 1.7,
                                     text_size = 3.2, step = 0.19, box = 0.12) {
  levs <- lapply(seq_along(col_palettes), function(k)
    intersect(names(col_palettes[[k]]), unique(as.character(meta[, names(col_palettes)[k]]))))
  title_gap <- 0.6 * step
  n_rows    <- sum(lengths(levs)) + length(levs)           # entries + titles
  height_in <- n_rows * step + (length(levs) - 1) * title_gap + 0.1

  p <- ggplot()
  y <- height_in - 0.05 - step / 2
  for (k in seq_along(col_palettes)) {
    p <- p + annotate("text", x = 0.02, y = y, label = names(col_palettes)[k],
                      hjust = 0, vjust = 0.5, fontface = "bold", size = text_size + 0.4)
    y <- y - step
    for (lv in levs[[k]]) {
      p <- p + annotate("rect", xmin = 0.02, xmax = 0.02 + box, ymin = y - box / 2,
                        ymax = y + box / 2, fill = col_palettes[[k]][[lv]],
                        color = "grey20", linewidth = 0.15) +
               annotate("text", x = 0.02 + box + 0.07, y = y, label = lv,
                        hjust = 0, vjust = 0.5, size = text_size)
      y <- y - step
    }
    y <- y - title_gap
  }
  p <- p +
    scale_x_continuous(limits = c(0, width_in), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, height_in), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", color = NA))
  list(plot = p, width = width_in, height = height_in)
}

tree_s1 <- upgma(bruvo_s1)
meta_s1 <- data.frame(Country  = df_fig1b[tree_s1$tip.label, "Pop_country"],
                      Genotype = df_fig1b[tree_s1$tip.label, "Pop_geno_map"],
                      Source   = provenance_category(
                        df_fig1b[tree_s1$tip.label,
                                 if ("Source" %in% names(df_fig1b)) "Source" else "source"]),
                      row.names = tree_s1$tip.label,
                      stringsAsFactors = FALSE)
fig1_palettes <- list(Country  = country_colors,
                      Genotype = map_geno_colors,
                      Source   = provenance_colors)
p_fig1_tree <- upgma_color_table_plot(tree_s1, meta_s1, fig1_palettes)
print(p_fig1_tree)
save_ggplot_multi(p_fig1_tree, "outputs/s1/Fig1_UPGMA_table",
                  width = 4, height = 10, dpi = 400, bg = "white")

fig1_legends <- upgma_color_legend_plot(meta_s1, fig1_palettes)
print(fig1_legends$plot)
save_ggplot_multi(fig1_legends$plot, "outputs/s1/Fig1_legends",
                  width = fig1_legends$width, height = fig1_legends$height,
                  dpi = 400, bg = "white")

# Statistical evaluation of group separation: flat AMOVA on the same
# Pop_geno_map grouping shown in the genotype-colored panel above (EU_2_A1 vs
# EC1 vs EU_13_A2 vs US7 vs US_1 vs PE_3 vs BR_1/EU_8_A1 vs Other).
amova_fig1b <- amova_analysis(
  gm          = s1$gm,
  df_in       = df_fig1b,
  sn          = s1$sn,
  formula_str = "~Pop_geno_map",
  strata_cols = c("Pop_geno_map"),
  prefix      = "S1_genotype_pcoa_merged",
  n_perm      = 999,
  outdir      = "outputs/s1"
)

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 8 — FIGURE: PCoA of EU_2_A1 (+PE_7.1) sub-population structure
#  Standalone block — does NOT depend on Section 5 / run_analysis(s2, ...)
#  having been run first.
#
#  Same PCoA (Bruvo distance, cmdscale) approach as Section 7, now restricted
#  to dataset 2 (EU_2_A1 + PE_7.1 — PE_7.1 already folded into EU_2_A1 via
#  Pop_geno, so s2 is effectively "all EU_2_A1"). Points are colored by
#  country (pie charts for MLGs spanning >1 country, same as Section 7).
#  Ellipses (linetype) group isolates into three a priori populations,
#  mirroring the DAPC-based structure discussed in Results:
#    - "Argentina + Uruguay" : the tight, low-diversity clonal group
#    - "Panama"              : the isolated, high-diversity group
#    - "Rest"                : all remaining EU_2_A1 countries (Bolivia,
#                              Peru, Brazil, Chile, etc.)
#  This grouping is local to this section (Pop_3pop, derived from s2$df) —
#  it does not touch s2$df itself or any other analysis.
# ══════════════════════════════════════════════════════════════════════════════
subpop_linetypes <- c(
  "Argentina + Uruguay" = "solid",
  "Panama"              = "dashed",
  "Rest"                = "dotted"
)

gc_base_s2 <- as.genclone(df2genind(s2$gm, sep = "/", ploidy = 3, ind.names = s2$sn))
bruvo_s2   <- bruvo.dist(gc_base_s2, replen = RL(gc_base_s2))

df_fig_eu2a1 <- s2$df
df_fig_eu2a1$Pop_3pop <- with(df_fig_eu2a1,
  ifelse(Pop_country %in% c("Argentina","Uruguay"), "Argentina + Uruguay",
  ifelse(Pop_country == "Panama", "Panama", "Rest"))
)

fig_eu2a1_pcoa <- pcoa_publication_plot(
  bruvo_dist            = bruvo_s2,
  df_in                 = df_fig_eu2a1,
  gc                    = gc_base_s2,
  point_group_col       = "Pop_country",
  point_colors_vec      = country_colors,
  point_legend_name     = "Country",
  point_style           = "pie",   # MLGs mixing >1 country render as pie charts
  point_alpha           = 0.70,
  pie_radius_range      = c(0.008, 0.025),   # ajustar a ojo según cómo se vea
  ellipse_group_col     = "Pop_3pop",
  ellipse_linetypes_vec = subpop_linetypes,
  filename              = "outputs/s2/Fig_EU2A1_PCoA_3subpops.png"
)

# Statistical evaluation: flat AMOVA on the three-population hypothesis
# (Argentina+Uruguay vs Panama vs Rest), matching exactly what the
# ellipses show above.
amova_eu2a1_subpop <- amova_analysis(
  gm          = s2$gm,
  df_in       = df_fig_eu2a1,
  sn          = s2$sn,
  formula_str = "~Pop_3pop",
  strata_cols = c("Pop_3pop"),
  prefix      = "S2_EU2A1_3subpop",
  n_perm      = 999,
  outdir      = "outputs/s2"
)

gc_s2_country <- as.genclone(df2genind(s2$gm, sep = "/", ploidy = 3,
                                       pop = s2$df$Pop_country, ind.names = s2$sn))
msn_s2 <- msn_plot(gc_s2_country, bruvo_s2, country_colors, seed = 12, nodelab = 5)
# 


set.seed(16)
plot_poppr_msn(gc_s2_country, msn_s2,
               palette = function(n) unname(country_colors[levels(pop(gc_s2_country))]),
               nodelab = 5, quantiles = FALSE,
               vertex.label.color = "transparent",
               gadj = 25, glim = c(0, 0.95))


pal_alpha <- function(n) {
  scales::alpha(unname(country_colors[levels(pop(gc_s2_country))]), alpha = 0.65)
}

set.seed(16)
plot_poppr_msn(gc_s2_country, msn_s2,
               palette = pal_alpha,
               nodelab = 5, quantiles = FALSE,
               vertex.label.color = "transparent",
               gadj = 25, glim = c(0, 0.95))

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 9 — PCoA restricted to the map genotypes (EU_2_A1, EC1, EU_13_A2,
#  US7, US_1) + PE_3 (PCoA-only, not on the map)
#  Standalone block — does NOT depend on any previous section having been run,
#  only on df1/df2 (Section 2) and pcoa_publication_plot() (Section 3).
#
#  Same genotype grouping used in the country-distribution map (Section 6,
#  make_eu2a1_map()'s cat_fn): EU_2_A1 (+ PE_7.1), EC1 (+ the EC_1 label as it
#  appears in df1/df2's raw Genotype column), EU_13_A2, US7 (+ the US_7 /
#  US_7_A2 label variants), and US_1 — all five are also on the country map.
#  PE_3 is added here too but is PCoA-only (not part of make_eu2a1_map()'s
#  category set / not on the map). Every other genotype is grouped as "Other"
#  and, in both plots below, rendered as a gray background point with no
#  population circle of its own — only the six named genotypes above get a
#  95%-CI ellipse, and that ellipse is always defined by genotype, never by
#  country (two isolates of the same genotype from different countries still
#  belong to the same population circle; a circle grouping isolates by
#  country instead would mix genotypes together, which is exactly what this
#  analysis is meant to avoid).
#
#  Country-label unification: the raw Country column (as read from
#  full_dataset_SSR.csv) is NOT yet pooled the way 00_prepare_filtered_dataset.R
#  does it for the narrower poppr_isolates_filtered.csv — it mixes real country
#  names with several different Euroblight reference-panel placeholder codes
#  that all mean the same non-field-collection origin: "ND", "GB", "UK", "ref",
#  "EN", "xx", "United Kingdom", "United Kingdom(Wales)", "Netherlands" and
#  "Belgium" are pooled here into a single "EU_controls" category (same label/
#  color already used everywhere else in this script). "USA" isolates (present
#  via the US7 lineage) get their own color below, since country_colors didn't
#  need one before this analysis. This pooling is scoped to a new column
#  (Pop_country_map) local to this section — it does not touch df1$Country,
#  s1$df$Pop_country, or any other analysis.
#
#  add_map_pcoa(s, dataset_label, file_prefix) works on either prepped subset
#  (s1 or s2, see Section 2) and produces both color variants below.
# ══════════════════════════════════════════════════════════════════════════════

eu_reference_codes <- c("ND", "GB", "UK", "ref", "EN", "xx",
                         "United Kingdom", "United Kingdom(Wales)",
                         "Netherlands", "Belgium")

# map_geno_colors / map_geno_linetypes and the categorization logic live in
# Section 1 (map_genotype_category()): the single shared genotype palette used
# across this whole script — see the comment there for what folds into what.

# country_colors extended with the one category this dataset needs that the
# base palette didn't have yet ("USA", present via the US7 lineage). Every
# unified country value that appears in the full dataset (Argentina, Bolivia,
# Brazil, Chile, Colombia, Costa Rica, EU_controls, Ecuador, Guatemala,
# Honduras, Panama, Peru, USA, Uruguay) already has a color in country_colors
# once "USA" is added, so every isolate in the country-colored PCoA gets its
# real country color — not just the six named map genotypes.
country_colors_map <- c(country_colors,
                         "USA" = "#b15928")

add_map_pcoa <- function(s, dataset_label, file_prefix, ellipses = TRUE) {
  d <- s$df

  # ── Unify country labels (see header comment above) ─────────────────────
  d$Country <- trimws(d$Country)
  d$Country[d$Country %in% eu_reference_codes] <- "EU_controls"

  # ── Map-genotype grouping (shared helper — see Section 1) ─────────────────
  d$Pop_geno_map <- map_genotype_category(d$Genotype)

  # Ellipses (population circles): only the six named genotypes above get one.
  # "Other" -> NA, so pcoa_publication_plot()'s table()-based group-size
  # filter drops it automatically (NA %in% valid_grps is FALSE, never TRUE).
  d$Pop_geno_map_ellipse <- ifelse(d$Pop_geno_map == "Other", NA, d$Pop_geno_map)

  # Country-coloring variant: every isolate gets its real (unified) country,
  # not just the six named map genotypes — the population circles (ellipses)
  # are still restricted to those six via Pop_geno_map_ellipse above, so
  # "who gets a circle" and "what color a point gets" stay independent.
  d$Pop_country_map <- d$Country

  cat(paste0("\n== Map-genotype PCoA -- ", dataset_label, " ==\n"))
  print(table(d$Pop_geno_map))
  cat("\nAll countries in this dataset (after unification):\n")
  print(table(d$Pop_country_map))
  cat("\nCountries among the named genotypes only:\n")
  print(table(d$Pop_country_map[d$Pop_geno_map != "Other"]))

  sn <- make.unique(d$Sample)
  rownames(d) <- sn
  gm <- d[, loci]
  rownames(gm) <- sn

  gc    <- as.genclone(df2genind(gm, sep = "/", ploidy = 3, ind.names = sn))
  bruvo <- bruvo.dist(gc, replen = RL(gc))

  # (1) Points colored by country (pie charts for MLGs spanning >1 country,
  # same convention as Figures 1b/EU2A1-subpop above); ellipses by genotype,
  # fixed color + linetype since points and ellipses use different variables.
  out_country <- pcoa_publication_plot(
    bruvo_dist            = bruvo,
    df_in                 = d,
    gc                    = gc,
    point_group_col       = "Pop_country_map",
    point_colors_vec      = country_colors_map,
    point_legend_name     = "Country",
    point_style           = "pie",
    ellipse_group_col     = if (ellipses) "Pop_geno_map_ellipse" else NULL,
    ellipse_linetypes_vec = map_geno_linetypes,
    ellipse_colors_vec    = NULL,
    filename              = paste0(file_prefix, "_PCoA_mapgenotypes_country.png")
  )

  # (2) Points AND ellipses colored by genotype (same palette -> ellipse line
  # always matches point fill); "Otros" points render gray, with no ellipse.
  out_genotype <- pcoa_publication_plot(
    bruvo_dist          = bruvo,
    df_in               = d,
    gc                  = gc,
    point_group_col     = "Pop_geno_map",
    point_colors_vec    = map_geno_colors,
    point_legend_name   = "Genotype",
    mixed_color         = "black",
    ellipse_group_col   = if (ellipses) "Pop_geno_map_ellipse" else NULL,
    ellipse_colors_vec  = map_geno_colors,
    filename            = paste0(file_prefix, "_PCoA_mapgenotypes_genotype.png")
  )

  invisible(list(df = d, country = out_country, genotype = out_genotype))
}

# Run on df1 (the full dataset, s1) — the only one of the two that actually
# contains all the named genotypes (df2/s2 is EU_2_A1(+PE_7.1) only, so EC1/
# EU_13_A2/US7/US_1/PE_3 would all fall into "Otros" there). add_map_pcoa(s2, ...)
# also works if that comparison is ever needed.
pcoa_map_df1 <- add_map_pcoa(s1, "df1 (full dataset)", "outputs/s1/Fig_mapgenotypes_df1",
                              ellipses = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 10 — LOCUS-LEVEL DIFFERENTIATION: which SSR locus/allele distinguishes
#  a given EU_2_A1 sub-population from the rest of the lineage?
#
#  Two questions, two matched tools:
#    (1) Panamanian EU_2_A1 isolates vs. the rest of EU_2_A1 (uses the Section 9
#        dataset, pcoa_map_df1$df). No prior discriminant analysis targets this
#        contrast, so locus_allele_freq_table() does a straightforward per-locus
#        allele-frequency comparison and flags any allele private to one group —
#        with only 12 SSR loci this is fully inspectable without a formal test.
#    (2) Argentina+Uruguay vs. the rest of EU_2_A1 within S2. Here a supervised
#        DAPC (Argentina+Uruguay vs. rest) was already fit and cross-validated
#        in dapc_cluster_validation() (val_s2, Section 5). dapc_locus_loadings()
#        refits a single DAPC at that already-validated optimal n.pca and plots
#        each allele's contribution to the discriminant axis (loadingplot()),
#        which ranks loci by discriminating power directly from that existing
#        analysis. locus_allele_freq_table() is also run on this contrast, as a
#        cross-check that the top DAPC-loading locus corresponds to a visibly
#        skewed or private allele and not an artifact of several weak loci.
# ══════════════════════════════════════════════════════════════════════════════

locus_allele_freq_table <- function(df_in, group_col, loci_vec = loci) {
  # For every locus, tabulate allele frequencies separately within each level
  # of group_col (the third-slot "0" padding used for diploid-looking calls is
  # dropped before tallying). Flags alleles private to a single group (present
  # there, absent from every other group) so a differential locus can be
  # spotted without needing a formal per-locus differentiation test.
  out <- list()
  for (l in loci_vec) {
    freqs <- lapply(split(df_in[[l]], df_in[[group_col]]), function(x) {
      al <- unlist(strsplit(x, "/"))
      al <- al[al != "0" & al != "" & !is.na(al)]
      if (length(al) == 0) return(NULL)
      round(table(al) / length(al), 3)
    })
    freqs <- freqs[!vapply(freqs, is.null, logical(1))]
    all_alleles <- unique(unlist(lapply(freqs, names)))
    private <- all_alleles[vapply(all_alleles, function(a) {
      present_in <- vapply(freqs, function(f) a %in% names(f), logical(1))
      sum(present_in) == 1
    }, logical(1))]
    out[[l]] <- list(freq = freqs, private_alleles = private)
  }
  out
}

print_locus_diff_report <- function(report) {
  for (l in names(report)) {
    cat("\n==", l, "==\n")
    print(report[[l]]$freq)
    if (length(report[[l]]$private_alleles) > 0) {
      cat("Private allele(s) at this locus:",
          paste(report[[l]]$private_alleles, collapse = ", "), "\n")
    }
  }
}

# --- (1) Panama vs. rest of EU_2_A1 (uses Section 9's pcoa_map_df1) ---------
d_panama            <- pcoa_map_df1$df
eu2a1_panama        <- d_panama[d_panama$Pop_geno_map == "EU_2_A1", ]
eu2a1_panama$grp    <- ifelse(eu2a1_panama$Country == "Panama",
                               "Panama", "Rest_EU2A1")

report_panama <- locus_allele_freq_table(eu2a1_panama, "grp")
cat("\n\n########## Panama vs. rest of EU_2_A1 ##########\n")
print_locus_diff_report(report_panama)

# --- (2a) Argentina+Uruguay vs. rest of EU_2_A1 — allele-frequency cross-check
s2_aruy         <- s2$df
s2_aruy$grp     <- ifelse(s2_aruy$Pop_country %in% c("Argentina", "Uruguay"),
                           "Arg_Uru", "Rest_EU2A1")

report_aruy <- locus_allele_freq_table(s2_aruy, "grp")
cat("\n\n########## Argentina+Uruguay vs. rest of EU_2_A1",
    "(allele-frequency check) ##########\n")
print_locus_diff_report(report_aruy)

# --- (2b) Argentina+Uruguay vs. rest of EU_2_A1 — DAPC loadings -------------
dapc_locus_loadings <- function(gen, pop_country, hyp_groups, n_pca, seed = 42) {
  hyp_label  <- paste(hyp_groups, collapse = " + ")
  pop_binary <- factor(ifelse(pop_country %in% hyp_groups, hyp_label, "Other"))
  pop(gen)   <- pop_binary

  set.seed(seed)
  dp <- dapc(gen, n.pca = n_pca, n.da = 1)
  print(loadingplot(dp$var.contr, axis = 1,
                     threshold = quantile(dp$var.contr, 0.9),
                     lab.jitter = 1,
                     main = paste0("Locus/allele contribution to ",
                                    hyp_label, " vs. rest discriminant axis")))
  invisible(dp)
}

n_pca_aruy <- val_s2$xval$`Number of PCs Achieving Highest Mean Success`
n_pca_aruy <- as.numeric(as.character(n_pca_aruy))[1]

dapc_aruy_loadings <- dapc_locus_loadings(
  gen         = gen_s2_base,
  pop_country = s2$df$Pop_country,
  hyp_groups  = c("Argentina", "Uruguay"),
  n_pca       = n_pca_aruy
)








# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 11 — DATA-DRIVEN POPULATION STRUCTURE OF EU_2_A1 (dataset S2)
#
#  Standalone block — depends only on s2 (Section 2) and the helpers of Section 3
#  (save_ggplot_multi, country_colors).
#
#  Goal: let the structure come out of the genotypes, without choosing groups.
#
#   (1) Unsupervised clustering: adegenet::find.clusters() = k-means on the principal
#       components of the allele frequencies (NOT a Bayesian method; the BIC is only
#       the score used to compare numbers of clusters K). It never sees the country.
#       Run on two versions of the data:
#         - all_isolates    : every isolate (identical clones are counted as many times as
#                             they were sampled, so the most sampled MLGs weigh heavily)
#         - clone_corrected : one isolate per MLG within each country
#       Reported: K chosen by "diffNgroup" over 30 random seeds, the BIC curve, how
#       reproducible the partition is across random starts (adjusted Rand index),
#       and the clusters x country tables for K = 2, 3, 4.
#   (2) Pairwise PhiST between countries (AMOVA on clone-corrected genotypes, Bruvo
#       distance, 9999 permutations, Holm-adjusted p) + hierarchical clustering of the
#       PhiST matrix: which countries group together, with country as the only
#       (sampling-unit) label and no pre-defined grouping of countries.
#
#  How to read it: the clusters are DESCRIBED by their composition. They must not be
#  used as groups in an AMOVA / cross-validation on the same data to "prove" the
#  structure: k-means maximizes between-cluster separation, so those tests are
#  significant by construction. (A cross-validation or AMOVA of a hypothesis stated
#  beforehand from geography is a different, legitimate analysis.)
#  Outputs: outputs/s2/unsupervised/
# ══════════════════════════════════════════════════════════════════════════════
suppressMessages(library(adegenet))
out11 <- "outputs/s2/unsupervised"
dir.create(out11, recursive = TRUE, showWarnings = FALSE)

# Adjusted Rand index between two partitions (1 = identical)
adjusted_rand <- function(a, b) {
  t  <- table(a, b); n <- sum(t); c2 <- function(x) x * (x - 1) / 2
  s  <- sum(c2(t)); ra <- sum(c2(rowSums(t))); cb <- sum(c2(colSums(t)))
  e  <- ra * cb / c2(n)
  (s - e) / ((ra + cb) / 2 - e)
}

# find.clusters() on a genind: K chosen by diffNgroup (k = NULL) or fixed K
fc_run <- function(g, seed, k = NULL, max_k = 8, n_pca = 50) {
  set.seed(seed)
  npc <- min(n_pca, nInd(g) - 1L)
  suppressMessages(
    if (is.null(k)) find.clusters(g, max.n.clust = max_k, n.pca = npc, choose.n.clust = FALSE,
                                  criterion = "diffNgroup", n.iter = 1e5, n.start = 10)
    else            find.clusters(g, n.clust = k, n.pca = npc, n.iter = 1e5, n.start = 10)
  )
}

# Full survey of one dataset: K by seed, BIC curve, stability, clusters x country
cluster_survey <- function(g, country, label, ks = 2:4, n_seeds = 30, n_starts = 10) {
  k_by_seed <- data.frame(dataset = label, seed = seq_len(n_seeds),
                          K = vapply(seq_len(n_seeds), function(s) nlevels(fc_run(g, s)$grp), 1L))
  bic <- fc_run(g, 42)$Kstat
  bic <- data.frame(dataset = label, K = seq_along(bic), BIC = as.numeric(bic))
  stab <- do.call(rbind, lapply(ks, function(k) {
    parts <- lapply(seq_len(n_starts), function(s) fc_run(g, s, k)$grp)
    ar <- combn(n_starts, 2, function(i) adjusted_rand(parts[[i[1]]], parts[[i[2]]]))
    data.frame(dataset = label, K = k, n_starts = n_starts,
               ARI_min = round(min(ar), 3), ARI_median = round(median(ar), 3))
  }))
  tabs <- lapply(ks, function(k) table(Cluster = fc_run(g, 42, k)$grp, Country = country))
  names(tabs) <- paste0("K", ks)
  list(k_by_seed = k_by_seed, bic = bic, stability = stab, tables = tabs)
}

# Pairwise PhiST (AMOVA, Bruvo distance) between all countries with >= min_n genotypes
pairwise_phist <- function(gc_cc, min_n = 5, n_perm = 9999, seed = 42) {
  cn   <- table(strata(gc_cc)$Country)
  keep <- names(cn[cn >= min_n])
  gc_cc <- gc_cc[strata(gc_cc)$Country %in% keep, ]
  # same loci for every pair: drop loci with > 5% missing data once, on the whole set
  gc_cc <- suppressMessages(missingno(gc_cc, type = "loci", cutoff = 0.05, quiet = TRUE))
  cat("PhiST loci kept:", paste(locNames(gc_cc), collapse = ", "), "\n")
  prs  <- combn(keep, 2, simplify = FALSE)
  res  <- do.call(rbind, lapply(prs, function(pr) {
    sub <- gc_cc[strata(gc_cc)$Country %in% pr, ]
    sub@strata$Country <- droplevels(factor(strata(sub)$Country))
    bd  <- suppressWarnings(bruvo.dist(sub, replen = RL(sub)))
    am  <- suppressWarnings(suppressMessages(
             poppr.amova(sub, ~Country, dist = bd, method = "ade4", missing = "ignore", quiet = TRUE)))
    set.seed(seed)
    rt  <- randtest(am, nrepet = n_perm)
    data.frame(country_1 = pr[1], country_2 = pr[2],
               n_1 = sum(strata(sub)$Country == pr[1]), n_2 = sum(strata(sub)$Country == pr[2]),
               PhiST = as.numeric(am$statphi$Phi), p = rt$pvalue)
  }))
  res$p_holm <- p.adjust(res$p, method = "holm")
  res
}

# ── Run ──────────────────────────────────────────────────────────────────────
cat("\n══ SECTION 11 — data-driven structure of EU_2_A1 ══\n")
pc_s2  <- s2$df$Pop_country
gen_u  <- df2genind(s2$gm, sep = "/", ploidy = 3, ind.names = s2$sn)
keep_u <- pc_s2 %in% names(table(pc_s2)[table(pc_s2) >= 3])        # as in dapc_cluster_validation()
mlg_u  <- mll(as.genclone(gen_u))
sets_u <- list(all_isolates    = which(keep_u),
               clone_corrected = which(keep_u & !duplicated(paste(mlg_u, pc_s2))))

surveys <- lapply(names(sets_u), function(nm) {
  idx <- sets_u[[nm]]
  cat("\n── ", nm, " (n = ", length(idx), ") ──\n", sep = "")
  sv <- cluster_survey(gen_u[idx, ], pc_s2[idx], nm)
  ks <- table(sv$k_by_seed$K)
  cat("K chosen by diffNgroup over 30 seeds: ",
      paste(names(ks), ks, sep = " -> ", collapse = " | "), "\n", sep = "")
  cat("BIC K=1..8: ", paste(round(sv$bic$BIC), collapse = " "), "\n", sep = "")
  print(sv$stability, row.names = FALSE)
  for (k in names(sv$tables)) { cat("\n", k, ": clusters x country\n", sep = ""); print(sv$tables[[k]]) }
  sv
})
names(surveys) <- names(sets_u)

write.csv(do.call(rbind, lapply(surveys, `[[`, "k_by_seed")), file.path(out11, "k_by_seed.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(surveys, `[[`, "bic")),       file.path(out11, "bic_by_k.csv"),   row.names = FALSE)
write.csv(do.call(rbind, lapply(surveys, `[[`, "stability")), file.path(out11, "cluster_stability_ari.csv"), row.names = FALSE)
k_summary <- do.call(rbind, lapply(surveys, function(sv) {
  t <- table(sv$k_by_seed$K); data.frame(dataset = sv$k_by_seed$dataset[1], K = as.integer(names(t)), n_seeds = as.integer(t))
}))
write.csv(k_summary, file.path(out11, "k_chosen_summary.csv"), row.names = FALSE)
long_tabs <- do.call(rbind, lapply(names(surveys), function(nm) do.call(rbind, lapply(names(surveys[[nm]]$tables), function(k) {
  d <- as.data.frame(surveys[[nm]]$tables[[k]]); data.frame(dataset = nm, K = k, d) }))))
write.csv(long_tabs, file.path(out11, "clusters_x_country_long.csv"), row.names = FALSE)
for (nm in names(surveys)) for (k in names(surveys[[nm]]$tables))
  write.csv(as.data.frame.matrix(surveys[[nm]]$tables[[k]]),
            file.path(out11, paste0("clusters_x_country_", nm, "_", k, ".csv")))

# ── Cross-validation of the a priori hypothesis Argentina+Uruguay vs the rest ──
# Same xvalDapc() set-up as dapc_cluster_validation() (100 replicates, 90/10 training/test
# split, up to 50 PCs). The groups here ARE defined by country (a priori), so this shows how
# separable they are, not that the structure exists. It is repeated on clone-corrected data
# (one isolate per MLG within each group) because, with all isolates, the isolates of the
# test set have identical clones in the training set.
xval_hypothesis <- function(gen, country, mlg, hyp = c("Argentina", "Uruguay"),
                            n_rep = 100, seed = 42) {
  hyp_label <- paste(hyp, collapse = " + ")
  grp <- factor(ifelse(country %in% hyp, hyp_label, "Other"))
  one <- function(idx, label) {
    x <- tab(gen[idx, ], freq = TRUE, NA.method = "mean")
    set.seed(seed)
    xv <- xvalDapc(x, grp[idx], n.pca.max = min(50, nrow(x) - 2L), training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE, n.rep = n_rep, xval.plot = FALSE)
    m  <- xv[["Mean Successful Assignment by Number of PCs of PCA"]]
    tb <- table(mlg[idx], grp[idx])
    data.frame(dataset = label, hypothesis = hyp_label,
               n_isolates = nrow(x), n_hypothesis_group = sum(grp[idx] == hyp_label),
               n_other = sum(grp[idx] == "Other"), n_MLGs = length(unique(mlg[idx])),
               n_MLGs_shared_between_groups = sum(tb[, 1] > 0 & tb[, 2] > 0),
               best_n_PCs = as.numeric(as.character(xv[["Number of PCs Achieving Highest Mean Success"]])),
               mean_correct_assignment_pct = round(100 * max(m), 1), n_replicates = n_rep)
  }
  rbind(one(seq_along(grp), "all_isolates"),
        one(which(!duplicated(paste(mlg, grp))), "clone_corrected (1 per MLG x group)"))
}
cat("\n── Cross-validation of Argentina+Uruguay vs the rest (all isolates and clone-corrected) ──\n")
xval_tab <- xval_hypothesis(gen_u, pc_s2, mlg_u)
write.csv(xval_tab, file.path(out11, "xval_hypothesis_ARUY.csv"), row.names = FALSE)
print(xval_tab, row.names = FALSE)

# BIC curves (K chosen most often marked)
bic_all <- do.call(rbind, lapply(surveys, `[[`, "bic"))
k_mode  <- do.call(rbind, lapply(surveys, function(sv) {
  t <- table(sv$k_by_seed$K); data.frame(dataset = sv$k_by_seed$dataset[1], K = as.integer(names(t)[which.max(t)]))
}))
p_bic <- ggplot(bic_all, aes(K, BIC)) +
  geom_line(color = "steelblue", linewidth = 0.8) + geom_point(color = "steelblue", size = 2.5) +
  geom_vline(data = k_mode, aes(xintercept = K), linetype = "dashed", color = "red") +
  facet_wrap(~dataset, scales = "free_y") + scale_x_continuous(breaks = 1:8) +
  labs(x = "Number of clusters (K)", y = "BIC",
       title = "k-means clustering of EU_2_A1: BIC by K (dashed: K most often chosen by diffNgroup)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
print(p_bic)
save_ggplot_multi(p_bic, file.path(out11, "bic_by_k"), width = 9, height = 3.8, dpi = 400)

# Pairwise PhiST between countries (clone-corrected)
cat("\n── Pairwise PhiST between countries (clone-corrected, 9999 permutations) ──\n")
gc_u <- as.genclone(gen_u)
strata(gc_u) <- data.frame(Country = pc_s2)
gc_u_cc <- clonecorrect(gc_u, strata = ~Country)
phist <- pairwise_phist(gc_u_cc, min_n = 5, n_perm = 9999)
write.csv(phist, file.path(out11, "pairwise_phiST_long.csv"), row.names = FALSE)
print(phist, row.names = FALSE)

cn <- sort(unique(c(phist$country_1, phist$country_2)))
Pm <- matrix(0, length(cn), length(cn), dimnames = list(cn, cn))
for (i in seq_len(nrow(phist))) { Pm[phist$country_1[i], phist$country_2[i]] <- phist$PhiST[i]
                                   Pm[phist$country_2[i], phist$country_1[i]] <- phist$PhiST[i] }
Pm[Pm < 0] <- 0                                           # a negative PhiST means "no differentiation"
write.csv(Pm, file.path(out11, "pairwise_phiST_matrix.csv"))
hc <- hclust(as.dist(Pm), method = "average")

hm <- expand.grid(A = cn, B = cn, stringsAsFactors = FALSE)
hm$PhiST <- Pm[cbind(hm$A, hm$B)]
hm$sig <- mapply(function(a, b) {
  if (a == b) return("")
  r <- phist[(phist$country_1 == a & phist$country_2 == b) | (phist$country_1 == b & phist$country_2 == a), ]
  if (r$p_holm < 0.05) "*" else "ns" }, hm$A, hm$B)
hm$A <- factor(hm$A, levels = cn[hc$order]); hm$B <- factor(hm$B, levels = rev(cn[hc$order]))
p_phi <- ggplot(hm[hm$A != hm$B, ], aes(A, B, fill = PhiST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(sprintf("%.2f", PhiST), sig)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, max(Pm)), name = expression(Phi[ST])) +
  labs(x = NULL, y = NULL, title = "EU_2_A1: pairwise PhiST between countries (clone-corrected)",
       subtitle = "* = significant after Holm correction (9999 permutations); ns = not significant\nCountries ordered by clustering. UK = European reference isolates") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
                                   plot.subtitle = element_text(size = 8))
print(p_phi)
save_ggplot_multi(p_phi, file.path(out11, "phiST_heatmap"), width = 7.5, height = 6.5, dpi = 400)

draw_hc <- function() { par(mar = c(2, 2, 3, 6)); plot(as.dendrogram(hc), horiz = TRUE,
  main = expression(paste("Countries clustered by pairwise ", Phi[ST], " (average linkage)"))) }
png(file.path(out11, "phiST_dendrogram.png"), width = 6, height = 4, units = "in", res = 300); draw_hc(); invisible(dev.off())
pdf(file.path(out11, "phiST_dendrogram.pdf"), width = 6, height = 4);                            draw_hc(); invisible(dev.off())
svg(file.path(out11, "phiST_dendrogram.svg"), width = 6, height = 4);                             draw_hc(); invisible(dev.off())
cat("\nSection 11 outputs written to ", out11, "\n", sep = "")

library(ggplot2)
library(sf)
library(rnaturalearth)
library(viridisLite)
library(ggspatial)
library(ggrepel)
library(cowplot)

isolates <- data.frame(
  Sample = c("Cba", "Tafi", "PiNSL", "Ar20079", "Ar21009", "Biocca", "Ar20112"),
  Lat    = c(-31.390527, -26.844469, -32.296175, -37.769626, -37.766904, -38.104092, -38.494389),
  Lon    = c(-63.866927, -65.752100, -65.995663, -58.306612, -58.295057, -57.670021, -58.655157)
)

# Response on cv. Innovator (Figure 4B): PiNSL-19, Ar21009, Biocca and Cba2
# trigger a hypersensitive response (HR = avirulent); Ar20079, Ar20112 and
# Tafi produce expanding lesions with sporulation (virulent).
isolates$Grupo <- factor(
  ifelse(isolates$Sample %in% c("PiNSL", "Biocca", "Cba", "Ar21009"),
         "Innovator_avirulent", "Innovator_virulent"),
  levels = c("Innovator_virulent", "Innovator_avirulent"),
  labels = c("Virulent on Innovator", "Avirulent on Innovator (HR)")
)

# Ar20079 (virulent) and Ar21009 (avirulent) come from the same site (their
# coordinates differ by ~0.003 deg), so at map scale one point hides the other.
# For drawing only, they are pulled apart east-west (Plot_Lon) so the two
# points sit side by side; Lat/Lon keep the real coordinates.
same_site_offset <- 0.25   # deg of longitude on each side
isolates$Plot_Lon <- isolates$Lon
isolates$Plot_Lon[isolates$Sample == "Ar20079"] <- isolates$Lon[isolates$Sample == "Ar20079"] - same_site_offset
isolates$Plot_Lon[isolates$Sample == "Ar21009"] <- isolates$Lon[isolates$Sample == "Ar21009"] + same_site_offset

grupo_colors <- setNames(
  c(viridis(1, begin = 0.25, end = 0.25), viridis(1, begin = 0.55, end = 0.55)),
  levels(isolates$Grupo)
)

arg <- ne_states(country = "Argentina", returnclass = "sf")

pad <- 1.2
bbox <- st_bbox(c(xmin = min(isolates$Lon) - pad, xmax = max(isolates$Lon) + pad,
                  ymin = min(isolates$Lat) - pad, ymax = max(isolates$Lat) + pad),
                crs = st_crs(4326))

# ── Main map: sin leyenda ni escala adentro ─────────────────────────────────
main_map <- ggplot() +
  geom_sf(data = arg, fill = "grey96", color = "grey55", linewidth = 0.25) +
  geom_point(data = isolates, aes(x = Plot_Lon, y = Lat, fill = Grupo),
             shape = 21, size = 4.2, color = "black", stroke = 0.6, alpha = 0.95) +
  geom_text_repel(data = isolates, aes(x = Plot_Lon, y = Lat, label = Sample),
                  size = 3.2, fontface = "italic", seed = 1,
                  segment.size = 0.3, min.segment.length = 0,
                  box.padding = 0.4, max.overlaps = Inf) +
  scale_fill_manual(values = grupo_colors, name = NULL) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1.0, "cm"), width = unit(1.0, "cm")) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_line(color = "grey90", linewidth = 0.2),
        legend.position = "none",
        panel.border = element_rect(color = "black", linewidth = 0.5))

# ── Leyenda, extraída como grob aparte ──────────────────────────────────────
legend_source <- main_map +
  geom_point(data = isolates, aes(x = Plot_Lon, y = Lat, fill = Grupo),
             shape = 21, size = 4.2, color = "black") +
  theme(legend.position = "right", legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = "white", color = "grey70"))
legend_grob <- get_legend(legend_source)

# ── Escala, calculada sobre el mismo bbox real (así la distancia es correcta)
# pero renderizada en un panel aparte para poder ubicarla afuera del mapa ────
scale_source <- ggplot() +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.6, text_cex = 0.9,
                   height = unit(0.3, "cm")) +
  theme_void()

# ── Inset locator ────────────────────────────────────────────────────────────
arg_country  <- ne_countries(country = "Argentina", returnclass = "sf")
locator_rect <- st_as_sfc(bbox)

inset_map <- ggplot() +
  geom_sf(data = arg_country, fill = "grey92", color = "grey50", linewidth = 0.2) +
  geom_sf(data = locator_rect, fill = NA, color = "red", linewidth = 0.6) +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white", color = NA))

# ── Composición: mapa a la izquierda; a su costado (sin superponerse al mapa)
# una columna con el inset arriba y la leyenda debajo ─────────────────────────
side_col <- plot_grid(inset_map, ggdraw() + draw_grob(legend_grob), NULL,
                      ncol = 1, rel_heights = c(0.45, 0.20, 0.35))
final_plot <- plot_grid(main_map, side_col, nrow = 1, rel_widths = c(0.64, 0.36))

final_plot

save_ggplot_multi(final_plot, "outputs/s2/mapa_argentina_innovator", width = 8, height = 7.8, dpi = 400)