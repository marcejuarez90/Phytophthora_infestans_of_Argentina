# ══════════════════════════════════════════════════════════════════════════════
#  make_figsup_sampling_years.R
#
#  Supplementary Figure S1: isolates per sampling year in the SSR isolate table
#  (inputs/full_dataset_SSR.csv).
#    A (top left)      all isolates, colored by genotype
#    B (top right)     all isolates, colored by country
#    C (bottom left)   EU_2_A1 (+ PE_7.1), colored by genotype
#    D (bottom right)  EU_2_A1 (+ PE_7.1), colored by country
#  Genotypes use map_genotype_category() / map_geno_colors, countries use
#  country_colors (Euroblight reference codes pooled as EU_controls, as in
#  Section 9 of the main script). Isolates without a sampling year are excluded
#  (their number is printed). Page 10 x 8.4 in, four equal panels.
#
#  Palettes and the genotype categorization are read from Section 1 of
#  haplotype_network_v4.R, so they stay identical to every other figure.
#
#  Run from this folder:   Rscript make_figsup_sampling_years.R
#  Output: one PDF and one PNG per panel in outputs/s1/.
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
library(grid)

# ── Section 1 of the main script (palettes + map_genotype_category()) ─────────
main_src <- readLines("haplotype_network_v4.R", encoding = "UTF-8", warn = FALSE)
s1_from  <- grep("SECTION 1 — GLOBAL PARAMETERS", main_src)
s1_to    <- grep("SECTION 2 — DATA LOADING", main_src)
stopifnot(length(s1_from) == 1, length(s1_to) == 1)
eval(parse(text = main_src[s1_from:(s1_to - 2)], encoding = "UTF-8"))

# Same pooling / extra color as Section 9 (eu_reference_codes, country_colors_map).
eu_reference_codes <- c("ND", "GB", "UK", "ref", "EN", "xx",
                        "United Kingdom", "United Kingdom(Wales)",
                        "Netherlands", "Belgium")
country_colors_map <- c(country_colors, "USA" = "#b15928")

# ── Data ───────────────────────────────────────────────────────────────────────
d <- read.csv("inputs/full_dataset_SSR.csv", stringsAsFactors = FALSE)
d$Year         <- as.integer(d$Year)
d$Pop_geno_map <- map_genotype_category(d$Genotype)
d$Country      <- trimws(d$Country)
d$Country[d$Country %in% eu_reference_codes] <- "EU_controls"

e <- d[d$Genotype %in% c("EU_2_A1", "PE_7.1"), ]

cat("All isolates:", nrow(d), "| without year (excluded):", sum(is.na(d$Year)), "\n")
cat("EU_2_A1     :", nrow(e), "| without year (excluded):", sum(is.na(e$Year)), "\n")
print(e[is.na(e$Year), c("Sample", "Genotype", "Country", "source")], row.names = FALSE)

d <- d[!is.na(d$Year), ]
e <- e[!is.na(e$Year), ]
year_levels <- as.character(seq(min(d$Year), max(d$Year)))
d$Year_lab  <- factor(d$Year, levels = year_levels)
e$Year_lab  <- factor(e$Year, levels = year_levels)

# ── Bar panels ─────────────────────────────────────────────────────────────────
PANEL_W <- 5; PANEL_H <- 4

year_bar_plot <- function(df, fill_col, palette, legend_title) {
  df$fill <- factor(df[[fill_col]], levels = intersect(names(palette), unique(df[[fill_col]])))
  every5  <- year_levels[as.integer(year_levels) %% 5 == 0]
  ggplot(df, aes(x = Year_lab, fill = fill)) +
    geom_bar(width = 0.8, colour = NA) +
    scale_x_discrete(drop = FALSE, breaks = every5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = palette, name = legend_title, drop = TRUE) +
    labs(x = "Sampling year", y = "Number of isolates",
         subtitle = sprintf("n = %d", nrow(df))) +
    theme_classic(base_size = 10) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.key.size = unit(0.32, "cm"), legend.text = element_text(size = 7.5),
          legend.title = element_text(size = 8.5, face = "bold"),
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 8.5, hjust = 1),
          plot.margin = margin(22, 8, 4, 8))   # top margin leaves room for the panel letter
}

panels <- list(
  A = year_bar_plot(d, "Pop_geno_map", map_geno_colors,    "Genotype"),
  B = year_bar_plot(d, "Country",      country_colors_map, "Country"),
  C = year_bar_plot(e, "Pop_geno_map", map_geno_colors,    "Genotype"),
  D = year_bar_plot(e, "Country",      country_colors_map, "Country")
)
# Lower-left corner of each panel (NPC of the page): top row = all isolates,
# bottom row = EU_2_A1.
panel_xy <- list(A = c(0, 0.5), B = c(0.5, 0.5), C = c(0, 0), D = c(0.5, 0))

# ── Each panel on its own ──────────────────────────────────────────────────────
# One file per panel: how they are arranged on the page is a presentation choice
# and is not fixed here.
dir.create("outputs/s1", showWarnings = FALSE, recursive = TRUE)
panel_file <- c(A = "sampling_years_all_by_genotype",
                B = "sampling_years_all_by_country",
                C = "sampling_years_EU2A1_by_genotype",
                D = "sampling_years_EU2A1_by_country")
for (k in names(panels)) {
  out <- file.path("outputs/s1", panel_file[[k]])
  ggsave(paste0(out, ".pdf"), panels[[k]], width = PANEL_W, height = PANEL_H)
  ggsave(paste0(out, ".png"), panels[[k]], width = PANEL_W, height = PANEL_H, dpi = 300)
  cat("Wrote", paste0(out, c(".pdf", ".png"), collapse = " , "), "\n")
}
