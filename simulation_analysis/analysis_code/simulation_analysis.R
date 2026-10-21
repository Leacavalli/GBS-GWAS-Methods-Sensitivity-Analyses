################################################################################
# Simulation analysis: power and false-positive rate of GBS GWAS under 13
# population structure correction methods (pyseer on simulated phenotypes)
#
# Set 1 SNPs: range of allele frequencies (positive-control power simulations)
# Set 2 SNPs: range of allele frequencies x homoplasy indices (HI)
# Null sets: permuted phenotypes
#
# Inputs (in simulation_analysis/inputs/):
#   pyseer_power_set1_AF.txt, pyseer_power_set2_AFxHI.txt,
#   pyseer_null_set1_AF.txt,  pyseer_null_set2_AFxHI.txt,
#   homoplasy_index_per_SNP.tsv
#
# Generates (in simulation_analysis/outputs/):
#   Figure_3_Power_by_EOD_LOD_ratio_and_AF_set1.tiff
#   Figure_4_Power_by_HI_set2.tiff
#   Figure_5_Power_by_HI_AF_EOD_LOD_ratio_set2.tiff
#   Table_2_Power_low_vs_high_HI.csv
#   Figure_S1_AF_vs_HI.tiff
#   Figure_S2_Crude_OR_by_AF.tiff
#   Figure_S5_Power_by_AF_effect_size_set1.tiff
#   Figure_S6_Power_by_HI_set1.tiff
#   Figure_S7_AF_by_HI_set1.tiff
#   Figure_S8_Power_by_AF_effect_size_HI_set1.tiff
#   Figure_S9_Null_simulation_checks.tiff
#   Figure_S10_QQplots_null_simulations.tiff
#
# Run from RStudio (Source) or from the command line:
#   Rscript simulation_analysis/analysis_code/simulation_analysis.R
################################################################################

rm(list = ls())

# 0. SET UP ---------------------------------------------------------------------

## Paths (relative to this script: <repo>/simulation_analysis/analysis_code/)
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
}
repo_path <- dirname(dirname(get_script_dir()))
path      <- file.path(repo_path, "simulation_analysis", "inputs")
outpath   <- file.path(repo_path, "simulation_analysis", "outputs")
dir.create(outpath, showWarnings = FALSE)

## Packages
use_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in c("tidyverse", "patchwork", "scales")) use_pkg(pkg)

## Helper: save TIFF with LZW compression (lossless, much smaller files)
save_tiff <- function(plot, filename, width, height) {
  ggsave(file.path(outpath, filename), plot, width = width, height = height,
         dpi = 300, device = "tiff", compression = "lzw")
}

## Constants
N_EOD     <- 623
N_LOD     <- 917
SIG_LEVEL <- 2.55e-6   # Bonferroni threshold of the empirical SNP-GWAS

sim_cols <- c("correction", "variant", "af", "filter-pvalue", "lrt-pvalue", "beta", "beta-std-err")
cols_set1 <- c("AF_range", "MUT", "EOD_LOD_ratio", "simulation", sim_cols)
cols_set2 <- c("AF_range", "HI_range", "MUT", "EOD_LOD_ratio", "simulation", sim_cols)

method_levels <- c("CC", "SC", "Mash distances",
                   "SNP-based kinship matrix", "Unitig-based kinship matrix",
                   "COG-based kinship matrix", "RAxML-NG", "IQ-TREE",
                   "FastTree SP", "FastTree DP FEM", "FastTree DP LMM",
                   "FastTree DP LMM (Gubbins)", "No correction")
ratio_levels <- c("1:2", "1:3", "2:1", "3:1")

method_colors <- c(
  "CC"                          = "#1b9e77",
  "SC"                          = "#d95f02",
  "Mash distances"              = "#7570b3",
  "SNP-based kinship matrix"    = "#e7298a",
  "Unitig-based kinship matrix" = "#66a61e",
  "COG-based kinship matrix"    = "#e6ab02",
  "RAxML-NG"                    = "#a6761d",
  "IQ-TREE"                     = "#666666",
  "FastTree SP"                 = "#f4a582",
  "FastTree DP FEM"             = "#92c5de",
  "FastTree DP LMM"             = "#d6604d",
  "FastTree DP LMM (Gubbins)"   = "#4393c3",
  "No correction"               = "#878787"
)

## Helpers
recode_correction <- function(df) {
  df |>
    mutate(correction = recode(correction,
                               "Mash"                   = "Mash distances",
                               "snp_kinship"            = "SNP-based kinship matrix",
                               "unitig_kinship"         = "Unitig-based kinship matrix",
                               "COG_kinship"            = "COG-based kinship matrix",
                               "RAxML"                  = "RAxML-NG",
                               "IQtree"                 = "IQ-TREE",
                               "FasttreeSP"             = "FastTree SP",
                               "FasttreeDP_FEM"         = "FastTree DP FEM",
                               "FasttreeDP_LMM"         = "FastTree DP LMM",
                               "FasttreeDP_LMM_Gubbins" = "FastTree DP LMM (Gubbins)",
                               "nocorrection"           = "No correction"),
           correction = factor(correction, levels = method_levels))
}

add_sig <- function(df) {
  df |> mutate(sig = ifelse(`lrt-pvalue` < SIG_LEVEL, 1, 0),
               sig = ifelse(is.na(sig), 0, sig))
}

HI_bin_from_HI <- function(HI) {
  case_when(HI >= 0.90             ~ "[0.90, 1.00]",
            HI >= 0.80 & HI < 0.90 ~ "[0.80, 0.90)",
            HI >= 0.70 & HI < 0.80 ~ "[0.70, 0.80)",
            HI >= 0.60 & HI < 0.70 ~ "[0.60, 0.70)",
            HI >= 0.45 & HI < 0.60 ~ "[0.45, 0.60)",
            HI <  0.10             ~ "[0.00, 0.10)")
}

read_sim <- function(file, cols) {
  df <- read.delim(file.path(path, file), header = TRUE, colClasses = "character")
  colnames(df) <- cols   # internal column names used below
  as_tibble(df) |> mutate(af = as.numeric(af), `lrt-pvalue` = as.numeric(`lrt-pvalue`),
                          beta = as.numeric(beta))
}

power_heatmap <- function(df, x, y, fill, legend_name, title, subtitle, y_labels = waiver(),
                          label_fmt = function(v) round(v, 2), text_size = 4) {
  ggplot(df, aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label_fmt({{ fill }})), color = "black", size = text_size, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "darkred", name = legend_name,
                        limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_discrete(labels = y_labels) +
    labs(x = "Population Structure Correction Method", title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
          panel.grid      = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}


# 1. LOAD DATA ------------------------------------------------------------------

## Homoplasy index of every SNP in the empirical VCF
homoplasy_results <- read.delim(file.path(path, "homoplasy_index_per_SNP.tsv")) |>
  mutate(variant = paste("AP018935.1", POS, REF, ALT, sep = "_"))

## Allele frequencies from the empirical SNP-GWAS (FastTree DP LMM)
SNPGWAS_data <- read.table(
  file.path(repo_path, "empirical_analysis/inputs/pyseer/population_structure_methods/FastTree_DP_LMM/SNPGWAS.txt"),
  sep = "\t", header = TRUE)

## Set 1: power simulations across allele frequencies
sim_results <- read_sim("pyseer_power_set1_AF.txt", cols_set1) |>
  mutate(EOD_ratio = as.numeric(str_split_fixed(EOD_LOD_ratio, ":", 2)[, 1]),
         LOD_ratio = as.numeric(str_split_fixed(EOD_LOD_ratio, ":", 2)[, 2]),
         EOD_1 = round(af * (N_EOD + N_LOD) * (EOD_ratio / (EOD_ratio + LOD_ratio))),
         EOD_0 = N_EOD - EOD_1,
         LOD_1 = round(af * (N_EOD + N_LOD) * (LOD_ratio / (EOD_ratio + LOD_ratio))),
         LOD_0 = N_LOD - LOD_1,
         OR    = round((LOD_1 / (EOD_1 + LOD_1)) / (LOD_0 / (EOD_0 + LOD_0)), 1)) |>
  add_sig() |>
  recode_correction() |>
  mutate(EOD_LOD_ratio = factor(EOD_LOD_ratio, levels = ratio_levels)) |>
  left_join(homoplasy_results |> select(variant, HI), by = "variant") |>
  mutate(HI_bin = factor(HI_bin_from_HI(HI),
                         levels = c("[0.00, 0.10)", "[0.60, 0.70)", "[0.70, 0.80)",
                                    "[0.80, 0.90)", "[0.90, 1.00]")))

## Set 2: power simulations across homoplasy indices
sim_results_set2 <- read_sim("pyseer_power_set2_AFxHI.txt", cols_set2) |>
  left_join(homoplasy_results |> select(variant, HI), by = "variant") |>
  add_sig() |>
  recode_correction() |>
  mutate(EOD_LOD_ratio = factor(EOD_LOD_ratio, levels = ratio_levels),
         HI_bin = recode(HI_range,
                         "00_10" = "[0.00, 0.10)", "45_55" = "[0.45, 0.60)",
                         "60_70" = "[0.60, 0.70)", "70_80" = "[0.70, 0.80)",
                         "80_90" = "[0.80, 0.90)", "90_100" = "[0.90, 1.00]"),
         HI_bin = factor(HI_bin, levels = c("[0.00, 0.10)", "[0.45, 0.60)", "[0.60, 0.70)",
                                            "[0.70, 0.80)", "[0.80, 0.90)", "[0.90, 1.00]")))

## Null (permuted) simulations
sim_results_set1_null <- read_sim("pyseer_null_set1_AF.txt", cols_set1) |>
  add_sig() |>
  recode_correction()

sim_results_set2_null <- read_sim("pyseer_null_set2_AFxHI.txt", cols_set2) |>
  add_sig() |>
  recode_correction() |>
  left_join(homoplasy_results |> select(variant, HI), by = "variant") |>
  mutate(HI_bin = factor(HI_bin_from_HI(HI),
                         levels = c("[0.00, 0.10)", "[0.45, 0.60)", "[0.60, 0.70)",
                                    "[0.70, 0.80)", "[0.80, 0.90)", "[0.90, 1.00]")))


# 2. FIGURE S1: ALLELE FREQUENCY VS HOMOPLASY INDEX -----------------------------

HI_df <- SNPGWAS_data |>
  select(variant, af) |>
  left_join(homoplasy_results |> select(variant, HI), by = "variant") |>
  filter(!is.na(af), !is.na(HI))

base_theme <- theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"))

p_scatter <- ggplot(HI_df, aes(x = af, y = HI, color = HI)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.7) +
  scale_color_gradient2(low = "steelblue", mid = "grey80", high = "darkred",
                        midpoint = 0.5, name = "Homoplasy\nIndex (HI)") +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Allele Frequency (af)", y = "Homoplasy Index (HI)",
       title = "Allele Frequency vs Homoplasy Index") +
  base_theme

p_af <- ggplot(HI_df, aes(x = af)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(x = "Allele Frequency (af)", y = "Density", title = "Allele Frequency") +
  base_theme

p_hi <- ggplot(HI_df, aes(x = HI)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "darkorange", color = "white", alpha = 0.8) +
  labs(x = "Homoplasy Index (HI)", y = "Density", title = "Homoplasy Index") +
  base_theme

FigureS1 <- p_scatter / (p_af | p_hi) +
  plot_annotation(subtitle = paste0("n = ", nrow(HI_df), " SNPs"),
                  theme = theme(plot.subtitle = element_text(size = 10, color = "grey40")))
save_tiff(FigureS1, "Figure_S1_AF_vs_HI.tiff", width = 10, height = 10)


# 3. FIGURE S2: EXPECTED CRUDE OR BY AF AND EOD:LOD RATIO -----------------------

ratio_labels <- c("EOD:LOD = 1:3 (LOD strongly enriched)",
                  "EOD:LOD = 1:2 (LOD moderately enriched)",
                  "EOD:LOD = 2:1 (EOD moderately enriched)",
                  "EOD:LOD = 3:1 (EOD strongly enriched)")

OR_df <- expand.grid(AF = seq(0.05, 0.95, 0.05), EOD_LOD_ratio = c(1/4, 1/3, 2/3, 3/4)) |>
  mutate(N_1     = round(AF * (N_EOD + N_LOD)),
         N_0     = (N_EOD + N_LOD) - N_1,
         N_EOD_1 = EOD_LOD_ratio * N_1,
         N_LOD_1 = N_1 - N_EOD_1,
         N_EOD_0 = N_EOD - N_EOD_1,
         N_LOD_0 = N_LOD - N_LOD_1,
         OR      = (N_LOD_1 / N_EOD_1) / (N_LOD_0 / N_EOD_0),
         ratio_label = factor(ratio_labels[match(EOD_LOD_ratio, c(1/4, 1/3, 2/3, 3/4))],
                              levels = ratio_labels)) |>
  filter(is.finite(OR))

FigureS2 <- ggplot(OR_df, aes(x = AF, y = OR, colour = ratio_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), labels = percent) +
  scale_colour_manual(values = setNames(c("#d73027", "#fc8d59", "#4575b4", "#313695"), ratio_labels)) +
  labs(x = "Allele Frequency", y = "Crude Odds Ratio (log scale)",
       colour = "EOD:LOD ratio among carriers",
       title = "Expected crude OR by allele frequency and EOD:LOD ratio among carriers") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
save_tiff(FigureS2, "Figure_S2_Crude_OR_by_AF.tiff", width = 8, height = 5)


# 4. FIGURE 3: POWER BY EOD:LOD RATIO (TOP) AND ALLELE FREQUENCY (BOTTOM), SET 1 --

Figure3_top <- sim_results |>
  group_by(correction, EOD_LOD_ratio) |>
  summarise(power = signif(mean(sig == 1), 2), .groups = "drop") |>
  power_heatmap(correction, EOD_LOD_ratio, power, "Average\nproportion\nsignificant",
                "GWAS Detection Power",
                "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across allele frequencies 0.05–0.50 in 0.05 steps (20 simulations each)") +
  labs(y = "EOD:LOD Ratio")

Figure3_bottom <- sim_results |>
  group_by(correction, AF_range) |>
  summarise(power = signif(mean(sig == 1), 2), .groups = "drop") |>
  power_heatmap(correction, AF_range, power, "Average\nproportion\nsignificant",
                "GWAS Detection Power",
                "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across all EOD:LOD ratios") +
  labs(y = "Allele Frequency")

Figure3 <- Figure3_top / Figure3_bottom + plot_layout(heights = c(1, 2.5), guides = "collect")
save_tiff(Figure3, "Figure_3_Power_by_EOD_LOD_ratio_and_AF_set1.tiff", width = 8, height = 12)


# 5. FIGURE S5: POWER BY AF X EOD:LOD RATIO, NO CORRECTION VS FASTTREE DP LMM ---

af_ratio_heatmap <- function(correction_name, title, text_size) {
  sim_results |>
    filter(correction == correction_name) |>
    group_by(AF_range, EOD_LOD_ratio) |>
    summarise(power = mean(sig == 1), .groups = "drop") |>
    ggplot(aes(x = AF_range, y = EOD_LOD_ratio, fill = power)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(100 * power), "%")),
              color = "black", size = text_size, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "darkred", name = "Power", labels = percent,
                        breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
    labs(x = "Allele Frequency Bin", y = "EOD:LOD Ratio Among Carriers",
         title = title, subtitle = "% of simulations reaching significance") +
    theme_minimal(base_size = 13) +
    theme(panel.grid      = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1),
          plot.title      = element_text(face = "bold", size = 13),
          plot.subtitle   = element_text(size = 11, color = "grey40"),
          legend.position = "right")
}

FigureS5 <- af_ratio_heatmap("No correction", "No population structure correction", 4.5) +
  af_ratio_heatmap("FastTree DP LMM", "FastTree DP LMM correction", 5) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "GWAS detection power by allele frequency and effect size",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_tiff(FigureS5, "Figure_S5_Power_by_AF_effect_size_set1.tiff", width = 20, height = 7)


# 6. FIGURE S6: POWER BY HOMOPLASY INDEX, SET 1 ---------------------------------

hi_axis_labels <- function(df) {
  n <- df |> filter(!is.na(HI_bin)) |> group_by(HI_bin) |>
    summarise(n_snps = n_distinct(MUT), .groups = "drop")
  setNames(paste0(n$HI_bin, "\n(n=", n$n_snps, " SNPs)"), n$HI_bin)
}

hi_heatmap <- function(df, value, legend_name, title, subtitle) {
  df |>
    filter(!is.na(HI_bin)) |>
    group_by(correction, HI_bin) |>
    summarise(value = round(mean({{ value }} == 1), 2), .groups = "drop") |>
    power_heatmap(correction, HI_bin, value, legend_name, title, subtitle,
                  y_labels = hi_axis_labels(df), label_fmt = function(v) sprintf("%.2f", v),
                  text_size = 3.5) +
    labs(y = "Homoplasy Index (HI) Bin") +
    theme(axis.text.y = element_text(size = 9), axis.title.y = element_text(size = 11),
          legend.title = element_text(size = 10))
}

hi_subtitle <- "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across all allele frequencies and EOD:LOD ratios | 6 quantile-based HI bins"

FigureS6 <- hi_heatmap(sim_results, sig, "Average\nproportion\nsignificant",
                       "GWAS Detection Power by Homoplasy Level", hi_subtitle)
save_tiff(FigureS6, "Figure_S6_Power_by_HI_set1.tiff", width = 8, height = 5)


# 7. FIGURE S7: ALLELE FREQUENCY DISTRIBUTION BY HI BIN, SET 1 ------------------

set.seed(2024)   # reproducible jitter
n_per_hibin <- sim_results |> group_by(HI_bin) |> summarise(n_snps = n_distinct(MUT), .groups = "drop")

FigureS7 <- sim_results |>
  group_by(MUT, HI_bin) |>
  summarise(af_mean = mean(af, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = HI_bin, y = af_mean)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "grey30", fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5, color = "grey20") +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)),
               vjust = -0.8, color = "red", fontface = "bold", size = 3.5) +
  geom_text(data = n_per_hibin, aes(x = HI_bin, y = 0.75, label = paste0("n=", n_snps)),
            inherit.aes = FALSE, size = 3.5, color = "grey30") +
  labs(x = "Homoplasy Index Bin", y = "Allele Frequency",
       title = "Allele frequency distribution by homoplasy index bin",
       subtitle = "Median allele frequency shown in red; n = number of SNPs per bin") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x        = element_text(angle = 45, hjust = 1),
        plot.title         = element_text(face = "bold"),
        plot.subtitle      = element_text(color = "grey40", size = 11))
save_tiff(FigureS7, "Figure_S7_AF_by_HI_set1.tiff", width = 8, height = 5)


# 8. FIGURE S8: POWER BY AF X EOD:LOD RATIO X HI, NO CORRECTION VS FASTTREE DP LMM

make_hi_af_heatmap <- function(correction_name, title) {
  n_per_bin <- sim_results |>
    filter(!is.na(HI_bin)) |>
    group_by(HI_bin) |>
    summarise(n_snps = n_distinct(MUT), .groups = "drop") |>
    mutate(HI_bin_label = paste0(HI_bin, "\n(n=", n_snps, " SNPs)"))

  sim_results |>
    filter(correction == correction_name, !is.na(HI_bin)) |>
    left_join(n_per_bin, by = "HI_bin") |>
    mutate(HI_bin_label = factor(HI_bin_label, levels = rev(n_per_bin$HI_bin_label))) |>
    group_by(HI_bin_label, AF_range, EOD_LOD_ratio) |>
    summarise(power = round(mean(sig == 1), 2), .groups = "drop") |>
    mutate(label_color = ifelse(power > 0.6, "white", "black")) |>
    ggplot(aes(x = AF_range, y = EOD_LOD_ratio, fill = power)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(100 * power), "%"), color = label_color),
              size = 3.5, fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradient(low = "white", high = "darkred", name = "Power", labels = percent,
                        breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
    facet_grid(HI_bin_label ~ ., switch = "y") +
    labs(x = "Allele Frequency Bin", y = "EOD:LOD Ratio Among Carriers",
         title = title, subtitle = "% of simulations reaching significance") +
    theme_minimal(base_size = 13) +
    theme(panel.grid        = element_blank(),
          axis.text.x       = element_text(angle = 45, hjust = 1),
          plot.title        = element_text(face = "bold", size = 13),
          plot.subtitle     = element_text(size = 11, color = "grey40"),
          strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
          strip.placement   = "outside",
          legend.position   = "right")
}

FigureS8 <- make_hi_af_heatmap("No correction", "No population structure correction") +
  make_hi_af_heatmap("FastTree DP LMM", "FastTree DP LMM correction") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "GWAS detection power by allele frequency, effect size, and homoplasy index",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_tiff(FigureS8, "Figure_S8_Power_by_AF_effect_size_HI_set1.tiff", width = 20, height = 10)


# 9. FIGURE 4: POWER BY HOMOPLASY INDEX, SET 2 ----------------------------------

Figure4 <- hi_heatmap(sim_results_set2, sig, "Average\nproportion\nsignificant",
                      "GWAS Detection Power by Homoplasy Level", hi_subtitle)
save_tiff(Figure4, "Figure_4_Power_by_HI_set2.tiff", width = 8, height = 5)


# 10. TABLE 2: POWER FOR LOW- VS HIGH-HOMOPLASY VARIANTS, SET 2 -----------------

Table2 <- sim_results_set2 |>
  filter(HI_bin %in% c("[0.00, 0.10)", "[0.90, 1.00]")) |>
  group_by(correction, HI_bin) |>
  summarise(power = round(mean(sig == 1), 2), .groups = "drop") |>
  mutate(HI_label = ifelse(HI_bin == "[0.00, 0.10)", "Low HI (0.00–0.10)", "High HI (0.90–1.00)")) |>
  select(correction, HI_label, power) |>
  pivot_wider(names_from = HI_label, values_from = power) |>
  mutate(Diff_num = `High HI (0.90–1.00)` - `Low HI (0.00–0.10)`) |>
  arrange(Diff_num) |>
  mutate(`Low HI (0.00–0.10)`  = sprintf("%.2f", `Low HI (0.00–0.10)`),
         `High HI (0.90–1.00)` = sprintf("%.2f", `High HI (0.90–1.00)`),
         Difference = ifelse(Diff_num >= 0, paste0("+", sprintf("%.2f", Diff_num)), sprintf("%.2f", Diff_num))) |>
  select(Method = correction, `Low HI (0.00–0.10)`, `High HI (0.90–1.00)`, Difference)
print(Table2)
write.csv(Table2, file.path(outpath, "Table_2_Power_low_vs_high_HI.csv"), row.names = FALSE, fileEncoding = "UTF-8")


# 11. FIGURE 5: POWER BY HI ACROSS AF X EOD:LOD RATIO, SET 2 --------------------

eod_lod_labels <- c("1:2" = "1:2 (LOD moderately enriched)", "1:3" = "1:3 (LOD strongly enriched)",
                    "2:1" = "2:1 (EOD moderately enriched)", "3:1" = "3:1 (EOD strongly enriched)")
af_labels <- c("0.05" = "AF 0.05", "0.10" = "AF 0.10", "0.15" = "AF 0.15",
               "0.20" = "AF 0.20", "0.25" = "AF 0.25", "0.40" = "AF 0.40")

Figure5 <- sim_results_set2 |>
  group_by(correction, HI_bin, AF_range, EOD_LOD_ratio) |>
  summarise(power = mean(sig == 1), .groups = "drop") |>
  mutate(EOD_LOD_label = factor(eod_lod_labels[as.character(EOD_LOD_ratio)], levels = eod_lod_labels),
         AF_label      = factor(af_labels[AF_range], levels = rev(af_labels))) |>
  ggplot(aes(x = HI_bin, y = power, color = correction, group = correction)) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "red", linewidth = 0.8) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(AF_label ~ EOD_LOD_label) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(x = "Homoplasy Index (HI) Bin", y = "Average Proportion Significant",
       title = "GWAS Detection Power Across AF x EOD:LOD",
       subtitle = "Each panel = one EOD:LOD ratio  |  averaged across all allele frequencies  |  x-axis: low → high homoplasy\nRed dotted line = 80% power threshold") +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position  = "right",
        legend.text      = element_text(size = 9),
        legend.title     = element_text(size = 10, face = "bold"))
save_tiff(Figure5, "Figure_5_Power_by_HI_AF_EOD_LOD_ratio_set2.tiff", width = 18, height = 10)


# 12. FIGURE S9: NULL SIMULATION CHECKS -----------------------------------------

fpr_subtitle <- "Proportion of null (permuted) simulations with lrt-pvalue < 2.55×10⁻⁶"

## (A) Crude OR under the null (uncorrected beta as proxy for crude log-OR)
crude_OR_null_df <- sim_results_set1_null |>
  filter(correction == "No correction", is.finite(beta)) |>
  mutate(OR = exp(beta))
n_per_afbin_null <- crude_OR_null_df |> group_by(AF_range) |>
  summarise(n_snps = n_distinct(MUT), .groups = "drop")

FigureS9_A <- ggplot(crude_OR_null_df, aes(x = AF_range, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.alpha = 0.3, fill = "grey90", colour = "grey30") +
  geom_text(data = n_per_afbin_null,
            aes(x = AF_range, y = max(crude_OR_null_df$OR, na.rm = TRUE) * 1.1, label = paste0("n=", n_snps)),
            inherit.aes = FALSE, size = 3, color = "grey30") +
  scale_y_log10() +
  labs(x = "Allele Frequency Bin", y = "Crude Odds Ratio (log scale)",
       title = "Crude OR distribution under the null (permuted labels)",
       subtitle = "Uncorrected estimates across null simulations; dashed line = OR = 1 (no association)\nn = number of unique SNPs contributing to each bin") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.x   = element_text(angle = 45, hjust = 1))

## (B) False-positive rate by allele frequency (set 1, permuted)
FigureS9_B <- sim_results_set1_null |>
  group_by(correction, AF_range) |>
  summarise(fpr = signif(mean(sig == 1), 2), .groups = "drop") |>
  power_heatmap(correction, AF_range, fpr, "Empirical\nfalse positive\nrate",
                "GWAS False Positive Rate by Allele Frequency",
                paste0(fpr_subtitle, ",\nby allele frequency bin | Bonferroni threshold under the null")) +
  labs(y = "Allele Frequency")

## (C) False-positive rate by homoplasy index (set 2, permuted)
FigureS9_C <- hi_heatmap(sim_results_set2_null, sig, "Empirical\nfalse positive\nrate",
                         "GWAS False Positive Rate by Homoplasy Level",
                         paste0(fpr_subtitle, " | 6 quantile-based HI bins"))

FigureS9 <- FigureS9_A / (FigureS9_B + FigureS9_C) +
  plot_annotation(tag_levels = "A",
                  title = "Null simulation checks: calibration of GWAS correction methods",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_tiff(FigureS9, "Figure_S9_Null_simulation_checks.tiff", width = 16, height = 12)


# 13. FIGURE S10: P-VALUE CALIBRATION UNDER THE NULL (QQ PLOTS) -----------------

## Pool Set 1 + Set 2 null permutations
qq_null_df <- bind_rows(sim_results_set1_null |> select(correction, `lrt-pvalue`),
                        sim_results_set2_null |> select(correction, `lrt-pvalue`)) |>
  filter(!is.na(`lrt-pvalue`), `lrt-pvalue` > 0, `lrt-pvalue` <= 1)

## Genomic inflation factor per method
lambda_df <- qq_null_df |>
  group_by(correction) |>
  summarise(n = n(),
            lambda = round(median(qchisq(1 - `lrt-pvalue`, df = 1)) / qchisq(0.5, df = 1), 3),
            .groups = "drop")
print(lambda_df)

FigureS10 <- qq_null_df |>
  group_by(correction) |>
  arrange(`lrt-pvalue`, .by_group = TRUE) |>
  mutate(exp_neglog = -log10((row_number() - 0.5) / n()),
         obs_neglog = -log10(`lrt-pvalue`)) |>
  ungroup() |>
  left_join(lambda_df |> select(correction, lambda), by = "correction") |>
  mutate(facet_label = paste0(correction, "  (λ = ", sprintf("%.2f", lambda), ")")) |>
  ggplot(aes(x = exp_neglog, y = obs_neglog)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_point(alpha = 0.4, size = 0.8, colour = "steelblue") +
  facet_wrap(~ facet_label, scales = "free") +
  labs(x = expression(Expected ~ -log[10](italic(p)) ~ "(uniform null)"),
       y = expression(Observed ~ -log[10](italic(p)) ~ "(null simulations)"),
       title = "P-value calibration under the null (QQ-plots)",
       subtitle = "Combined Set 1 + Set 2 null permutations | red dashed line = perfect calibration | λ = genomic inflation factor") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        strip.text    = element_text(size = 8, face = "bold"))
save_tiff(FigureS10, "Figure_S10_QQplots_null_simulations.tiff", width = 16, height = 12)

sessionInfo()
