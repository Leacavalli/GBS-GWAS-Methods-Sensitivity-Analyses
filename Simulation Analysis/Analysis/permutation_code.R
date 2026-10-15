# Clean environment
rm( list = ls() )

# Load libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(tidyr)
library(tidyverse)

# set local path
path <- "~/GitHub/GBS-GWAS-Methods-Sensitivity-Analyses/"

# ── Data loading and processing ──────────────────────────────────
sim_results_set1_permuted <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/sim_GWAS_results_set1_permuted.txt"),
  header = F)[-1,]
colnames(sim_results_set1_permuted)<- c(
  "AF_range", "MUT", "EOD_LOD_ratio", "simulation",
  "correction", "variant", "af", "filter-pvalue",
  "lrt-pvalue", "beta", "beta-std-err"
)

sim_results_set2_permuted <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/sim_GWAS_results_set2_permuted.txt"),
  header =  F)[-1,] 
colnames(sim_results_set2_permuted) <- c("AF_range","HI_range",  "MUT", "EOD_LOD_ratio", "simulation",
                                    "correction", "variant", "af", "filter-pvalue",
                                    "lrt-pvalue", "beta", "beta-std-err")

# ── Homoplasy lookup (same source used for the power analysis) ──
# NOTE: if this object isn't already in your environment from the power-
# analysis script, uncomment the lines below to rebuild it.
homoplasy_results <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/homoplasy_results_SNPs.tsv")
) |>
  mutate(variant = paste("AP018935.1", POS, REF, ALT, sep = "_"))

# ── Method factor levels (same order used throughout the power analysis) ──
method_levels <- c("CC", "SC", "Mash distances",
                   "SNP-based kinship matrix", "Unitig-based kinship matrix",
                   "COG-based kinship matrix", "RAxML-NG", "IQ-TREE",
                   "FastTree SP", "FastTree DP FEM", "FastTree DP LMM",
                   "FastTree DP LMM (Gubbins)", "No correction")

recode_correction <- function(df) {
  df |>
    mutate(
      correction = recode(correction,
                          "CC"                     = "CC",
                          "SC"                     = "SC",
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
                          "nocorrection"           = "No correction"
      ),
      correction = factor(correction, levels = method_levels)
    )
}

# ── Shared heatmap theme (matches power-analysis figures) ───────
heatmap_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  )

# ── Colour palette per method (reused from the power-analysis script) ──
method_colors <- c(
  "CC"                        = "#1b9e77",
  "SC"                        = "#d95f02",
  "Mash distances"            = "#7570b3",
  "SNP-based kinship matrix"  = "#e7298a",
  "Unitig-based kinship matrix" = "#66a61e",
  "COG-based kinship matrix"  = "#e6ab02",
  "RAxML-NG"                  = "#a6761d",
  "IQ-TREE"                   = "#666666",
  "FastTree SP"               = "#f4a582",
  "FastTree DP FEM"           = "#92c5de",
  "FastTree DP LMM"           = "#d6604d",
  "FastTree DP LMM (Gubbins)" = "#4393c3",
  "No correction"             = "#878787"
)



#####################################################
#  Set 1 (permuted) — FPR by allele frequency        #
#  Null equivalent of Figure 3 (bottom panel)/Fig. 4  #
#####################################################

sim_results_set1_null <- as_tibble(sim_results_set1_permuted) |>
  filter(AF_range != "freq") |>
  mutate(
    af  = as.numeric(af),
    sig = ifelse(`lrt-pvalue` < 2.55e-6, 1, 0),
    sig = ifelse(is.na(sig), 0, sig)
  ) |>
  recode_correction() |>
  left_join(homoplasy_results |> select(variant, HI), by = "variant")

# ── FPR summary: method × allele frequency bin ───────────────────
fpr_af_df <- sim_results_set1_null |>
  group_by(correction, AF_range) |>
  summarise(fpr = signif(mean(sig == 1), 2), .groups = "drop")

FigureS_FPR_AF <- fpr_af_df |>
  ggplot(aes(x = correction, y = AF_range, fill = fpr)) +
  geom_tile(color = "white", linewidth = 0.5)  +
  geom_text(aes(label = sprintf("%.1f", fpr)),
            color = "black", size = 3.5, fontface = "bold")  +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Empirical\nfalse positive\nrate",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Allele Frequency",
    title    = "GWAS False Positive Rate\nby Allele Frequency",
    subtitle = "Proportion of null (permuted) simulations with lrt-pvalue < 2.55×10⁻⁶,\nby allele frequency bin | Bonferroni threshold under the null"
  ) +
  heatmap_theme

FigureS_FPR_AF

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_FPR_AF.tiff"), plot = FigureS_FPR_AF, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_FPR_AF.png"),  plot = FigureS_FPR_AF, width = 8, height = 5, dpi = 300)



#####################################################
#  Diagnostic: crude OR under the null, by AF bin    #
#  (companion to Figure 2 / Figure S2)               #
#####################################################

# Uses the uncorrected ("No correction") beta as a proxy for crude log-OR.
# Under the null (permuted labels), this should be centred on OR = 1 across
# allele frequency bins, unlike the theoretical bias curves in Figure 2/S2
# (which apply to the artificially enriched EOD:LOD ratios used in the
# power simulations, not the true 0.69:1 ratio used here).

crude_OR_null_df <- sim_results_set1_null |>
  filter(correction == "No correction", is.finite(beta)) |>
  mutate(OR = exp(beta))

# ── n unique SNPs contributing per AF bin ─────────────────────────
n_per_afbin_null <- crude_OR_null_df |>
  group_by(AF_range) |>
  summarise(n_snps = n_distinct(MUT), .groups = "drop")

FigureS_CrudeOR_Null <- crude_OR_null_df |>
  ggplot(aes(x = AF_range, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.alpha = 0.3, fill = "grey90", colour = "grey30") +
  geom_text(
    data = n_per_afbin_null,
    aes(x = AF_range, y = max(crude_OR_null_df$OR, na.rm = TRUE) * 1.1,
        label = paste0("n=", n_snps)),
    inherit.aes = FALSE, size = 3, color = "grey30"
  ) +
  scale_y_log10() +
  labs(
    x     = "Allele Frequency Bin",
    y     = "Crude Odds Ratio (log scale)",
    title = "Crude OR distribution under the null (permuted labels)",
    subtitle = "Uncorrected estimates across null simulations; dashed line = OR = 1 (no association)\nn = number of unique SNPs contributing to each bin"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )

FigureS_CrudeOR_Null

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_CrudeOR_Null.tiff"), plot = FigureS_CrudeOR_Null, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_CrudeOR_Null.png"),  plot = FigureS_CrudeOR_Null, width = 8, height = 5, dpi = 300)



#####################################################
#  Set 2 (permuted) — FPR by homoplasy index         #
#  Null equivalent of Figure 4 / Figure 9            #
#####################################################

sim_results_set2_null <- as_tibble(sim_results_set2_permuted|>
                                     mutate(str_remove(correction, "ratio_null/"))) |>
  filter(AF_range != "freq") |>
  mutate(
    af  = as.numeric(af),
    sig = ifelse(`lrt-pvalue` < 2.55e-6, 1, 0),
    sig = ifelse(is.na(sig), 0, sig)
  ) |>
  recode_correction() |>
  left_join(homoplasy_results |> select(variant, HI), by = "variant") |>
  mutate(HI_bin = case_when(
    HI >= 0.90              ~ "[0.90, 1.00]",
    HI >= 0.80 & HI < 0.90  ~ "[0.80, 0.90)",
    HI >= 0.70 & HI < 0.80  ~ "[0.70, 0.80)",
    HI >= 0.60 & HI < 0.70  ~ "[0.60, 0.70)",
    HI >= 0.45 & HI < 0.60  ~ "[0.45, 0.60)",
    HI <  0.10              ~ "[0.00, 0.10)"
  )) |>
  mutate(HI_bin = factor(HI_bin, levels = c(
    "[0.00, 0.10)",
    "[0.45, 0.60)",
    "[0.60, 0.70)",
    "[0.70, 0.80)",
    "[0.80, 0.90)",
    "[0.90, 1.00]"
  )))

# ── n SNPs per HI bin for axis labels ─────────────────────────────
n_per_hibin_null <- sim_results_set2_null |>
  filter(!is.na(HI_bin)) |>
  group_by(HI_bin) |>
  summarise(n_snps = n_distinct(MUT), .groups = "drop")

hi_labels_null <- setNames(
  paste0(levels(sim_results_set2_null$HI_bin), "\n(n=", n_per_hibin_null$n_snps, " SNPs)"),
  levels(sim_results_set2_null$HI_bin)
)

# ── FPR summary: method × HI bin ─────────────────────────────────
fpr_HI_df <- sim_results_set2_null |>
  filter(!is.na(HI_bin)) |>
  group_by(correction, HI_bin) |>
  summarise(fpr = round(mean(sig == 1), 2), .groups = "drop")

FigureS_FPR_HI <- fpr_HI_df |>
  ggplot(aes(x = correction, y = HI_bin, fill = fpr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", fpr)),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Empirical\nfalse positive\nrate",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = hi_labels_null) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Homoplasy Index (HI) Bin",
    title    = "GWAS False Positive Rate\nby Homoplasy Level",
    subtitle = "Proportion of null (permuted) simulations with lrt-pvalue < 2.55×10⁻⁶ | 6 quantile-based HI bins"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 9),
    axis.title.y    = element_text(size = 11),
    legend.position = "right",
    legend.title    = element_text(size = 10)
  )

FigureS_FPR_HI

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_FPR_HI.tiff"), plot = FigureS_FPR_HI, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_FPR_HI.png"),  plot = FigureS_FPR_HI, width = 8, height = 5, dpi = 300)



#####################################################
#  Table: FPR at low vs. high HI per method          #
#  Null equivalent of Table 2                        #
#####################################################


TableS_FPR_lowhighHI <- fpr_HI_df |>
  filter(HI_bin %in% c("[0.00, 0.10)", "[0.90, 1.00]")) |>
  mutate(HI_label = ifelse(HI_bin == "[0.00, 0.10)",
                           "Low HI (0.00–0.10)",
                           "High HI (0.90–1.00)")) |>
  select(correction, HI_label, fpr) |>
  pivot_wider(names_from = HI_label, values_from = fpr) |>
  mutate(
    Diff_num             = `High HI (0.90–1.00)` - `Low HI (0.00–0.10)`,
    `Low HI (0.00–0.10)`  = sprintf("%.2f", `Low HI (0.00–0.10)`),
    `High HI (0.90–1.00)` = sprintf("%.2f", `High HI (0.90–1.00)`)
  ) |>
  arrange(Diff_num) |>
  mutate(Difference = ifelse(Diff_num >= 0,
                             paste0("+", sprintf("%.2f", Diff_num)),
                             sprintf("%.2f", Diff_num))) |>
  select(-Diff_num) |>
  rename(Method = correction)

TableS_FPR_lowhighHI



#####################################################
#  Line plot: FPR vs. HI bin, one line per method    #
#  Null equivalent of Figure 5 / Figure 10           #
#  (no EOD:LOD facet — the null model uses only the  #
#   true observed 0.69:1 ratio, see Figure S2)       #
#####################################################

fpr_HI_line_df <- sim_results_set2_null |>
  filter(!is.na(HI_bin)) |>
  group_by(correction, HI_bin) |>
  summarise(fpr = mean(sig == 1), .groups = "drop")

FigureS_FPR_HI_line <- fpr_HI_line_df |>
  ggplot(aes(x = HI_bin, y = fpr, color = correction, group = correction)) +
  geom_hline(yintercept = 2.55e-6 * 100, linetype = "dotted",
             color = "red", linewidth = 0.5) +  # nominal alpha reference (~0%)
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x        = "Homoplasy Index (HI) Bin",
    y        = "Empirical False Positive Rate",
    title    = "GWAS False Positive Rate Across Homoplasy Index",
    subtitle = "Null model uses the single true observed EOD:LOD ratio (0.69:1) — no ratio facet"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold")
  )

FigureS_FPR_HI_line

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_FPR_HI_line.tiff"), plot = FigureS_FPR_HI_line, width = 10, height = 6, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_FPR_HI_line.png"),  plot = FigureS_FPR_HI_line, width = 10, height = 6, dpi = 300)


FigureS_NullChecks <- (FigureS_CrudeOR_Null) /
  (FigureS_FPR_AF+ FigureS_FPR_HI) +
  plot_annotation(
    tag_levels = "A",
    title = "Null simulation checks: calibration of GWAS correction methods",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

FigureS_NullChecks




#####################################################
#  P-value calibration: QQ-plots under the null      #
#  (full-distribution check, not just threshold FPR) #
#####################################################

# Combine both null simulation sets (Set 1 + Set 2) so every SNP/permutation
# contributes to each method's calibration check.
qq_null_df <- bind_rows(
  sim_results_set1_null |> select(correction, `lrt-pvalue`),
  sim_results_set2_null |> select(correction, `lrt-pvalue`)
) |>
  filter(!is.na(`lrt-pvalue`), `lrt-pvalue` > 0, `lrt-pvalue` <= 1)

# ── Build expected-vs-observed quantiles per method ──────────────
qq_plot_df <- qq_null_df |>
  group_by(correction) |>
  arrange(`lrt-pvalue`, .by_group = TRUE) |>
  mutate(
    n           = n(),
    rank        = row_number(),
    expected_p  = (rank - 0.5) / n,
    obs_neglog  = -log10(`lrt-pvalue`),
    exp_neglog  = -log10(expected_p)
  ) |>
  ungroup()

# ── Genomic inflation factor (lambda) per method ─────────────────
# lambda = median(observed chi-sq) / median(expected chi-sq under null)
# lambda ~ 1  -> well calibrated
# lambda > 1  -> anti-conservative (too many small p-values; inflated FPR)
# lambda < 1  -> conservative (too few small p-values; underpowered test)
lambda_df <- qq_null_df |>
  group_by(correction) |>
  summarise(
    n      = n(),
    lambda = round(
      median(qchisq(1 - `lrt-pvalue`, df = 1), na.rm = TRUE) /
        qchisq(0.5, df = 1),
      3
    ),
    .groups = "drop"
  )

lambda_df

# ── QQ-plot, faceted by correction method ─────────────────────────
FigureS_QQ_Null <- qq_plot_df |>
  left_join(lambda_df |> select(correction, lambda), by = "correction") |>
  mutate(facet_label = paste0(correction, "  (λ = ", sprintf("%.2f", lambda), ")")) |>
  ggplot(aes(x = exp_neglog, y = obs_neglog)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_point(alpha = 0.4, size = 0.8, colour = "steelblue") +
  facet_wrap(~ facet_label, scales = "free") +
  labs(
    x     = expression(Expected~-log[10](italic(p))~"(uniform null)"),
    y     = expression(Observed~-log[10](italic(p))~"(null simulations)"),
    title = "P-value calibration under the null (QQ-plots)",
    subtitle = "Combined Set 1 + Set 2 null permutations | red dashed line = perfect calibration | λ = genomic inflation factor"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text    = element_text(size = 8, face = "bold")
  )

FigureS_QQ_Null

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_QQ_Null.tiff"), plot = FigureS_QQ_Null, width = 16, height = 12, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/FigureS_QQ_Null.png"),  plot = FigureS_QQ_Null, width = 16, height = 12, dpi = 300)
