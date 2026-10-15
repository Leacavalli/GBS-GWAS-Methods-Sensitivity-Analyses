
# Clean environment
rm( list = ls() )

# Load libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

# set local path 
path <- "~/GitHub/GBS-GWAS-Methods-Sensitivity-Analyses/"



#########################
#       Figure 1        #
#########################

# ── Data loading ─────────────────────────────────────────────────
homoplasy_results <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/homoplasy_results_SNPs.tsv")
) |>
  mutate(variant = paste("AP018935.1", POS, REF, ALT, sep = "_"))

SNPGWAS_data <- read.table(
  paste0(path, "Empirical Analysis/inputs/pyseer/Population_Structure_Methods/FastTree_DP_LMM/SNPGWAS.txt"),
  sep = "\t", header = TRUE
) |>
  mutate(
    POS    = as.integer(str_split(as.character(variant), "_", simplify = TRUE)[, 2]),
    OR     = exp(beta),
    OR_pos = ifelse(OR < 1, "OR<1", "OR>1")
  )

# ── Join AF and HI ───────────────────────────────────────────────
HI_df <- left_join(
  SNPGWAS_data     |> select(variant, af),
  homoplasy_results |> select(variant, HI),
  by = "variant"
) |>
  filter(!is.na(af), !is.na(HI))

n_snps <- nrow(HI_df)

# ── Shared theme ─────────────────────────────────────────────────
base_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ── Panel A: AF vs HI scatter ────────────────────────────────────
p_scatter <- HI_df |>
  ggplot(aes(x = af, y = HI, color = HI)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black",
              linetype = "dashed", linewidth = 0.7) +
  scale_color_gradient2(low = "steelblue", mid = "grey80", high = "darkred",
                        midpoint = 0.5, name = "Homoplasy\nIndex (HI)") +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Allele Frequency (af)", y = "Homoplasy Index (HI)",
       title = "Allele Frequency vs Homoplasy Index") +
  base_theme

# ── Panel B: AF distribution ─────────────────────────────────────
p_af <- HI_df |>
  ggplot(aes(x = af)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue", color = "white", alpha = 0.8) +
  labs(x = "Allele Frequency (af)", y = "Density",
       title = "Allele Frequency") +
  base_theme

# ── Panel C: HI distribution ─────────────────────────────────────
p_hi <- HI_df |>
  ggplot(aes(x = HI)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "darkorange", color = "white", alpha = 0.8) +
  labs(x = "Homoplasy Index (HI)", y = "Density",
       title = "Homoplasy Index") +
  base_theme

# ── Combine ──────────────────────────────────────────────────────
Figure1 <- p_scatter / (p_af | p_hi) +
  plot_annotation(
    subtitle = paste0("n = ", n_snps, " SNPs"),
    theme    = theme(plot.subtitle = element_text(size = 10, color = "grey40"))
  )

Figure1

# ggsave(filename = paste0(path, "Simulation Analysis/Analysis/outputs/Figure1.tiff"),plot = Figure1,width = 10, height = 10,dpi = 300)
# ggsave(filename = paste0(path, "Simulation Analysis/Analysis/outputs/Figure1.png"),plot     = Figure1, width    = 10,height   = 10, dpi      = 300)

#########################
#       Figure 2        #
#########################

# ── Parameters ───────────────────────────────────────────────────
N_EOD   <- 623
N_LOD   <- 917
N_Total <- N_EOD + N_LOD

# ── Compute expected OR ──────────────────────────────────────────
OR_df <- expand.grid(
  AF             = seq(0.05, 0.95, 0.05),
  EOD_LOD_ratio  = c(1/4, 1/3, 2/3, 3/4)
) |>
  mutate(
    N_1      = round(AF * N_Total),
    N_0      = N_Total - N_1,
    N_EOD_1  = EOD_LOD_ratio * N_1,
    N_LOD_1  = N_1  - N_EOD_1,
    N_EOD_0  = N_EOD - N_EOD_1,
    N_LOD_0  = N_LOD - N_LOD_1,
    OR       = (N_LOD_1 / N_EOD_1) / (N_LOD_0 / N_EOD_0),
    ratio_label = factor(case_when(
      EOD_LOD_ratio == 1/4 ~ "EOD:LOD = 1:3 (LOD strongly enriched)",
      EOD_LOD_ratio == 1/3 ~ "EOD:LOD = 1:2 (LOD moderately enriched)",
      EOD_LOD_ratio == 2/3 ~ "EOD:LOD = 2:1 (EOD moderately enriched)",
      EOD_LOD_ratio == 3/4 ~ "EOD:LOD = 3:1 (EOD strongly enriched)"
    ), levels = c(
      "EOD:LOD = 1:3 (LOD strongly enriched)",
      "EOD:LOD = 1:2 (LOD moderately enriched)",
      "EOD:LOD = 2:1 (EOD moderately enriched)",
      "EOD:LOD = 3:1 (EOD strongly enriched)"
    ))
  ) |>
  filter(is.finite(OR))

# ── Figure 2 ─────────────────────────────────────────────────────
Figure2 <- OR_df |>
  ggplot(aes(x = AF, y = OR, colour = ratio_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), labels = scales::percent) +
  scale_colour_manual(
    values = c(
      "EOD:LOD = 1:3 (LOD strongly enriched)"    = "#d73027",
      "EOD:LOD = 1:2 (LOD moderately enriched)"  = "#fc8d59",
      "EOD:LOD = 2:1 (EOD moderately enriched)"  = "#4575b4",
      "EOD:LOD = 3:1 (EOD strongly enriched)"    = "#313695"
    )
  ) +
  labs(
    x      = "Allele Frequency",
    y      = "Crude Odds Ratio (log scale)",
    colour = "EOD:LOD ratio among carriers",
    title  = "Expected crude OR by allele frequency and EOD:LOD ratio among carriers"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
Figure2

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure2.tiff"),plot = Figure2, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure2.png"),plot = Figure2, width = 8, height = 5, dpi = 300)


#########################
#  Set 1 - Figures 3-8  #
#########################


# ── Data loading and processing ──────────────────────────────────
sim_results_raw <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/sim_GWAS_results_set1.txt"),
  header = FALSE
)

colnames(sim_results_raw) <- c("AF_range", "MUT", "EOD_LOD_ratio", "simulation",
                               "correction", "variant", "af", "filter-pvalue",
                               "lrt-pvalue", "beta", "beta-std-err")

# ── Data loading and processing ──────────────────────────────────
sim_results_raw_highfreq <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/sim_GWAS_results_set1_highfreq.txt"),
  header = FALSE
)

colnames(sim_results_raw_highfreq) <- c("AF_range", "MUT", "EOD_LOD_ratio", "simulation",
                               "correction", "variant", "af", "filter-pvalue",
                               "lrt-pvalue", "beta", "beta-std-err")


sim_results_raw_full <-rbind(sim_results_raw, sim_results_raw_highfreq)


# ── Method and EOD:LOD factor levels ────────────────────────────
method_levels <- c("CC", "SC", "Mash distances",
                   "SNP-based kinship matrix", "Unitig-based kinship matrix",
                   "COG-based kinship matrix", "RAxML-NG", "IQ-TREE",
                   "FastTree SP", "FastTree DP FEM", "FastTree DP LMM",
                   "FastTree DP LMM (Gubbins)", "No correction")

ratio_levels <- c("1:2", "1:3", "2:1", "3:1")

sim_results <- as_tibble(sim_results_raw_full) |>
  filter(AF_range != "freq") |>
  mutate(
    EOD_ratio = as.numeric(sapply(strsplit(EOD_LOD_ratio, ":"), "[", 1)),
    LOD_ratio = as.numeric(sapply(strsplit(EOD_LOD_ratio, ":"), "[", 2))
  ) |>
  mutate(
    EOD_1 = round(af * (917 + 623) * (EOD_ratio / (EOD_ratio + LOD_ratio))),
    EOD_0 = 623 - EOD_1,
    LOD_1 = round(af * (917 + 623) * (LOD_ratio / (EOD_ratio + LOD_ratio))),
    LOD_0 = 917 - LOD_1,
    OR    = round((LOD_1 / (EOD_1 + LOD_1)) / (LOD_0 / (EOD_0 + LOD_0)), 1),
    sig   = ifelse(`lrt-pvalue` < 2.55e-6, 1, 0)
  )|>
  mutate(sig= ifelse(is.na(sig), 0, sig))|>
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
                        "FasttreeSP"         = "FastTree SP",
                        "FasttreeDP_FEM"         = "FastTree DP FEM",
                        "FasttreeDP_LMM"         = "FastTree DP LMM",
                        "FasttreeDP_LMM_Gubbins" = "FastTree DP LMM (Gubbins)",
                        "nocorrection"           = "No correction"
    ),
    correction    = factor(correction, levels = method_levels),
    EOD_LOD_ratio = factor(EOD_LOD_ratio, levels = ratio_levels)
  ) |> 
  left_join(homoplasy_results |> select(variant, HI))|>
  mutate(HI_bin = case_when(
    HI >= 0.90              ~ "[0.90, 1.00]",
    HI >= 0.80 & HI < 0.90 ~ "[0.80, 0.90)",
    HI >= 0.70 & HI < 0.80 ~ "[0.70, 0.80)",
    HI >= 0.60 & HI < 0.70 ~ "[0.60, 0.70)",
    HI >= 0.45 & HI < 0.60 ~ "[0.45, 0.60)",
    HI <  0.10              ~ "[0.00, 0.10)"
  )) |>
  mutate(HI_bin = factor(HI_bin, levels = c(
    "[0.00, 0.10)",
    "[0.60, 0.70)",
    "[0.70, 0.80)",
    "[0.80, 0.90)",
    "[0.90, 1.00]"
  )))


# ── Figure 3: overall power heatmap - EOD:LOD Ratio ─────────────────────────────
power_EOD_LOD_ratio_df <- sim_results |>
  group_by(correction, EOD_LOD_ratio) |>
  summarise(power = signif(mean(sig == 1), 2), .groups = "drop")

Figure3 <- power_EOD_LOD_ratio_df |>
  ggplot(aes(x = correction, y = EOD_LOD_ratio, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(power, 2)),
            color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Average\nproportion\nsignificant",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "EOD:LOD Ratio",
    title    = "GWAS Detection Power",
    subtitle = "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across allele frequencies 0.05–0.50 in 0.05 steps (20 simulations each)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  )
Figure3

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure3.tiff"),plot = Figure3, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure3.png"),plot = Figure3, width = 8, height = 5, dpi = 300)


# ── Figure 4: overall power heatmap  - Allele Frequency" ─────────────────────────────
power_af_df <- sim_results |>
  group_by(correction, AF_range) |>
  summarise(power = signif(mean(sig == 1), 2), .groups = "drop")

Figure4 <- power_af_df |>
  ggplot(aes(x = correction, y = AF_range, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(power, 2)),
            color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Average\nproportion\nsignificant",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Allele Frequency",
    title    = "GWAS Detection Power",
    subtitle = "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across all EOD:LOD ratios"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  )
Figure4

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure4.tiff"),plot = Figure4, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure4.png"),plot = Figure4, width = 8, height = 5, dpi = 300)


# ── Figure 5 : AF x EOD:LOD ratio power heatmap  - No correction vs Fasttree DP LMM " ─────────────────────────────
# ── shared theme ────────────────────────────────────────────────
heatmap_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 11, color = "grey40"),
    legend.position  = "right"
  )

# ── no correction (binary fill) ─────────────────────────────────
p1 <- sim_results |>
  filter(correction == "No correction") |>
  group_by(AF_range, EOD_LOD_ratio) |>
  summarise(power = mean(sig == 1), OR = first(OR), .groups = "drop") |>
  ggplot(aes(x = AF_range, y = EOD_LOD_ratio, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(100 * power), "%")),
            color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name = "Power",
                      labels = scales::percent,
                      breaks = c(0, 0.25, 0.5, 0.75, 1),
                      limits = c(0, 1)) +
  labs(
    x        = "Allele Frequency Bin",
    y        = "EOD:LOD Ratio Among Carriers",
    title    = "No population structure correction",
    subtitle = "% of simulations reaching significance"
  ) +
  heatmap_theme

# ── FastTree DP LMM ─────────────────────────────────────────────
p2 <- sim_results |>
  filter(correction == "FastTree DP LMM") |>
  group_by(AF_range, EOD_LOD_ratio) |>
  summarise(power = mean(sig == 1), OR = first(OR), .groups = "drop") |>
  ggplot(aes(x = AF_range, y = EOD_LOD_ratio, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(100 * power), "%")),
            color = "black", size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name = "Power",
                      labels = scales::percent,
                      breaks = c(0, 0.25, 0.5, 0.75, 1),
                      limits = c(0, 1)) +
  labs(
    x        = "Allele Frequency Bin",
    y        = "EOD:LOD Ratio Among Carriers",
    title    = "FastTree DP LMM correction",
    subtitle = "% of simulations reaching significance"
  ) +
  heatmap_theme

# ── combine ─────────────────────────────────────────────────────
library(patchwork)
Figure5 <- p1 + p2 + plot_layout(guides = "collect") +
  plot_annotation(
    title   = "GWAS detection power by allele frequency and effect size",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )
Figure5

ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure5.tiff"),plot = Figure5, width = 20, height = 7, dpi = 300)
ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure5.png"),plot = Figure5, width = 20, height = 7, dpi = 300)


# ── Figure 6 : overall power heatmap  - Homoplasy Index ─────────────────────────────
# ── n SNPs per HI bin for y-axis labels ─────────────────────────
n_per_hibin <- sim_results |>
  group_by(HI_bin) |>
  summarise(n_snps = n_distinct(MUT), .groups = "drop")

hi_labels <- setNames(
  paste0(levels(sim_results$HI_bin), "\n(n=", n_per_hibin$n_snps, " SNPs)"),
  levels(sim_results$HI_bin)
)

# ── Summarise power ──────────────────────────────────────────────
power_HI_df <- sim_results |>
  filter(!is.na(HI_bin)) |>
  group_by(correction, HI_bin) |>
  summarise(power = round(mean(sig == 1), 2), .groups = "drop")

# ── Figure 6 ─────────────────────────────────────────────────────
Figure6 <- power_HI_df |>
  ggplot(aes(x = correction, y = HI_bin, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", power)),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Average\nproportion\nsignificant",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = hi_labels) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Homoplasy Index (HI) Bin",
    title    = "GWAS Detection Power by Homoplasy Level",
    subtitle = "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across all allele frequencies and EOD:LOD ratios | 6 quantile-based HI bins"
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
Figure6

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure6.tiff"),plot = Figure6, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure6.png"),plot = Figure6, width = 8, height = 5, dpi = 300)


# ── Figure 7 ─────────────────────────────────────────────────────
Figure7 <- sim_results |> 
  select(MUT, af, HI_bin)|>
  group_by(MUT, HI_bin)|>
  summarise(af_mean= mean(af, na.rm=T))|>
  ggplot(aes(x = HI_bin, y = af_mean)) +
  geom_boxplot(width = 0.15, outlier.shape = NA,
               alpha = 0.7, color = "grey30", fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5, color = "grey20") +
  stat_summary(fun = median, geom = "text",
               aes(label = round(after_stat(y), 2)),
               vjust = -0.8, color = "red", fontface = "bold", size = 3.5) +
  # n labels at top
  geom_text(
    data = n_per_hibin,
    aes(x = HI_bin, y = 0.75,
        label = paste0("n=", n_snps)),
    inherit.aes = FALSE, size = 3.5, color = "grey30"
  ) +
  labs(
    x        = "Homoplasy Index Bin",
    y        = "Allele Frequency",
    title    = "Allele frequency distribution by homoplasy index bin",
    subtitle = "Median allele frequency shown in red; n = number of SNPs per bin"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey40", size = 11)
  )
Figure7

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure7.tiff"),plot = Figure7, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure7.png"),plot = Figure7, width = 8, height = 5, dpi = 300)


# ── Figure 8 ─────────────────────────────────────────────────────
# ── Shared elements ──────────────────────────────────────────────
heatmap_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 11, color = "grey40"),
    strip.text.y     = element_text(size = 9, face = "bold"),
    legend.position  = "right"
  )

shared_fill <- scale_fill_gradient(
  low    = "white", high = "darkred",
  name   = "Power",
  labels = scales::percent,
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1)
)

make_hi_af_heatmap <- function(correction_name, title) {
  
  # n SNPs per HI bin
  n_per_bin <- sim_results |>
    filter(!is.na(HI_bin)) |>
    group_by(HI_bin) |>
    summarise(n_snps = n_distinct(MUT), .groups = "drop") |>
    mutate(HI_bin_label = paste0(HI_bin, "\n(n=", n_snps, " SNPs)"))
  
  # join labels and reorder factor (reversed)
  sim_results |>
    filter(correction == correction_name, !is.na(HI_bin)) |>
    left_join(n_per_bin, by = "HI_bin") |>
    mutate(HI_bin_label = factor(HI_bin_label,
                                 levels = rev(n_per_bin$HI_bin_label))) |>
    group_by(HI_bin_label, AF_range, EOD_LOD_ratio) |>
    summarise(power = round(mean(sig == 1), 2), .groups = "drop") |>
    mutate(label_color = ifelse(power > 0.6, "white", "black")) |>
    ggplot(aes(x = AF_range, y = EOD_LOD_ratio, fill = power)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(100 * power), "%"),
                  color = label_color),
              size = 3.5, fontface = "bold") +
    scale_color_identity() +
    shared_fill +
    facet_grid(HI_bin_label ~ ., switch = "y") +
    labs(
      x        = "Allele Frequency Bin",
      y        = "EOD:LOD Ratio Among Carriers",
      title    = title,
      subtitle = "% of simulations reaching significance"
    ) +
    heatmap_theme +
    theme(
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.placement   = "outside"
    )
}

# ── Build panels ─────────────────────────────────────────────────
p1 <- make_hi_af_heatmap("No correction",   "No population structure correction")
p2 <- make_hi_af_heatmap("FastTree DP LMM", "FastTree DP LMM correction")

# ── Combine ──────────────────────────────────────────────────────
Figure8 <- p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "GWAS detection power by allele frequency, effect size, and homoplasy index",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
Figure8

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure8.tiff"),plot = Figure8, width = 20, height = 10, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure8.png"),plot = Figure8, width = 20, height = 10, dpi = 300)

#########################
#  Set 2 - Figures 9-11 #
#########################


# ── Data loading and processing ──────────────────────────────────
sim_results_set2_raw <- read.delim(
  paste0(path, "Simulation Analysis/Analysis/inputs/sim_GWAS_results_set2.txt"),
  header = FALSE) [-1,]

colnames(sim_results_set2_raw) <- c("AF_range","HI_range",  "MUT", "EOD_LOD_ratio", "simulation",
                               "correction", "variant", "af", "filter-pvalue",
                               "lrt-pvalue", "beta", "beta-std-err")

sim_results_set2_df <- sim_results_set2_raw |> 
  left_join(homoplasy_results |> select(variant, HI)) |>
  mutate(EOD_LOD_ratio = factor(EOD_LOD_ratio, levels = c("1:2", "1:3", "2:1", "3:1")))|>  
  mutate(sig   = ifelse(`lrt-pvalue` < 2.55e-6, 1, 0)
  )|>
  mutate(sig= ifelse(is.na(sig), 0, sig))|>
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
                        "FasttreeSP"         = "FastTree SP",
                        "FasttreeDP_FEM"         = "FastTree DP FEM",
                        "FasttreeDP_LMM"         = "FastTree DP LMM",
                        "FasttreeDP_LMM_Gubbins" = "FastTree DP LMM (Gubbins)",
                        "nocorrection"           = "No correction"
    ),
    correction    = factor(correction, levels = method_levels),
    EOD_LOD_ratio = factor(EOD_LOD_ratio, levels = ratio_levels)
  ) |> 
  mutate(HI_bin = case_when(
    HI_range == "90_100" ~ "[0.90, 1.00]",
    HI_range == "80_90" ~ "[0.80, 0.90)",
    HI_range == "70_80" ~ "[0.70, 0.80)",
    HI_range == "60_70" ~ "[0.60, 0.70)",
    HI_range == "45_55" ~ "[0.45, 0.60)",
    HI_range == "00_10" ~ "[0.00, 0.10)"
  ))|>
  mutate(HI_bin = factor(HI_bin, levels = c(
    "[0.00, 0.10)",
    "[0.45, 0.60)",
    "[0.60, 0.70)",
    "[0.70, 0.80)",
    "[0.80, 0.90)",
    "[0.90, 1.00]"
  )))


# ── Figure 9 : overall power heatmap  - Homoplasy Index ─────────────────────────────
# ── n SNPs per HI bin for y-axis labels ─────────────────────────
n_per_hibin_set2 <- sim_results_set2_df |>
  group_by(HI_bin) |>
  summarise(n_snps = n_distinct(MUT), .groups = "drop")

hi_labels_set2 <- setNames(
  paste0(levels(sim_results_set2_df$HI_bin), "\n(n=", n_per_hibin_set2$n_snps, " SNPs)"),
  levels(sim_results_set2_df$HI_bin)
)

# ── Summarise power ──────────────────────────────────────────────
power_HI_set2_df <- sim_results_set2_df |>
  filter(!is.na(HI_bin)) |>
  group_by(correction, HI_bin) |>
  summarise(power = round(mean(sig == 1), 2), .groups = "drop")

# ── Figure 9 ─────────────────────────────────────────────────────
Figure9 <- power_HI_set2_df |>
  ggplot(aes(x = correction, y = HI_bin, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", power)),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Average\nproportion\nsignificant",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = hi_labels_set2) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Homoplasy Index (HI) Bin",
    title    = "GWAS Detection Power by Homoplasy Level",
    subtitle = "Colour = mean proportion of simulations with lrt-pvalue < 2.55×10⁻⁶\naveraged across all allele frequencies and EOD:LOD ratios | 6 quantile-based HI bins"
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
Figure9

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure9.tiff"),plot = Figure9, width = 8, height = 5, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure9.png"),plot = Figure9, width = 8, height = 5, dpi = 300)


# ── Table 1: Low vs High HI power difference per method ─────────
table1 <- power_HI_set2_df |>
  filter(HI_bin %in% c("[0.00, 0.10)", "[0.90, 1.00]")) |>
  mutate(HI_label = ifelse(HI_bin == "[0.00, 0.10)",
                           "Low HI (0.00–0.10)",
                           "High HI (0.90–1.00)")) |>
  select(correction, HI_label, power) |>
  tidyr::pivot_wider(names_from = HI_label, values_from = power) |>
  mutate(
    Diff_num   = `High HI (0.90–1.00)` - `Low HI (0.00–0.10)`,
    `Low HI (0.00–0.10)`  = sprintf("%.2f", `Low HI (0.00–0.10)`),
    `High HI (0.90–1.00)` = sprintf("%.2f", `High HI (0.90–1.00)`)
  ) |>
  arrange(Diff_num) |>
  mutate(Difference = ifelse(Diff_num >= 0,
                             paste0("+", sprintf("%.2f", Diff_num)),
                             sprintf("%.2f", Diff_num))) |>
  select(-Diff_num) |>
  rename(Method = correction)

table1




# ── Summarise power by method, HI bin, EOD:LOD ratio ────────────
power_HI_AF_EOD_LOD_ratio_df <- sim_results_set2_df |>
  group_by(correction, HI_bin, AF_range , EOD_LOD_ratio) |>
  summarise(power = mean(sig == 1), .groups = "drop") |>
  mutate(
    EOD_LOD_label = recode(EOD_LOD_ratio,
                           "1:2" = "1:2 (LOD moderately enriched)",
                           "1:3" = "1:3 (LOD strongly enriched)",
                           "2:1" = "2:1 (EOD moderately enriched)",
                           "3:1" = "3:1 (EOD strongly enriched)"
    ),
    EOD_LOD_label = factor(EOD_LOD_label, levels = c(
      "1:2 (LOD moderately enriched)",
      "1:3 (LOD strongly enriched)",
      "2:1 (EOD moderately enriched)",
      "3:1 (EOD strongly enriched)"
    ))
  )|>
  mutate(AF_label = recode(AF_range,
                           "0.05" = "AF 0.05",
                           "0.10" = "AF 0.10",
                           "0.15" = "AF 0.15",
                           "0.20" = "AF 0.20",
                           "0.25" = "AF 0.25",
                           "0.40" = "AF 0.40"),
         AF_label = factor(AF_label, levels = rev(c("AF 0.05","AF 0.10","AF 0.15","AF 0.20","AF 0.25","AF 0.40")))
  )

# ── Colour palette per method ────────────────────────────────────
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

# ── n SNPs per HI bin for facet labels ──────────────────────────
n_per_hibin_set2 <- sim_results_set2_df |>
  filter(!is.na(HI_bin)) |>
  group_by(HI_bin) |>
  summarise(n_snps = n_distinct(MUT), .groups = "drop")

hi_facet_labels <- setNames(
  paste0(levels(sim_results_set2_df$HI_bin), "\n(n=", n_per_hibin_set2$n_snps, " SNPs)"),
  levels(sim_results_set2_df$HI_bin)
)

# ── Figure 11 ────────────────────────────────────────────────────
Figure10 <- power_HI_AF_EOD_LOD_ratio_df |>
  ggplot(aes(x = HI_bin, y = power, color = correction, group = correction)) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "red",
             linewidth = 0.8) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(AF_label ~ EOD_LOD_label) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x        = "Homoplasy Index (HI) Bin",
    y        = "Average Proportion Significant",
    title    = "GWAS Detection Power Across AF x EOD:LOD",
    subtitle = "Each panel = one EOD:LOD ratio  |  averaged across all allele frequencies  |  x-axis: low → high homoplasy\nRed dotted line = 80% power threshold"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position  = "right",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold")
  )

Figure10

# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure10.tiff"),plot = Figure10, width = 18, height = 10, dpi = 300)
# ggsave(paste0(path, "Simulation Analysis/Analysis/outputs/Figure10.png"),plot = Figure10, width = 18, height = 10, dpi = 300)







power_HI_AF_EOD_LOD_ratio_df |>
  filter(correction =="FastTree DP LMM")|>
  ggplot(aes(x = HI_bin, y = power, color = correction, group = correction)) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "red",
             linewidth = 0.8) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(AF_label ~ EOD_LOD_label) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x        = "Homoplasy Index (HI) Bin",
    y        = "Average Proportion Significant",
    title    = "GWAS Detection Power Across AF x EOD:LOD",
    subtitle = "Each panel = one EOD:LOD ratio  |  averaged across all allele frequencies  |  x-axis: low → high homoplasy\nRed dotted line = 80% power threshold"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position  = "right",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold")
  )
















################################
# Permutation Analysis
# False positives 
################################



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
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(fpr, 2)),
            color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Empirical\nfalse positive\nrate",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Allele Frequency",
    title    = "GWAS False Positive Rate\by Allele Frequency",
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
  geom_text(aes(label = sprintf("%.2f", fpr)),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "darkred",
                      name   = "Empirical\nfalse positive\nrate",
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = hi_labels_null) +
  labs(
    x        = "Population Structure Correction Method",
    y        = "Homoplasy Index (HI) Bin",
    title    = "GWAS False Positive Rate\by Homoplasy Level",
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
