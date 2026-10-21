################################################################################
# Empirical analysis: sensitivity of GBS GWAS results to population structure
# correction methods (pyseer outputs from 1,540 U.S. neonatal invasive isolates)
#
# Generates (in empirical_analysis/outputs/):
#   Figure_2_QQplots_population_structure_methods.tiff
#   Table_1_N_significant_hits_per_method.csv
#   Figure_S3_Phylogenies.tiff
#   Table_S1_Recombination_rates_per_SC.csv
#   Figure_S4_Recombination_phylogenies.tiff
#   Table_S2_GWAS_hits_with_vs_without_recombination.csv
#   Supplementary_Data_pyseer_significant_hits.xlsx
#
# Run from RStudio (Source) or from the command line:
#   Rscript empirical_analysis/analysis_code/empirical_analysis.R
################################################################################

rm(list = ls())

# 0. SET UP ---------------------------------------------------------------------

## Paths (relative to this script: <repo>/empirical_analysis/analysis_code/)
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
}
dir_path <- dirname(get_script_dir())
path     <- file.path(dir_path, "inputs")
outpath  <- file.path(dir_path, "outputs")
dir.create(outpath, showWarnings = FALSE)

## Packages (CRAN + Bioconductor)
use_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in c("tidyverse", "patchwork", "openxlsx", "ape", "phytools",
              "treeio", "ggtree", "ggpubr")) use_pkg(pkg)

## Helper: save TIFF with LZW compression (lossless, much smaller files)
save_tiff <- function(plot, filename, width, height) {
  ggsave(file.path(outpath, filename), plot, width = width, height = height,
         dpi = 300, device = "tiff", compression = "lzw")
}

## Population structure correction methods (folder names in inputs/pyseer/)
methods <- tribble(
  ~method,                            ~model, ~table_label,                                               ~fig_label,
  "CC",                               "FEM",  "CC",                                                       "CC",
  "SC",                               "FEM",  "SC",                                                       "SC",
  "Mash_distances",                   "FEM",  "Mash distances",                                           "Mash\nDistances",
  "SNP_kinship",                      "LMM",  "SNP-based kinship matrix",                                 "SNP\nKinship",
  "Unitig_kinship",                   "LMM",  "Unitig-based kinship matrix",                              "Unitig\nKinship",
  "COG_kinship",                      "LMM",  "COG-based kinship matrix",                                 "COG\nKinship",
  "RAxML",                            "LMM",  "SNP-based Phylogeny (RAxML-NG)",                           "RAxML",
  "IQTREE",                           "LMM",  "SNP-based Phylogeny (IQ-TREE)",                            "IQtree",
  "FastTree_SP",                      "LMM",  "Whole Genome Phylogeny (FastTree SP)",                     "Fasttree SP",
  "FastTree_DP_FEM",                  "FEM",  "Whole Genome Phylogeny (FastTree DP)",                     "Fasttree DP\n(FEM)",
  "FastTree_DP_LMM",                  "LMM",  "Whole Genome Phylogeny (FastTree DP)",                     "Fasttree DP\n(LMM)",
  "FastTree_DP_LMM_noRecombination",  "LMM",  "Recombination-free Whole Genome Phylogeny (FastTree DP)",  "Fasttree DP\n(Recombination-free)"
)
gwas_types <- c("SNPGWAS", "DBGWAS", "PanGWAS")

## Bonferroni thresholds: 0.05 / number of tested variants (af 0.05-0.95);
## for DBGWAS, the number of unique unitig patterns (87,899)
N_DBGWAS_PATTERNS <- 87899
sig_threshold <- function(data, gwas) {
  if (gwas == "DBGWAS") 0.05 / N_DBGWAS_PATTERNS else 0.05 / nrow(data)
}

## Helpers
read_pyseer <- function(file) {
  read.table(file, sep = "\t", header = TRUE) |>
    filter(af >= 0.05, af < 0.95)
}

add_beta_ci <- function(df) {
  df |>
    mutate(beta_lower = beta - 1.96 * beta.std.err,
           beta_upper = beta + 1.96 * beta.std.err,
           beta_CI    = sprintf("%.2f (%.2f; %.2f)", beta, beta_lower, beta_upper),
           P_value    = formatC(lrt.pvalue, format = "e", digits = 2))
}

lambda_gc <- function(df) {
  p <- df$lrt.pvalue
  p <- p[!is.na(p) & p > 0]
  round(median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1), 2)
}

QQ_plot <- function(data) {
  df <- data.frame(expected = -log10(ppoints(nrow(data))),
                   observed = -log10(sort(data$lrt.pvalue)))
  lim <- max(c(df$expected, df$observed))
  ggplot(df, aes(x = expected, y = observed)) +
    geom_point() +
    geom_abline(slope = 1, color = "red", lty = "dashed") +
    labs(x = "-log10(Expected P)", y = "-log10(Observed P)") +
    theme_minimal() +
    xlim(0, lim) + ylim(0, lim)
}


# 1. LOAD PYSEER RESULTS --------------------------------------------------------

gwas_results <- list()
for (m in methods$method) {
  for (g in gwas_types) {
    data <- read_pyseer(file.path(path, "pyseer/population_structure_methods", m, paste0(g, ".txt")))
    gwas_results[[m]][[g]] <- list(
      data    = data,
      sighits = data |> filter(lrt.pvalue < sig_threshold(data, g)),
      lambda  = lambda_gc(data)
    )
  }
}

## Genomic inflation factors (shown on Figure 2)
lambdaGC_df <- map_dfr(methods$method, function(m) {
  tibble(Method = m,
         SNPGWAS = gwas_results[[m]]$SNPGWAS$lambda,
         PanGWAS = gwas_results[[m]]$PanGWAS$lambda,
         DBGWAS  = gwas_results[[m]]$DBGWAS$lambda)
})
print(lambdaGC_df)


# 2. FIGURE 2: QQ PLOTS PER POPULATION STRUCTURE METHOD -------------------------

add_lambda <- function(p, value) {
  p + annotate("text", x = Inf, y = -Inf,
               label = sprintf('lambda[GC] == "%.2f"', value), parse = TRUE,
               hjust = 1.1, vjust = -0.5, size = 5)
}

qq_plots <- list()
for (g in gwas_types) {
  for (m in methods$method) {
    res <- gwas_results[[m]][[g]]
    qq_plots[[length(qq_plots) + 1]] <- add_lambda(QQ_plot(res$data), res$lambda)
  }
}

text_plot <- function(label, angle = 0) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 5,
             fontface = "bold", angle = angle) +
    theme_void()
}
title_row   <- wrap_plots(lapply(methods$fig_label, text_plot), ncol = nrow(methods))
left_column <- plot_spacer() /
  wrap_plots(lapply(c("SNP-GWAS", "DBGWAS", "PanGWAS"), text_plot, angle = 90), ncol = 1) +
  plot_layout(heights = c(1, 3))
main_column <- title_row / wrap_plots(qq_plots, ncol = nrow(methods), nrow = 3) +
  plot_layout(heights = c(0.25, 3))

Figure2 <- wrap_plots(list(left_column, main_column), ncol = 2, widths = c(0.25, 12))
save_tiff(Figure2, "Figure_2_QQplots_population_structure_methods.tiff", width = 30, height = 10)


# 3. TABLE 1: NUMBER OF SIGNIFICANT HITS PER METHOD -----------------------------

Table1 <- methods |>
  mutate(`SNP-GWAS` = map_int(method, ~ nrow(gwas_results[[.x]]$SNPGWAS$sighits)),
         PanGWAS    = map_int(method, ~ nrow(gwas_results[[.x]]$PanGWAS$sighits)),
         DBGWAS     = map_int(method, ~ nrow(gwas_results[[.x]]$DBGWAS$sighits))) |>
  select(Method = table_label, Model = model, `SNP-GWAS`, PanGWAS, DBGWAS)
print(Table1)
write.csv(Table1, file.path(outpath, "Table_1_N_significant_hits_per_method.csv"), row.names = FALSE)


# 4. FIGURE S3: PHYLOGENIES BUILT WITH FOUR METHODS -----------------------------

## Metadata
serotypes   <- read.table(file.path(path, "typing/Serotypes.txt"), header = TRUE)
CCs         <- read.csv(file.path(path, "typing/CC.txt"))
SCs         <- read.csv(file.path(path, "typing/SC.csv")) |> rename(Sample = Taxon, SC = Cluster)
typing_data <- serotypes |> left_join(CCs) |> left_join(SCs)

SC_colors <- c("#dbdb8d", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
               "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#f7b6d2", "#c5b0d5", "#c49c94", "black")

## Trees (midpoint rooted)
trees <- list(
  "RaxML\nSNP-Phylogeny"                 = "phylogenies/RAxML.tre",
  "IQ-tree\nSNP-Phylogeny"               = "phylogenies/IQTREE.tre",
  "FastTree SP\nWhole-genome Phylogeny"  = "phylogenies/FastTree_SP.tre",
  "FastTree DP\nWhole-genome Phylogeny"  = "phylogenies/FastTree_DP.tre"
) |>
  purrr::map(~ midpoint.root(read.tree(file.path(path, .x))))

## SCs with <=5 isolates grouped as "other"
fewer5_SC <- which(table(typing_data$SC) <= 5)
metadata_tree_S3 <- typing_data |>
  mutate(SC2 = ifelse(SC %in% fewer5_SC, "other", SC),
         SC2 = factor(SC2, levels = c(1:12, "other"))) |>
  rename(label = Sample) |>
  filter(label %in% trees[[1]]$tip.label)

plot_tree <- function(tree, metadata, colors) {
  tree_metadata <- full_join(keep.tip(tree, metadata$label), metadata, by = "label")
  ggtree(tree_metadata, ladderize = TRUE, lwd = 0.1) +
    geom_tippoint(aes(color = as.factor(SC2)), size = 1.5) +
    scale_colour_manual(name = "SC", values = colors)
}
tree_theme <- theme(legend.position = "none",
                    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
legend_theme <- theme(legend.position = "bottom", legend.direction = "horizontal",
                      legend.box = "horizontal",
                      legend.text = element_text(size = 15), legend.title = element_text(size = 15))

S3_plots <- imap(trees, ~ plot_tree(.x, metadata_tree_S3, SC_colors))
S3_legend <- as_ggplot(get_legend(
  S3_plots[[1]] + legend_theme + guides(color = guide_legend(nrow = 2, override.aes = list(size = 5)))
))
FigureS3 <- wrap_plots(imap(S3_plots, ~ .x + tree_theme + ggtitle(.y)), ncol = 2, nrow = 2) /
  S3_legend + plot_layout(heights = c(4, 1))
save_tiff(FigureS3, "Figure_S3_Phylogenies.tiff", width = 10, height = 15)


# 5. TABLE S1: RECOMBINATION RATES PER SEQUENCE CLUSTER -------------------------

## r/m       = SNPs inside recombinations / SNPs outside recombinations
## rho/theta = recombination blocks / SNPs outside recombinations
recomb_rates <- function(file) {
  d <- read.delim(file.path(path, "gubbins", file))
  c(r.m       = round(sum(d$Number.of.SNPs.Inside.Recombinations) /
                        sum(d$Number.of.SNPs.Outside.Recombination), 2),
    rho.theta = round(sum(d$Number.of.Recombination.Blocks) /
                        sum(d$Number.of.SNPs.Outside.Recombination), 2))
}

## SC-specific reference genomes used for the per-SC alignments
SC_references <- c(SC1 = "B105 (CP021773.1)",      SC2 = "MA07 (CP134688.1)",
                   SC3 = "HU-GS5823 (AP018935.1)", SC4 = "Sag153 (CP036376.1)",
                   SC5 = "NGBS357 (CP012503.1)",   SC6 = "BCJB1835 (CP101996.1)",
                   SC7 = "BJ01 (CP059383.1)",      SC8 = "GBS1-NY (CP007570.1)",
                   SC9 = "GBSIR0001 (CP120953.1)", SC10 = "FDAARGOS_512 (CP033822.1)")

TableS1 <- map_dfr(c("Fulltree", paste0("SC", 1:10)), function(sc) {
  ref_AP <- if (sc == "Fulltree") recomb_rates("gubbins_fulltree.per_branch_statistics.csv")
            else recomb_rates(paste0("gubbins_", sc, "_out.per_branch_statistics.csv"))
  ref_SC <- if (sc == "Fulltree") c(r.m = NA, rho.theta = NA)
            else recomb_rates(paste0("gubbins_REFd_", sc, "_out.per_branch_statistics.csv"))
  tibble(SC = sc,
         `r/m (AP018935.1)`        = ref_AP[["r.m"]],
         `rho/theta (AP018935.1)`  = ref_AP[["rho.theta"]],
         `r/m (SC reference)`      = ref_SC[["r.m"]],
         `rho/theta (SC reference)` = ref_SC[["rho.theta"]],
         Reference = if (sc == "Fulltree") NA else SC_references[[sc]])
})
print(TableS1)
write.csv(TableS1, file.path(outpath, "Table_S1_Recombination_rates_per_SC.csv"), row.names = FALSE)


# 6. FIGURE S4: FASTTREE DP PHYLOGENY WITH VS WITHOUT RECOMBINATION -------------

## SCs with <=7 isolates grouped as "other"
fewer7_SC <- which(table(SCs$SC) <= 7)
metadata_tree_S4 <- SCs |>
  mutate(SC2 = ifelse(SC %in% fewer7_SC, "other", SC),
         SC2 = factor(SC2, levels = c(1:12, "other"))) |>
  rename(label = Sample)
SC_colors_S4 <- setNames(c(SC_colors[1:12], "black"), c(1:12, "other"))

## Branch lengths of the FastTree tree rescaled to number of substitutions
## (alignment length: 2,231,314 bp); Gubbins tree is already in substitutions
N_SITES <- 2231314
Fasttree_DP <- midpoint.root(read.tree(file.path(path, "phylogenies/FastTree_DP.tre")))
Fasttree_DP$edge.length <- Fasttree_DP$edge.length * N_SITES
Gubbins_tree <- midpoint.root(read.tree(file.path(path, "gubbins/recombination_free.tre")))

S4_plots <- list(
  "FastTree DP\nWhole-genome Phylogeny\n(With Recombination)" = Fasttree_DP,
  "FastTree DP\nWhole-genome Phylogeny\nRecombination-free"   = Gubbins_tree
) |>
  purrr::map(~ plot_tree(.x, metadata_tree_S4, SC_colors_S4) +
        geom_treescale(x = 0, y = 1400, width = 2000, color = "red",
                       fontsize = 6, linesize = 2, offset = 20))

S4_legend <- as_ggplot(get_legend(
  S4_plots[[1]] + legend_theme + guides(color = guide_legend(nrow = 2, override.aes = list(size = 5)))
))
FigureS4 <- wrap_plots(imap(S4_plots, ~ .x + tree_theme + ggtitle(.y)), ncol = 2, nrow = 1) /
  S4_legend + plot_layout(heights = c(4, 1))
save_tiff(FigureS4, "Figure_S4_Recombination_phylogenies.tiff", width = 15, height = 10)


# 7. TABLE S2: SIGNIFICANT HITS WITH VS WITHOUT RECOMBINATION -------------------

## SNP functional annotations (SnpEff)
SNP_ann <- read.csv(file.path(path, "pyseer/population_structure_methods/FastTree_DP_LMM/SNP_annotations.csv"),
                    header = FALSE,
                    col.names = c("POS", "REF", "ALT", "FUNCLASS", "IMPACT", "gene", "nt_change", "aa_change")) |>
  mutate(Variant  = paste0(REF, POS, ALT),
         aa_change = str_remove(aa_change, "p."),
         type = case_when(IMPACT == "HIGH"     ~ "Nonsense",
                          IMPACT == "MODERATE" ~ "Missense",
                          IMPACT == "LOW"      ~ "Synonymous",
                          TRUE                 ~ "non-coding"),
         Annotation = ifelse(type == "non-coding",
                             paste0("Mutation: non-coding; gene: ", gene),
                             paste0("Mutation: ", type, " (", aa_change, "); gene: ", gene))) |>
  distinct(Variant, .keep_all = TRUE) |>
  select(Variant, Annotation)

## Unitig annotations (pyseer annotate_hits output)
read_unitig_ann <- function(method) {
  read.delim(file.path(path, "pyseer/population_structure_methods", method, "DBGWAS_annotation.txt"),
             sep = "\t", header = FALSE,
             col.names = colnames(gwas_results[[method]]$DBGWAS$data)) |>
    mutate(ref_or_contig = ifelse(str_extract(notes, "^[^:]+") == "AP018935.1", "ref", "unaligned contig"),
           coord_range   = str_extract(notes, "(?<=:)[0-9]+-[0-9]+")) |>
    separate(notes, into = c("contig_position", "upstream", "gene", "downstream"),
             sep = ";", remove = FALSE, fill = "right", extra = "drop") |>
    mutate(Annotation = paste0("Coord:", coord_range, " (", ref_or_contig, "), gene:", gene)) |>
    distinct(variant, .keep_all = TRUE) |>
    select(variant, Annotation)
}

## COG annotations (manually curated)
COG_ann <- c(arnC       = "Undecaprenyl-phosphate 4-deoxy-4-formamido-L-arabinose transferase, glycosyltransferase family 2",
             group_1981 = "Glycosyltransferase (glyG)")

## Unitigs are numbered in the order they appear in the recombination-containing
## analysis; unitigs only found in the recombination-free analysis are numbered after
DB_with    <- gwas_results$FastTree_DP_LMM$DBGWAS$sighits
DB_without <- gwas_results$FastTree_DP_LMM_noRecombination$DBGWAS$sighits
unitig_ids <- c(DB_with$variant, setdiff(DB_without$variant, DB_with$variant))
unitig_ids <- setNames(paste("Unitig", seq_along(unitig_ids)), unitig_ids)

hits_table <- function(method, label) {
  res <- gwas_results[[method]]
  snp <- res$SNPGWAS$sighits |> add_beta_ci() |>
    mutate(Variant = str_split_fixed(variant, "_", 4) |> (\(x) paste0(x[, 3], x[, 2], x[, 4]))()) |>
    left_join(SNP_ann, by = "Variant") |>
    mutate(GWAS = "SNP-GWAS")
  pan <- res$PanGWAS$sighits |> add_beta_ci() |>
    mutate(Variant = variant, Annotation = unname(COG_ann[variant]), GWAS = "PanGWAS")
  db <- res$DBGWAS$sighits |> add_beta_ci() |>
    left_join(read_unitig_ann(method), by = "variant") |>
    mutate(Variant = unname(unitig_ids[variant]), GWAS = "DBGWAS") |>
    arrange(parse_number(Variant))
  bind_rows(snp, pan, db) |>
    mutate(Phylogeny = label, AF = sprintf("%.3f", af)) |>
    select(Phylogeny, GWAS, Variant, AF, `β (95% CI)` = beta_CI, `p-value` = P_value,
           `Functional Annotation` = Annotation)
}

TableS2 <- bind_rows(
  hits_table("FastTree_DP_LMM",                 "Original (recombination-containing) phylogeny"),
  hits_table("FastTree_DP_LMM_noRecombination", "Recombination-free phylogeny")
) |>
  ## * = variant significant with only one of the two phylogenies
  add_count(GWAS, Variant, name = "n_analyses") |>
  mutate(Variant = ifelse(n_analyses == 1, paste0(Variant, "*"), Variant)) |>
  select(-n_analyses)
print(as_tibble(TableS2), n = Inf, width = Inf)
write.csv(TableS2, file.path(outpath, "Table_S2_GWAS_hits_with_vs_without_recombination.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# 8. SUPPLEMENTARY DATA: ALL SIGNIFICANT PYSEER HITS PER METHOD ------------------

sheet_names <- c(CC = "CCs (FEM)", SC = "SCs (FEM)", Mash_distances = "Mash dist (FEM)",
                 RAxML = "RaxML (LMM)", IQTREE = "IQtree (LMM)", FastTree_SP = "FastTree SP (LMM)",
                 FastTree_DP_LMM = "FastTree DP (LMM)", FastTree_DP_FEM = "FastTree DP (FEM)",
                 SNP_kinship = "SNP kinship (LMM)", Unitig_kinship = "Unitig kinship (LMM)",
                 COG_kinship = "COG kinship (LMM)", FastTree_DP_LMM_noRecombination = "No recomb (LMM)")
keep_cols <- c("variant", "af", "filter.pvalue", "lrt.pvalue", "beta", "beta.std.err", "variant_h2", "notes")

wb <- createWorkbook()
for (m in names(sheet_names)) {
  for (g in gwas_types) {
    sheet <- paste(g, sheet_names[[m]])
    addWorksheet(wb, sheet)
    writeData(wb, sheet, gwas_results[[m]][[g]]$sighits |> select(any_of(keep_cols)))
  }
}

## Alternative pangenome tools (Panaroo is the default used above)
pangenome_tools <- c("PanGWAS Roary"           = "PanGWAS_ROARY.txt",
                     "PanGWAS Panaroo + CLARC" = "PanGWAS_PANAROO_CLARC.txt",
                     "PanGWAS Roary + CLARC"   = "PanGWAS_ROARY_CLARC.txt")
for (sheet in names(pangenome_tools)) {
  data <- read_pyseer(file.path(path, "pyseer/pangenome_methods", pangenome_tools[[sheet]]))
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data |> filter(lrt.pvalue < 0.05 / nrow(data)) |> select(any_of(keep_cols)))
}
saveWorkbook(wb, file.path(outpath, "Supplementary_Data_pyseer_significant_hits.xlsx"), overwrite = TRUE)

sessionInfo()
