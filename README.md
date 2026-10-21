# GBS-GWAS-Methods-Sensitivity-Analyses

Code and data to reproduce the analyses in *"Benchmarking Population Structure Correction Methods in Bacterial Genome Wide Association Study: A Simulation-Based Study Using an Adaptable, Open-source Pipeline"* (Balasubramanian, Cavalli et al.). The paper assesses how sensitive group B *Streptococcus* (GBS) GWAS results are to population structure correction methods and pangenome tools, using:

- **Empirical analyses:** pyseer SNP-, Pan- and DBGWAS on 1,540 U.S. neonatal invasive GBS isolates (EOD vs LOD/VLOD), run under 12 population structure correction strategies.
- **Simulation analyses:** power and false-positive rates of the same correction methods, estimated on simulated phenotypes across allele frequencies, effect sizes and homoplasy indices.

The upstream GWAS pipeline (read QC, assembly, phylogenies, pyseer) is available at <https://github.com/Leacavalli/GBS-GWAS-pipeline>.

## Repository structure

```
empirical_analysis/
  analysis_code/empirical_analysis.R   # Figure 2, Table 1, Figures S3–S4, Tables S1–S2, Supplementary Data
  inputs/
    phylogenies/                       # RAxML-NG, IQ-TREE, FastTree SP/DP trees
    gubbins/                           # Gubbins recombination-free tree + per-branch statistics (full tree, per SC)
    typing/                            # serotypes, clonal complexes, PopPUNK sequence clusters
    pyseer/
      population_structure_methods/<method>/{SNPGWAS,PanGWAS,DBGWAS}.txt
      pangenome_methods/               # PanGWAS with Roary, Roary+CLARC, Panaroo+CLARC
  outputs/
simulation_analysis/
  simulation_code/                     # SLURM scripts: homoplasy index + simulated phenotypes (see README there)
  analysis_code/simulation_analysis.R  # Figures 3–5, Table 2, Figures S1–S2, S5–S10
  inputs/                              # compiled pyseer results on simulated / permuted phenotypes, homoplasy index
  outputs/
```

## Reproducing the figures and tables

Run each script from RStudio (open it, then Source), or from the repository root:

```bash
Rscript empirical_analysis/analysis_code/empirical_analysis.R
Rscript simulation_analysis/analysis_code/simulation_analysis.R
```

Paths are resolved relative to each script. Missing packages are installed on first run.

- **Empirical analysis:** tidyverse, patchwork, openxlsx, ape, phytools, ggpubr, plus Bioconductor treeio and ggtree.
- **Simulation analysis:** tidyverse, patchwork, scales. It also reads the empirical FastTree DP LMM SNP-GWAS results.

Figures are saved as 300 dpi TIFF files with LZW compression.

| Manuscript item | Output file | Script |
|---|---|---|
| Figure 1 | — (pipeline schematic, not code-generated) | — |
| Figure 2 | `empirical_analysis/outputs/Figure_2_QQplots_population_structure_methods.tiff` | empirical |
| Table 1 | `empirical_analysis/outputs/Table_1_N_significant_hits_per_method.csv` | empirical |
| Figure 3 | `simulation_analysis/outputs/Figure_3_Power_by_EOD_LOD_ratio_and_AF_set1.tiff` | simulation |
| Figure 4 | `simulation_analysis/outputs/Figure_4_Power_by_HI_set2.tiff` | simulation |
| Figure 5 | `simulation_analysis/outputs/Figure_5_Power_by_HI_AF_EOD_LOD_ratio_set2.tiff` | simulation |
| Table 2 | `simulation_analysis/outputs/Table_2_Power_low_vs_high_HI.csv` | simulation |
| Figure S1 | `simulation_analysis/outputs/Figure_S1_AF_vs_HI.tiff` | simulation |
| Figure S2 | `simulation_analysis/outputs/Figure_S2_Crude_OR_by_AF.tiff` | simulation |
| Figure S3 | `empirical_analysis/outputs/Figure_S3_Phylogenies.tiff` | empirical |
| Table S1 | `empirical_analysis/outputs/Table_S1_Recombination_rates_per_SC.csv` | empirical |
| Figure S4 | `empirical_analysis/outputs/Figure_S4_Recombination_phylogenies.tiff` | empirical |
| Table S2 | `empirical_analysis/outputs/Table_S2_GWAS_hits_with_vs_without_recombination.csv` | empirical |
| Figure S5 | `simulation_analysis/outputs/Figure_S5_Power_by_AF_effect_size_set1.tiff` | simulation |
| Figure S6 | `simulation_analysis/outputs/Figure_S6_Power_by_HI_set1.tiff` | simulation |
| Figure S7 | `simulation_analysis/outputs/Figure_S7_AF_by_HI_set1.tiff` | simulation |
| Figure S8 | `simulation_analysis/outputs/Figure_S8_Power_by_AF_effect_size_HI_set1.tiff` | simulation |
| Figure S9 | `simulation_analysis/outputs/Figure_S9_Null_simulation_checks.tiff` | simulation |
| Figure S10 | `simulation_analysis/outputs/Figure_S10_QQplots_null_simulations.tiff` | simulation |
| Supplementary data | `empirical_analysis/outputs/Supplementary_Data_pyseer_significant_hits.xlsx` (all significant pyseer hits per method and pangenome tool) | empirical |

## Population structure correction methods

| Folder (`inputs/pyseer/population_structure_methods/`) | Method | Model |
|---|---|---|
| `CC` | Clonal complexes | FEM |
| `SC` | PopPUNK sequence clusters | FEM |
| `Mash_distances` | Mash distances | FEM |
| `SNP_kinship`, `Unitig_kinship`, `COG_kinship` | Kinship matrices from SNPs, unitigs or COGs | LMM |
| `RAxML`, `IQTREE` | SNP-based phylogenies | LMM |
| `FastTree_SP` | Whole-genome phylogeny, FastTree single precision | LMM |
| `FastTree_DP_FEM`, `FastTree_DP_LMM` | Whole-genome phylogeny, FastTree double precision | FEM / LMM |
| `FastTree_DP_LMM_noRecombination` | Recombination-free (Gubbins) FastTree DP phylogeny | LMM |
