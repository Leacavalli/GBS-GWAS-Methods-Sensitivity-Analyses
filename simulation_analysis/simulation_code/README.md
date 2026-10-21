# Simulation code

Cluster (SLURM) scripts that produce the simulated phenotypes and homoplasy index used by `simulation_analysis/analysis_code/simulation_analysis.R`.
Paths at the top of each script point to the Hanage lab FASRC storage (`/n/hanage_gbs_gwas_l3/...`). Edit them before running elsewhere.

| File | Purpose | Main output |
|---|---|---|
| `homoplasy.sh` | Builds a conda env (`homoplasy_env`, Python 3.10). Runs an embedded Python script that reads the merged Snippy VCF and the FastTree DP phylogeny. It computes, for every SNP, the Fitch parsimony score, the minimum number of changes, the consistency index (CI) and the homoplasy index (HI = 1 − CI). | `homoplasy_results.tsv`, saved in this repo as `simulation_analysis/inputs/homoplasy_results_SNPs.tsv` |
| `run_simulation.sh` | **Set 1** power simulations, pure bash/awk. Selects 10 SNPs per allele-frequency bin, from 5% to 95% in 5% steps. The 5% bin is restricted to 6–10% and the 95% bin to 90–94%, so that SNPs pass pyseer's `--min-af`/`--max-af` 0.05/0.95 filters. For each SNP it generates 20 phenotype files at each of 4 EOD:LOD ratios among carriers (1:3, 1:2, 2:1, 3:1), keeping the global totals fixed at 623 EOD and 917 LOD. | `03_phenotype_files/freq_XXpct/<SNP>/ratio_EODx_LODy/sim_NNN.pheno` |
| `run_simulation_homoplasy.sh` | **Set 2** power simulations. Same phenotype generator as above, but it randomly samples 3 SNPs from each allele-frequency bin (5, 10, 15, 20, 25 and 40%) × HI bin (0–0.10, 0.45–0.55, 0.60–0.70, 0.70–0.80, 0.80–0.90, 0.90–1.00) combination. Reads the HI table produced by `homoplasy.sh`. | `03_phenotype_files/freq_XXpct/hi_XX_XX/<SNP>/ratio_EODx_LODy/sim_NNN.pheno` |
| `code from other github - CHECK WITH RUCHITA/run_simulation.sh` | Earlier version of `run_simulation.sh`. It selects 1 SNP per bin over 5–50% only (800 phenotype files), and its SLURM header is commented out. | — |
| `code from other github - CHECK WITH RUCHITA/validate_simulations.sh` | QC of a simulation output directory: file counts, file structure, sample names vs VCF, fixed EOD/LOD totals, recomputed allele frequencies, carrier EOD:LOD ratios, non-carrier totals and variability between simulations. | `<sim_dir>/validation/validation.log`, `errors.txt` |

## To check with Ruchita

1. **The "other github" folder.** Is the older `run_simulation.sh` still needed, or can it be deleted? It is superseded by `run_simulation.sh` here.
2. **`validate_simulations.sh` uses the old design.** It hard-codes `EXP_FILES=800`, which is 10 SNPs × 4 ratios × 20 simulations. The current Set 1 design is 19 bins × 10 SNPs × 4 × 20 = 15,200 files, and Set 2 differs again. Was it re-run on the final simulations? If so, should it be updated and moved next to the other scripts?
3. **Missing steps.** No script in this repo covers:
   - running pyseer on the simulated phenotypes for the 13 correction methods;
   - generating the permuted (null) phenotypes;
   - compiling the `sim_GWAS_results_set*.txt` tables in `simulation_analysis/inputs/`.

   Should these scripts be added, or should we point to where they live?
4. **Personal details in SBATCH headers.** `run_simulation.sh` and `run_simulation_homoplasy.sh` include `--mail-user=rbalasubramanian@g.harvard.edu`. Keep it, or replace it with a placeholder before the repo is made public?
