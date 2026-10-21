#!/bin/bash
# =============================================================================
#  GBS-GWAS NULL Simulation Analysis — Pure Bash + AWK
#
#  Companion to run_simulation.sh / run_simulation_homoplasy.sh (power sims)
#  and a drop-in input for run_pyseer_simulation.sh (unmodified except for
#  the SIMDIR= line) / its HI-aware equivalent.
#
#  Per coauthor feedback: for null simulations, the EOD:LOD ratio among
#  reference-allele carriers should equal the EOD:LOD ratio among
#  alternative-allele carriers, i.e. the TRUE observed ratio in this dataset
#  (623:917 = 0.69:1). This is achieved by taking the exact same code/paths
#  as the power simulations, but replacing the ratio-based phenotype
#  assignment with a random permutation of the TRUE observed phenotype
#  vector (equivalent to `PHENOTYPE = sample(PHENOTYPE)` in R). Permuting
#  labels breaks any genotype-phenotype association while preserving the
#  global 623/917 totals, so carriers and non-carriers get the same
#  expected EOD:LOD ratio by construction — that is the null model.
#
#  DESIGN: output layout is a byte-for-byte match to the power-sim SIMDIRs
#  so run_pyseer_simulation.sh (or the HI-aware equivalent) runs UNCHANGED
#  against it — only the ratio_EODx_LODy folders are replaced with a single
#  ratio_null folder per SNP (no EOD:LOD ratio argument for the null model).
#
#    Set #1 (AF-based, mirrors run_simulation.sh):
#      <OUTDIR1>/02_selected_snps/selected_snps.tsv   (copied from --set1-snps)
#      <OUTDIR1>/03_phenotype_files/<freq_bin>/<VID>/ratio_null/sim_NNN.pheno
#
#    Set #2 (HI-based, mirrors run_simulation_homoplasy.sh):
#      <OUTDIR2>/02_selected_snps/selected_snps.tsv   (copied from --set2-snps)
#      <OUTDIR2>/03_phenotype_files/<freq_bin>/<hi_bin>/<VID>/ratio_null/sim_NNN.pheno
#
#  For each selected SNP, N_SIMS independent random permutations of the true
#  phenotype vector are written. Run these through the SAME pyseer /
#  correction-method code used for the power simulations. The proportion of
#  a SNP's null sims with LRT p-value below the Bonferroni threshold
#  (2.55e-6) is its empirical false-positive rate.
# =============================================================================

#SBATCH --job-name=gwas_null_sim
#SBATCH --partition=fasse
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/gwas_null_sim_%j.out
#SBATCH --error=logs/gwas_null_sim_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rbalasubramanian@g.harvard.edu

set -euo pipefail

# =============================================================================
# ████  USER INPUTS  ████
# =============================================================================
VCF_FILE="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/\
Nextflow_pipeline_EOD_LOD_VLOD_CDC/11.SNIPPY_MULTI/merged.vcf"

# True (observed) EOD/LOD phenotype, NOT a simulated ratio-based one.
# Must have a header row and, for every sample in the VCF, a sample-ID
# column and a binary 0/1 EOD(=0)/LOD-VLOD(=1) column.
# e.g. /n/hanage_gbs_gwas_l3/Lab/.../0.Raw_data/phenotypes_CDC.txt
PHENO_FILE="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/\
Nextflow_pipeline_EOD_LOD_VLOD_CDC/Files/phenotypes.txt"
PHENO_SAMPLE_COL="samples"
PHENO_BINARY_COL="binary"

# Point these at your FINAL, production selected_snps.tsv files — the exact
# ones consumed by run_simulation.sh (Set #1) / run_simulation_homoplasy.sh
# (Set #2) for the power simulations, e.g.:
#   gwas_simulations_20260511/02_selected_snps/selected_snps.tsv
#   gwas_simulations_hi_20260511/02_selected_snps/selected_snps.tsv
# Leave a path blank ("") to skip that set.
SET1_SNPS="/n/hanage_gbs_gwas_l3/Lab/Variant_Simulation/gwas_simulations_20260520/02_selected_snps/selected_snps.tsv"
SET2_SNPS="/n/hanage_gbs_gwas_l3/Lab/Variant_Simulation/gwas_simulations_hi_20260511/02_selected_snps/selected_snps.tsv"

OUTDIR1="gwas_null_simulations_$(date +%Y%m%d)"       # Set #1 output (AF-based)
OUTDIR2="gwas_null_simulations_hi_$(date +%Y%m%d)"    # Set #2 output (HI-based)

TOT_EOD=623
TOT_LOD=917
N_SIMS=20
SEED=2024
FMT="tab"

# =============================================================================
# ████  STOP EDITING  ████
# =============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)                VCF_FILE="$2";        shift 2 ;;
        --pheno-file)         PHENO_FILE="$2";       shift 2 ;;
        --pheno-sample-col)   PHENO_SAMPLE_COL="$2"; shift 2 ;;
        --pheno-binary-col)   PHENO_BINARY_COL="$2"; shift 2 ;;
        --set1-snps)          SET1_SNPS="$2";        shift 2 ;;
        --set2-snps)          SET2_SNPS="$2";        shift 2 ;;
        --outdir1)             OUTDIR1="$2";          shift 2 ;;
        --outdir2)             OUTDIR2="$2";          shift 2 ;;
        --eod)                 TOT_EOD="$2";          shift 2 ;;
        --lod)                 TOT_LOD="$2";          shift 2 ;;
        --nsims)                N_SIMS="$2";           shift 2 ;;
        --seed)                 SEED="$2";             shift 2 ;;
        --fmt)                  FMT="$2";              shift 2 ;;
        -h|--help)              head -40 "$0";         exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

TOTAL=$((TOT_EOD + TOT_LOD))

AWK=awk
command -v gawk &>/dev/null && AWK=gawk

mkdir -p logs "00_true_phenotype"

LOG="logs/run_null_$(date +%Y%m%d_%H%M%S).log"

log()  { printf "[%s]  %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR:   $*"; exit 1; }

log "======================================================================"
log "  GBS-GWAS NULL Simulation Analysis  (Pure Bash + AWK)"
log "======================================================================"
log "  VCF file         : $VCF_FILE"
log "  Phenotype file   : $PHENO_FILE  (cols: $PHENO_SAMPLE_COL / $PHENO_BINARY_COL)"
log "  Set #1 SNP list  : ${SET1_SNPS:-<skipped>}"
log "  Set #1 outdir    : $OUTDIR1"
log "  Set #2 SNP list  : ${SET2_SNPS:-<skipped>}"
log "  Set #2 outdir    : $OUTDIR2"
log "  EOD (0) / LOD(1) : $TOT_EOD / $TOT_LOD  (total n=$TOTAL)"
log "  Null sims/SNP    : $N_SIMS  (independent random permutations of the true"
log "                     phenotype vector — no EOD:LOD ratio argument; the"
log "                     null ratio among carriers = overall ratio by design)"
log "  Seed             : $SEED"
log "  AWK binary       : $AWK"
log "======================================================================"

[[ -f "$VCF_FILE" ]]   || die "VCF not found: $VCF_FILE"
[[ -f "$PHENO_FILE" ]] || die "Phenotype file not found: $PHENO_FILE"
[[ -n "$SET1_SNPS" || -n "$SET2_SNPS" ]] || die "Provide at least one of --set1-snps / --set2-snps"
[[ -z "$SET1_SNPS" || -f "$SET1_SNPS" ]] || die "Set #1 SNP list not found: $SET1_SNPS"
[[ -z "$SET2_SNPS" || -f "$SET2_SNPS" ]] || die "Set #2 SNP list not found: $SET2_SNPS"

N_VCF_SAMPLES=$(grep -m1 "^#CHROM" "$VCF_FILE" \
    | tr '\t' '\n' | tail -n +10 | wc -l | tr -d ' ')
log "VCF sample count  : $N_VCF_SAMPLES  (expected $TOTAL)"
[[ "$N_VCF_SAMPLES" -eq "$TOTAL" ]] \
    || warn "Sample count mismatch — got $N_VCF_SAMPLES expected $TOTAL"

# =============================================================================
# STEP 1 — EXTRACT SAMPLE NAMES (same order pyseer/VCF uses)
# =============================================================================
log ""
log "── Step 1: Extracting sample names ──────────────────────────────────"

SAMPLES_FILE="00_true_phenotype/samples.txt"

grep -m1 "^#CHROM" "$VCF_FILE" \
    | tr '\t' '\n' \
    | tail -n +10 \
    > "$SAMPLES_FILE"

N_SAMP=$(wc -l < "$SAMPLES_FILE" | tr -d ' ')
log "  $N_SAMP sample names written → 00_true_phenotype/samples.txt"

# =============================================================================
# STEP 2 — BUILD THE TRUE PHENOTYPE VECTOR (in VCF sample order)
#
# This is the REAL observed EOD/LOD status per isolate — not a simulated
# ratio-based assignment. It is what gets permuted for the null model.
# Shared by both sets, built once.
# =============================================================================
log ""
log "── Step 2: Building true phenotype vector ───────────────────────────"

TRUE_PHENO_FILE="00_true_phenotype/true_phenotype.tsv"

$AWK -F'\t' \
    -v samples_file="$SAMPLES_FILE" \
    -v sample_col="$PHENO_SAMPLE_COL" \
    -v binary_col="$PHENO_BINARY_COL" \
    -v out_file="$TRUE_PHENO_FILE" \
    -v tot_eod="$TOT_EOD" \
    -v tot_lod="$TOT_LOD" \
'
BEGIN {
    ns = 0
    while ((getline s < samples_file) > 0) { samples[ns] = s; ns++ }
    close(samples_file)
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == sample_col) sc = i
        if ($i == binary_col) bc = i
    }
    if (!sc || !bc) {
        printf "ERROR: could not find columns \"%s\"/\"%s\" in phenotype file header\n",
            sample_col, binary_col > "/dev/stderr"
        exit 1
    }
    next
}
{
    pheno_lookup[$sc] = $bc
}
END {
    print "sample\tphenotype" > out_file
    n_eod = 0; n_lod = 0; n_missing = 0
    for (j = 0; j < ns; j++) {
        s = samples[j]
        if (!(s in pheno_lookup)) {
            printf "MISSING phenotype for sample: %s\n", s > "/dev/stderr"
            n_missing++
            continue
        }
        p = pheno_lookup[s] + 0
        print s "\t" p >> out_file
        if (p == 0) n_eod++; else n_lod++
    }
    close(out_file)
    printf "n_eod=%d n_lod=%d n_missing=%d\n", n_eod, n_lod, n_missing > "/dev/stderr"
    if (n_missing > 0) exit 1
    if (n_eod != tot_eod || n_lod != tot_lod) {
        printf "ASSERTION FAILED: expected EOD=%d LOD=%d, got EOD=%d LOD=%d\n",
            tot_eod, tot_lod, n_eod, n_lod > "/dev/stderr"
        exit 1
    }
}
' "$PHENO_FILE" 2> "logs/step2_stderr.log" \
    || die "True phenotype construction failed — see logs/step2_stderr.log"

cat "logs/step2_stderr.log" | tee -a "$LOG" > /dev/null
log "  True phenotype vector → 00_true_phenotype/true_phenotype.tsv"
log "  Confirmed EOD=$TOT_EOD / LOD=$TOT_LOD across all $N_SAMP VCF samples"

# =============================================================================
# Shared AWK permutation writer.
# Given a fixed set of samples (in VCF order) with true phenotype labels,
# writes N_SIMS files, each a fresh Fisher-Yates shuffle of the labels —
# exactly equivalent to `PHENOTYPE = sample(PHENOTYPE)` in R.
# =============================================================================
write_null_pheno_files () {
    local SNP_DIR="$1"
    local COMBO_SEED="$2"

    mkdir -p "$SNP_DIR"

    $AWK \
        -v true_pheno_file="$TRUE_PHENO_FILE" \
        -v n_sims="$N_SIMS" \
        -v combo_seed="$COMBO_SEED" \
        -v snp_dir="$SNP_DIR" \
        -v fmt="$FMT" \
        -v tot_eod="$TOT_EOD" \
        -v tot_lod="$TOT_LOD" \
    '
    BEGIN {
        ns = 0
        getline hdr < true_pheno_file
        while ((getline line < true_pheno_file) > 0) {
            split(line, f, "\t")
            samples[ns] = f[1]
            true_pheno[ns] = f[2] + 0
            ns++
        }
        close(true_pheno_file)

        for (sim = 1; sim <= n_sims; sim++) {

            perm_seed = (combo_seed + sim * 999983) % 16777216
            srand(perm_seed)

            for (j = 0; j < ns; j++) pheno[j] = true_pheno[j]
            for (j = ns - 1; j > 0; j--) {
                k = int(rand() * (j + 1))
                tmp = pheno[j]; pheno[j] = pheno[k]; pheno[k] = tmp
            }

            eod_n = 0; lod_n = 0
            for (j = 0; j < ns; j++) {
                if (pheno[j] == 0) eod_n++; else lod_n++
            }
            if (eod_n != tot_eod || lod_n != tot_lod) {
                printf "ASSERTION FAILED sim=%d EOD=%d LOD=%d\n",
                    sim, eod_n, lod_n > "/dev/stderr"
                exit 1
            }

            outfile = snp_dir "/" sprintf("sim_%03d.pheno", sim)
            if (fmt == "tab") {
                print "sample\tphenotype" > outfile
                for (j = 0; j < ns; j++)
                    print samples[j] "\t" pheno[j] > outfile
            } else {
                print "FID\tIID\tPHENOTYPE" > outfile
                for (j = 0; j < ns; j++)
                    print samples[j] "\t" samples[j] "\t" (pheno[j]+1) > outfile
            }
            close(outfile)
        }
    }
    ' /dev/null
}

# =============================================================================
# STEP 3a — SET #1 (AF-based) — mirrors run_simulation.sh layout exactly
#
# selected_snps.tsv columns (10):
#   target_pct  snp_rank  variant_id  chrom  pos  ref  alt
#   n_carriers  actual_pct  delta_pct
#
# Output:
#   <OUTDIR1>/02_selected_snps/selected_snps.tsv        (copied verbatim)
#   <OUTDIR1>/03_phenotype_files/<freq_bin>/<VID>/ratio_null/sim_NNN.pheno
#
# This is a straight drop-in for run_pyseer_simulation.sh: point SIMDIR at
# <OUTDIR1> and run it unmodified. Its Step 1 (per-SNP SNP.vcf creation)
# reads the same 10-column selected_snps.tsv and writes SNP.vcf into
# <freq_bin>/<VID>/ — one level above ratio_null/, exactly as it expects.
# =============================================================================
if [[ -n "$SET1_SNPS" ]]; then
    log ""
    log "── Step 3a [Set #1 / AF-based]: Generating null phenotype files ────"

    mkdir -p "$OUTDIR1/logs" "$OUTDIR1/02_selected_snps" "$OUTDIR1/03_phenotype_files"
    cp "$SET1_SNPS" "$OUTDIR1/02_selected_snps/selected_snps.tsv"

    N_ROWS=$(( $(wc -l < "$SET1_SNPS") - 1 ))
    log "  SNP list          : $SET1_SNPS"
    log "  SNPs              : $N_ROWS"
    log "  Files to write    : $N_ROWS × $N_SIMS = $(( N_ROWS * N_SIMS ))"

    while IFS=$'\t' read -r TARGET RANK VID CHROM POS REF ALT N_CARR ACT DELTA; do

        FREQ_LABEL=$(printf "freq_%02dpct" "$TARGET")
        SNP_DIR="$OUTDIR1/03_phenotype_files/$FREQ_LABEL/$VID/ratio_null"

        COMBO_HASH=$(printf "set1|%s" "$VID" | cksum | cut -d' ' -f1)
        COMBO_SEED=$(( (SEED + COMBO_HASH) % 16777216 ))

        write_null_pheno_files "$SNP_DIR" "$COMBO_SEED"

    done < <(tail -n +2 "$SET1_SNPS")

    N1=$(find "$OUTDIR1/03_phenotype_files" -name "*.pheno" 2>/dev/null | wc -l | tr -d ' ')
    log "  ✓  Set #1 done — $N1 null phenotype files → $OUTDIR1/03_phenotype_files/"
fi

# =============================================================================
# STEP 3b — SET #2 (HI-based) — mirrors run_simulation_homoplasy.sh layout
#
# selected_snps.tsv columns (12):
#   target_pct  hi_bin  snp_rank  variant_id  chrom  pos  ref  alt
#   n_carriers  actual_pct  delta_pct  actual_hi
#
# Output:
#   <OUTDIR2>/02_selected_snps/selected_snps.tsv                  (copied verbatim)
#   <OUTDIR2>/03_phenotype_files/<freq_bin>/<hi_bin>/<VID>/ratio_null/sim_NNN.pheno
#
# NOTE: run_pyseer_simulation.sh as written only parses the 4-level Set #1
# layout (freq_bin/VID/ratio/sim); it does not handle the extra hi_bin
# level here. Point whichever HI-aware pyseer script your collaborator
# uses for the power-sim HI runs at <OUTDIR2> instead — this layout matches
# run_simulation_homoplasy.sh's 03_phenotype_files/ output exactly.
# =============================================================================
if [[ -n "$SET2_SNPS" ]]; then
    log ""
    log "── Step 3b [Set #2 / HI-based]: Generating null phenotype files ────"

    mkdir -p "$OUTDIR2/logs" "$OUTDIR2/02_selected_snps" "$OUTDIR2/03_phenotype_files"
    cp "$SET2_SNPS" "$OUTDIR2/02_selected_snps/selected_snps.tsv"

    N_ROWS=$(( $(wc -l < "$SET2_SNPS") - 1 ))
    log "  SNP list          : $SET2_SNPS"
    log "  SNPs              : $N_ROWS"
    log "  Files to write    : $N_ROWS × $N_SIMS = $(( N_ROWS * N_SIMS ))"

    while IFS=$'\t' read -r TARGET HI_BIN RANK VID CHROM POS REF ALT N_CARR ACT DELTA ACTUAL_HI; do

        FREQ_LABEL=$(printf "freq_%02dpct" "$TARGET")
        SNP_DIR="$OUTDIR2/03_phenotype_files/$FREQ_LABEL/$HI_BIN/$VID/ratio_null"

        COMBO_HASH=$(printf "set2|%s|%s" "$VID" "$HI_BIN" | cksum | cut -d' ' -f1)
        COMBO_SEED=$(( (SEED + COMBO_HASH) % 16777216 ))

        write_null_pheno_files "$SNP_DIR" "$COMBO_SEED"

    done < <(tail -n +2 "$SET2_SNPS")

    N2=$(find "$OUTDIR2/03_phenotype_files" -name "*.pheno" 2>/dev/null | wc -l | tr -d ' ')
    log "  ✓  Set #2 done — $N2 null phenotype files → $OUTDIR2/03_phenotype_files/"
fi

# =============================================================================
# STEP 4 — SPOT CHECK
# =============================================================================
log ""
log "── Step 4: Spot-checking 10 random null phenotype files per set ────"

PASS=0; FAIL=0
for D in "$OUTDIR1" "$OUTDIR2"; do
    [[ -d "$D/03_phenotype_files" ]] || continue
    while IFS= read -r PFILE; do
        [[ -f "$PFILE" ]] || { log "  ✗  MISSING: $PFILE"; FAIL=$((FAIL+1)); continue; }
        READ_EOD=$(awk -F'\t' 'NR>1 && $2==0{c++} END{print c+0}' "$PFILE")
        READ_LOD=$(awk -F'\t' 'NR>1 && $2==1{c++} END{print c+0}' "$PFILE")
        LABEL="${PFILE#$D/}"
        if [[ "$READ_EOD" -eq "$TOT_EOD" && "$READ_LOD" -eq "$TOT_LOD" ]]; then
            log "  ✓  $LABEL  EOD=$READ_EOD LOD=$READ_LOD"
            PASS=$((PASS+1))
        else
            log "  ✗  $LABEL  GOT EOD=$READ_EOD LOD=$READ_LOD  WANT $TOT_EOD/$TOT_LOD"
            FAIL=$((FAIL+1))
        fi
    done < <(find "$D/03_phenotype_files" -name "*.pheno" 2>/dev/null | sort -R | head -10)
done

# Two random files for the SAME SNP should differ (permutation actually happened)
for D in "$OUTDIR1" "$OUTDIR2"; do
    [[ -d "$D/03_phenotype_files" ]] || continue
    DUP_CHECK_DIR=$(find "$D/03_phenotype_files" -type d -name "ratio_null" 2>/dev/null | head -1)
    if [[ -n "$DUP_CHECK_DIR" ]]; then
        F1="$DUP_CHECK_DIR/sim_001.pheno"
        F2="$DUP_CHECK_DIR/sim_002.pheno"
        if [[ -f "$F1" && -f "$F2" ]]; then
            if diff -q "$F1" "$F2" > /dev/null 2>&1; then
                warn "  sim_001 and sim_002 are IDENTICAL in $DUP_CHECK_DIR — permutation may not be working"
            else
                log "  ✓  sim_001 != sim_002 in $DUP_CHECK_DIR (permutation confirmed to vary)"
            fi
        fi
    fi
done

# =============================================================================
# FINAL SUMMARY
# =============================================================================
FINAL_N1=0; FINAL_N2=0
[[ -d "$OUTDIR1/03_phenotype_files" ]] && FINAL_N1=$(find "$OUTDIR1/03_phenotype_files" -name "*.pheno" 2>/dev/null | wc -l | tr -d ' ')
[[ -d "$OUTDIR2/03_phenotype_files" ]] && FINAL_N2=$(find "$OUTDIR2/03_phenotype_files" -name "*.pheno" 2>/dev/null | wc -l | tr -d ' ')

log ""
log "======================================================================"
log "  FINAL SUMMARY"
log "======================================================================"
log "  Set #1 null phenotype files : $FINAL_N1   ($OUTDIR1)"
log "  Set #2 null phenotype files : $FINAL_N2   ($OUTDIR2)"
log "  Spot check                  : $PASS passed  /  $FAIL failed"
log ""
log "  Next step:"
log "    Set #1 → sbatch-edit run_pyseer_simulation.sh: SIMDIR=\"$OUTDIR1\", resubmit unchanged."
log "    Set #2 → point your HI-aware pyseer runner at SIMDIR=\"$OUTDIR2\" the same way."
log "  Each SNP's summary rows will show ratio == \"ratio_null\" — the proportion"
log "  of a SNP's $N_SIMS null sims with LRT p-value < 2.55e-6 is its"
log "  empirical false-positive rate under that correction method."

if [[ "$FAIL" -eq 0 && "$(( FINAL_N1 + FINAL_N2 ))" -gt 0 ]]; then
    log ""
    log "  ✅  All checks passed — output is ready"
else
    [[ "$(( FINAL_N1 + FINAL_N2 ))" -eq 0 ]] && log "  ❌  No files written — check $LOG"
    [[ "$FAIL" -gt 0 ]] && log "  ❌  Spot-check failures — check $LOG"
fi

log ""
log "======================================================================"
log "  DONE  $(date)"
log "======================================================================"

