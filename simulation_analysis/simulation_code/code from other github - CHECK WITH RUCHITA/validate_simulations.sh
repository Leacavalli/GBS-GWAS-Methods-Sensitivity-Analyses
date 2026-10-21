#!/bin/bash
# =============================================================================
#  validate_simulations.sh
#  Comprehensive QC for GBS-GWAS simulation output
#
#  CHECKS PERFORMED
#  ────────────────
#  1. File count              → exactly 800 .pheno files exist
#  2. File structure          → correct number of rows and columns per file
#  3. Sample names            → match VCF exactly, no duplicates
#  4. Global EOD/LOD totals   → EOD=623 LOD=917 in EVERY file
#  5. Allele frequency        → re-compute from VCF and compare to reported
#  6. Carrier EOD/LOD ratio   → check ratio among SNP carriers in every file
#  7. Non-carrier totals      → correct EOD/LOD fill in non-carriers
#  8. Simulation variability  → confirm 20 sims differ from each other
#
#  USAGE
#  ─────
#  bash validate_simulations.sh <sim_output_dir> [path/to/merged.vcf]
#
#  Examples:
#    bash validate_simulations.sh gwas_simulations_20250408
#    bash validate_simulations.sh gwas_simulations_20250408 /path/to/merged.vcf
#
# =============================================================================

set -euo pipefail

# =============================================================================
# INPUTS
# =============================================================================
SIMDIR="${1:-}"
VCF_OVERRIDE="${2:-}"

[[ -z "$SIMDIR" ]] && {
    echo "Usage: bash validate_simulations.sh <sim_output_dir> [path/to/merged.vcf]"
    echo ""
    echo "Examples:"
    echo "  bash validate_simulations.sh gwas_simulations_20250408"
    echo "  bash validate_simulations.sh gwas_simulations_20250408 \\"
    echo "    /n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/\\"
    echo "    Nextflow_pipeline_EOD_LOD_VLOD_CDC/11.SNIPPY_MULTI/merged.vcf"
    exit 1
}

[[ -d "$SIMDIR" ]] || { echo "ERROR: Directory not found: $SIMDIR"; exit 1; }

# =============================================================================
# FIXED EXPECTED VALUES
# =============================================================================
EXP_EOD=623
EXP_LOD=917
EXP_TOTAL=$(( EXP_EOD + EXP_LOD ))    # 1540
EXP_FILES=800                           # 10 SNPs × 4 ratios × 20 sims
EXP_SIMS=20                             # simulations per (SNP × ratio)

# prefer gawk for better arithmetic
AWK=awk
command -v gawk &>/dev/null && AWK=gawk

# =============================================================================
# DIRECTORIES AND LOGGING
# =============================================================================
VAL_DIR="$SIMDIR/validation"
mkdir -p "$VAL_DIR"

LOG="$VAL_DIR/validation.log"
ERRORS="$VAL_DIR/errors.txt"

> "$LOG"
> "$ERRORS"

# ── logging helpers ────────────────────────────────────────────────────────────
log()     { printf "[%s]  %s\n"  "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG"; }
pass()    { printf "  ✓  %s\n"   "$*"                       | tee -a "$LOG"; }
fail()    { printf "  ✗  %s\n"   "$*"                       | tee -a "$LOG"
            echo "$*" >> "$ERRORS"; }
warn()    { printf "  ⚠  %s\n"   "$*"                       | tee -a "$LOG"; }
section() { printf "\n── %s ──\n" "$*"                      | tee -a "$LOG"; }

TOTAL_PASS=0
TOTAL_FAIL=0
ok()  { TOTAL_PASS=$(( TOTAL_PASS + 1 )); pass "$*"; }
err() { TOTAL_FAIL=$(( TOTAL_FAIL + 1 )); fail "$*"; }

log "======================================================================"
log "  GBS-GWAS Simulation Validation"
log "  Directory : $SIMDIR"
log "  Expected  : EOD=$EXP_EOD  LOD=$EXP_LOD  Total=$EXP_TOTAL"
log "  Expected  : $EXP_FILES phenotype files  ($EXP_SIMS sims per combo)"
log "======================================================================"

# =============================================================================
# LOCATE VCF
# =============================================================================
# Priority: 1) explicit CLI argument  2) extract from run log  3) skip CHECK 5
if [[ -n "$VCF_OVERRIDE" ]]; then
    VCF_PATH="$VCF_OVERRIDE"
    log "VCF path (from argument) : $VCF_PATH"
else
    # Use awk -F' : ' to avoid grep -oP lookbehind (not supported on Rocky8)
    VCF_PATH=$(awk -F' : ' '/VCF[[:space:]]/{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' \
               "$SIMDIR/logs/run.log" 2>/dev/null || true)
    log "VCF path (from log)      : ${VCF_PATH:-NOT FOUND}"
fi

SKIP_FREQ_CHECK=0
if [[ -z "$VCF_PATH" || ! -f "$VCF_PATH" ]]; then
    warn "VCF file not found — CHECK 5 (allele frequency) will be skipped"
    warn "To enable CHECK 5 re-run with the VCF as a second argument:"
    warn "  bash validate_simulations.sh $SIMDIR /path/to/merged.vcf"
    SKIP_FREQ_CHECK=1
else
    log "VCF confirmed            : $VCF_PATH"
fi

# =============================================================================
# CHECK 1 — REQUIRED OUTPUT FILES AND DIRECTORIES EXIST
# =============================================================================
section "CHECK 1: Required output files and directories exist"

for REQUIRED in \
    "$SIMDIR/01_allele_frequencies/all_variant_frequencies.tsv" \
    "$SIMDIR/02_selected_snps/selected_snps.tsv" \
    "$SIMDIR/03_phenotype_files" \
    "$SIMDIR/tmp/samples.txt" \
    "$SIMDIR/logs/run.log"; do

    if [[ -f "$REQUIRED" || -d "$REQUIRED" ]]; then
        ok "Found : $(basename $REQUIRED)"
    else
        err "MISSING : $REQUIRED"
    fi
done

# =============================================================================
# CHECK 2 — TOTAL PHENOTYPE FILE COUNT + BREAKDOWN
# =============================================================================
section "CHECK 2: Phenotype file count"

N_FILES=$(find "$SIMDIR/03_phenotype_files" -name "*.pheno" 2>/dev/null \
          | wc -l | tr -d ' ')

log "  Found    : $N_FILES .pheno files"
log "  Expected : $EXP_FILES"

if [[ "$N_FILES" -eq "$EXP_FILES" ]]; then
    ok "File count correct: $N_FILES"
else
    err "File count WRONG: got $N_FILES expected $EXP_FILES"
fi

# Breakdown: freq_bin × ratio → should each have exactly 20 files
log ""
log "  Breakdown by frequency bin and ratio:"
printf "  %-18s  %-16s  %s\n" "freq_bin" "ratio" "n_files" | tee -a "$LOG"
printf "  %-18s  %-16s  %s\n" "────────" "─────" "───────" | tee -a "$LOG"

N_COMBO_FAIL=0
while IFS= read -r FBIN_DIR; do
    FBNAME=$(basename "$FBIN_DIR")
    while IFS= read -r RDIR; do
        RNAME=$(basename "$RDIR")
        N=$(find "$RDIR" -name "*.pheno" 2>/dev/null | wc -l | tr -d ' ')
        STATUS="✓"
        [[ "$N" -ne "$EXP_SIMS" ]] && { STATUS="✗ WRONG"; N_COMBO_FAIL=$(( N_COMBO_FAIL + 1 )); }
        printf "  %-18s  %-16s  %s %s\n" "$FBNAME" "$RNAME" "$N" "$STATUS" \
            | tee -a "$LOG"
    done < <(find "$FBIN_DIR" -maxdepth 2 -type d -name "ratio_*" 2>/dev/null | sort)
done < <(find "$SIMDIR/03_phenotype_files" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | sort)

if [[ "$N_COMBO_FAIL" -eq 0 ]]; then
    ok "All combinations have exactly $EXP_SIMS simulation files"
else
    err "$N_COMBO_FAIL combinations have wrong simulation counts"
fi

# =============================================================================
# CHECK 3 — SAMPLE NAMES MATCH VCF
# =============================================================================
section "CHECK 3: Sample names match expected list"

SAMPLES_FILE="$SIMDIR/tmp/samples.txt"
N_SAMP=$(wc -l < "$SAMPLES_FILE" | tr -d ' ')

log "  Samples in tmp/samples.txt : $N_SAMP  (expected $EXP_TOTAL)"

if [[ "$N_SAMP" -eq "$EXP_TOTAL" ]]; then
    ok "Sample count correct: $N_SAMP"
else
    err "Sample count WRONG: got $N_SAMP expected $EXP_TOTAL"
fi

# Compare sample column of first .pheno file against samples.txt
FIRST_PHENO=$(find "$SIMDIR/03_phenotype_files" -name "sim_001.pheno" \
              2>/dev/null | sort | head -1)

if [[ -n "$FIRST_PHENO" ]]; then

    awk 'NR>1{print $1}' "$FIRST_PHENO" > "$VAL_DIR/pheno_samples.txt"
    N_PHENO_SAMP=$(wc -l < "$VAL_DIR/pheno_samples.txt" | tr -d ' ')

    log "  Rows in test .pheno file : $N_PHENO_SAMP"

    DIFF=$(diff "$SAMPLES_FILE" "$VAL_DIR/pheno_samples.txt" | wc -l | tr -d ' ')
    if [[ "$DIFF" -eq 0 ]]; then
        ok "Sample names and order in .pheno files match VCF exactly"
    else
        err "Sample name MISMATCH: $DIFF diff lines between samples.txt and .pheno"
        diff "$SAMPLES_FILE" "$VAL_DIR/pheno_samples.txt" \
            | head -20 | sed 's/^/    /' | tee -a "$LOG"
    fi

    # Duplicate sample check
    N_UNIQ=$(sort "$VAL_DIR/pheno_samples.txt" | uniq | wc -l | tr -d ' ')
    if [[ "$N_UNIQ" -eq "$N_PHENO_SAMP" ]]; then
        ok "No duplicate sample names ($N_UNIQ unique)"
    else
        err "DUPLICATE sample names: $N_PHENO_SAMP rows but only $N_UNIQ unique names"
    fi

else
    warn "No sim_001.pheno files found — cannot check sample names"
fi

# =============================================================================
# CHECK 4 — GLOBAL EOD/LOD TOTALS IN EVERY FILE
# =============================================================================
section "CHECK 4: Global EOD/LOD totals in every phenotype file"
log "  Checking all $N_FILES files for EOD=$EXP_EOD and LOD=$EXP_LOD ..."
log "  (This may take a few minutes)"

TOTALS_REPORT="$VAL_DIR/global_totals_check.tsv"
printf "file\tn_rows\tn_eod\tn_lod\tn_other\tstatus\n" > "$TOTALS_REPORT"

N_TOT_PASS=0
N_TOT_FAIL=0
N_TOT_CHECKED=0

while IFS= read -r PFILE; do

    N_TOT_CHECKED=$(( N_TOT_CHECKED + 1 ))

    RESULT=$($AWK -F'\t' '
    NR == 1 { next }
    {
        n_rows++
        if      ($2 == 0) n_eod++
        else if ($2 == 1) n_lod++
        else              n_other++
    }
    END { printf "%d\t%d\t%d\t%d", n_rows+0, n_eod+0, n_lod+0, n_other+0 }
    ' "$PFILE")

    N_ROW=$(echo "$RESULT" | cut -f1)
    N_EOD=$(echo "$RESULT" | cut -f2)
    N_LOD=$(echo "$RESULT" | cut -f3)
    N_OTH=$(echo "$RESULT" | cut -f4)

    RELPATH="${PFILE#$SIMDIR/}"

    if [[ "$N_EOD" -eq "$EXP_EOD" && \
          "$N_LOD" -eq "$EXP_LOD" && \
          "$N_OTH" -eq 0 && \
          "$N_ROW" -eq "$EXP_TOTAL" ]]; then
        STATUS="PASS"
        N_TOT_PASS=$(( N_TOT_PASS + 1 ))
    else
        STATUS="FAIL"
        N_TOT_FAIL=$(( N_TOT_FAIL + 1 ))
        echo "GLOBAL_TOTAL_FAIL: $RELPATH rows=$N_ROW eod=$N_EOD lod=$N_LOD other=$N_OTH" \
            >> "$ERRORS"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$RELPATH" "$N_ROW" "$N_EOD" "$N_LOD" "$N_OTH" "$STATUS" \
        >> "$TOTALS_REPORT"

done < <(find "$SIMDIR/03_phenotype_files" -name "*.pheno" 2>/dev/null | sort)

log "  Checked  : $N_TOT_CHECKED files"
log "  Passed   : $N_TOT_PASS"
log "  Failed   : $N_TOT_FAIL"

if [[ "$N_TOT_FAIL" -eq 0 ]]; then
    ok "ALL $N_TOT_CHECKED files have correct global totals: EOD=$EXP_EOD LOD=$EXP_LOD"
else
    err "$N_TOT_FAIL / $N_TOT_CHECKED files have WRONG global totals"
    log "  First few failures:"
    grep "GLOBAL_TOTAL_FAIL" "$ERRORS" | head -5 | sed 's/^/    /' | tee -a "$LOG"
    log "  Full list: $TOTALS_REPORT"
fi

# =============================================================================
# CHECK 5 — ALLELE FREQUENCY RE-VERIFICATION
# =============================================================================
section "CHECK 5: Re-verify allele frequencies against VCF"

if [[ "$SKIP_FREQ_CHECK" -eq 1 ]]; then
    warn "SKIPPED — no VCF path available"
    warn "Re-run: bash validate_simulations.sh $SIMDIR /path/to/merged.vcf"
else
    log "  Re-computing allele frequency directly from VCF for each selected SNP..."
    log "  Comparing recomputed values to those stored in selected_snps.tsv"

    SEL_FILE="$SIMDIR/02_selected_snps/selected_snps.tsv"

    FREQ_RECHECK="$VAL_DIR/frequency_recheck.tsv"
    printf "target_pct\tvariant_id\treported_pct\trecomputed_pct\tdelta_pp\t\
n_carriers_reported\tn_carriers_recomputed\tstatus\n" \
        > "$FREQ_RECHECK"

    N_FREQ_PASS=0
    N_FREQ_FAIL=0

    # Read each selected SNP from the TSV (skip header)
    tail -n +2 "$SEL_FILE" | while IFS=$'\t' read -r \
        TARGET VID CHROM POS REF ALT N_CARR ACT DELTA; do

        log "  Checking ${TARGET}%  →  $VID  (reported=${ACT}%  n_carriers=${N_CARR})"

        # Re-compute from VCF using awk — no grep -oP needed
        RECOMP=$($AWK -F'\t' \
            -v chrom="$CHROM" \
            -v pos="$POS"     \
            -v ref="$REF"     \
            -v alt="$ALT"     \
        '
        /^##/ { next }
        /^#CHROM/ { n_samples = NF - 9; next }
        {
            # match on all four fields for exact identification
            if ($1 != chrom || $2 != pos || $4 != ref || $5 != alt) next

            # locate GT index in FORMAT
            gt_pos = 1
            n_fmt  = split($9, fmt_arr, ":")
            for (f = 1; f <= n_fmt; f++) {
                if (fmt_arr[f] == "GT") { gt_pos = f; break }
            }

            # count ALT carriers
            alt_count = 0
            for (i = 10; i <= NF; i++) {
                n_f = split($i, sf, ":")
                gt  = (n_f >= gt_pos) ? sf[gt_pos] : sf[1]
                gsub(/\|/, "/", gt)
                n_al = split(gt, al, "/")
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "1") { alt_count++; break }
                }
            }

            printf "%.6f\t%d\t%d", alt_count / n_samples * 100, alt_count, n_samples
            exit
        }
        ' "$VCF_PATH" 2>/dev/null || echo "NOTFOUND")

        if [[ "$RECOMP" == "NOTFOUND" || -z "$RECOMP" ]]; then
            err "  Could not locate $VID in VCF — check CHROM/POS/REF/ALT"
            printf "%s\t%s\t%s\tNOTFOUND\tNA\t%s\tNA\tFAIL\n" \
                "$TARGET" "$VID" "$ACT" "$N_CARR" >> "$FREQ_RECHECK"
            N_FREQ_FAIL=$(( N_FREQ_FAIL + 1 ))
            continue
        fi

        RECOMP_PCT=$(  echo "$RECOMP" | cut -f1)
        RECOMP_CARR=$( echo "$RECOMP" | cut -f2)
        RECOMP_TOTAL=$(echo "$RECOMP" | cut -f3)

        # Delta between reported and recomputed (use awk for floating point)
        DELTA_PP=$($AWK "BEGIN{
            d = $ACT - $RECOMP_PCT
            if (d < 0) d = -d
            printf \"%.6f\", d
        }")

        DIFF_CARR=$(( N_CARR - RECOMP_CARR ))
        [[ "$DIFF_CARR" -lt 0 ]] && DIFF_CARR=$(( -DIFF_CARR ))

        # Pass if delta < 0.001 pp and carrier count is exact
        FREQ_OK=$($AWK "BEGIN{print ($DELTA_PP < 0.001) ? \"yes\" : \"no\"}")

        if [[ "$FREQ_OK" == "yes" && "$DIFF_CARR" -eq 0 ]]; then
            STATUS="PASS"
            N_FREQ_PASS=$(( N_FREQ_PASS + 1 ))
            log "  ✓  ${TARGET}%  $VID"
            log "       reported   = ${ACT}%   n_carriers = ${N_CARR}"
            log "       recomputed = ${RECOMP_PCT}%   n_carriers = ${RECOMP_CARR}"
            log "       delta      = ${DELTA_PP} pp"
        else
            STATUS="FAIL"
            N_FREQ_FAIL=$(( N_FREQ_FAIL + 1 ))
            err "  ${TARGET}%  $VID  FREQUENCY MISMATCH"
            log "       reported   = ${ACT}%   n_carriers = ${N_CARR}"
            log "       recomputed = ${RECOMP_PCT}%   n_carriers = ${RECOMP_CARR}"
            log "       delta      = ${DELTA_PP} pp"
            echo "FREQ_FAIL: $VID reported=$ACT% recomputed=$RECOMP_PCT% delta=${DELTA_PP}pp" \
                >> "$ERRORS"
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TARGET" "$VID" "$ACT" "$RECOMP_PCT" "$DELTA_PP" \
            "$N_CARR" "$RECOMP_CARR" "$STATUS" \
            >> "$FREQ_RECHECK"

    done

    log "  Frequency checks  passed : $N_FREQ_PASS"
    log "  Frequency checks  failed : $N_FREQ_FAIL"

    if [[ "$N_FREQ_FAIL" -eq 0 ]]; then
        ok "All allele frequencies verified correctly against VCF"
    else
        err "$N_FREQ_FAIL allele frequency mismatches — see validation/frequency_recheck.tsv"
    fi

    log "  Full results → validation/frequency_recheck.tsv"
fi

# =============================================================================
# CHECK 6 — CARRIER EOD/LOD RATIO IN EVERY FILE
# =============================================================================
section "CHECK 6: EOD/LOD ratio among SNP carriers (all files)"
log "  For every .pheno file:"
log "    1. Identifies SNP carriers using tmp/carriers/<VID>.txt"
log "    2. Counts how many carriers have EOD (0) vs LOD (1)"
log "    3. Verifies ratio matches the intended EOD:LOD ratio"

# Expected ratio fractions
declare -A RATIO_EOD_P RATIO_LOD_P
RATIO_EOD_P["EOD1_LOD3"]=1; RATIO_LOD_P["EOD1_LOD3"]=3
RATIO_EOD_P["EOD1_LOD2"]=1; RATIO_LOD_P["EOD1_LOD2"]=2
RATIO_EOD_P["EOD2_LOD1"]=2; RATIO_LOD_P["EOD2_LOD1"]=1
RATIO_EOD_P["EOD3_LOD1"]=3; RATIO_LOD_P["EOD3_LOD1"]=1

CARRIER_REPORT="$VAL_DIR/carrier_ratio_check.tsv"
printf "file\tsnp_id\tratio\tn_carriers\tn_eod_carriers\tn_lod_carriers\t\
obs_eod_frac\texp_eod_frac\tdelta\texp_eod_carriers\tstatus\n" \
    > "$CARRIER_REPORT"

N_CARR_PASS=0
N_CARR_FAIL=0
N_CARR_CHECKED=0

while IFS= read -r FBIN_DIR; do
    FBNAME=$(basename "$FBIN_DIR")

    while IFS= read -r SNP_DIR; do
        VID=$(basename "$SNP_DIR")

        # Carrier index file for this SNP
        CARR_FILE="$SIMDIR/tmp/carriers/${VID}.txt"
        if [[ ! -f "$CARR_FILE" ]]; then
            warn "No carrier file for $VID — skipping carrier ratio checks"
            continue
        fi
        N_CARR=$(wc -l < "$CARR_FILE" | tr -d ' ')

        while IFS= read -r RATIO_DIR; do
            RNAME=$(basename "$RATIO_DIR" | sed 's/^ratio_//')

            EP="${RATIO_EOD_P[$RNAME]:-}"
            LP="${RATIO_LOD_P[$RNAME]:-}"
            if [[ -z "$EP" ]]; then
                warn "Unrecognised ratio directory: $RNAME"
                continue
            fi
            DENOM=$(( EP + LP ))

            # Expected EOD carrier count (same rounding rule as simulation)
            EXP_EOD_C=$($AWK "BEGIN{printf \"%d\", int($N_CARR * $EP / $DENOM + 0.5)}")
            EXP_EOD_FRAC=$($AWK "BEGIN{printf \"%.6f\", $EP / $DENOM}")

            # Check every simulation file for this combination
            while IFS= read -r PFILE; do

                N_CARR_CHECKED=$(( N_CARR_CHECKED + 1 ))
                RELPATH="${PFILE#$SIMDIR/}"

                RESULT=$($AWK -F'\t' \
                    -v carr_file="$CARR_FILE" \
                '
                BEGIN {
                    while ((getline idx < carr_file) > 0)
                        carrier_idx[idx+0] = 1
                    close(carr_file)
                }
                NR == 1 { next }
                {
                    idx = NR - 2
                    if (idx in carrier_idx) {
                        if ($2 == 0) eod_c++
                        else         lod_c++
                    }
                }
                END { printf "%d\t%d", eod_c+0, lod_c+0 }
                ' "$PFILE")

                OBS_EOD_C=$(echo "$RESULT" | cut -f1)
                OBS_LOD_C=$(echo "$RESULT" | cut -f2)
                OBS_TOTAL=$(( OBS_EOD_C + OBS_LOD_C ))

                OBS_FRAC=$($AWK "BEGIN{printf \"%.6f\", $OBS_EOD_C / $N_CARR}")
                DELTA_FRAC=$($AWK "BEGIN{
                    d = $OBS_FRAC - $EXP_EOD_FRAC
                    if (d < 0) d = -d
                    printf \"%.6f\", d
                }")

                DIFF_C=$(( OBS_EOD_C - EXP_EOD_C ))
                [[ "$DIFF_C" -lt 0 ]] && DIFF_C=$(( -DIFF_C ))

                # Pass: carrier total correct AND EOD count within ±1 (rounding tolerance)
                if [[ "$OBS_TOTAL" -eq "$N_CARR" && "$DIFF_C" -le 1 ]]; then
                    STATUS="PASS"
                    N_CARR_PASS=$(( N_CARR_PASS + 1 ))
                elif [[ "$OBS_TOTAL" -ne "$N_CARR" ]]; then
                    STATUS="FAIL_CARRIER_TOTAL"
                    N_CARR_FAIL=$(( N_CARR_FAIL + 1 ))
                    echo "CARRIER_TOTAL_FAIL: $RELPATH got=$OBS_TOTAL expected=$N_CARR" \
                        >> "$ERRORS"
                else
                    STATUS="FAIL_RATIO"
                    N_CARR_FAIL=$(( N_CARR_FAIL + 1 ))
                    echo "CARRIER_RATIO_FAIL: $RELPATH $RNAME eod_c=$OBS_EOD_C expected=$EXP_EOD_C" \
                        >> "$ERRORS"
                fi

                printf "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s\n" \
                    "$RELPATH" "$VID" "$RNAME" \
                    "$N_CARR" "$OBS_EOD_C" "$OBS_LOD_C" \
                    "$OBS_FRAC" "$EXP_EOD_FRAC" "$DELTA_FRAC" \
                    "$EXP_EOD_C" "$STATUS" \
                    >> "$CARRIER_REPORT"

            done < <(find "$RATIO_DIR" -name "*.pheno" 2>/dev/null | sort)

        done < <(find "$SNP_DIR" -maxdepth 1 -type d -name "ratio_*" 2>/dev/null | sort)

    done < <(find "$FBIN_DIR" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | sort)

done < <(find "$SIMDIR/03_phenotype_files" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | sort)

log "  Files checked : $N_CARR_CHECKED"
log "  Passed        : $N_CARR_PASS"
log "  Failed        : $N_CARR_FAIL"

if [[ "$N_CARR_FAIL" -eq 0 ]]; then
    ok "ALL $N_CARR_CHECKED carrier ratio checks passed"
else
    err "$N_CARR_FAIL carrier ratio failures — see validation/carrier_ratio_check.tsv"
    grep "CARRIER_" "$ERRORS" | head -5 | sed 's/^/    /' | tee -a "$LOG"
fi

log "  Full results → validation/carrier_ratio_check.tsv"

# =============================================================================
# CHECK 7 — NON-CARRIER EOD/LOD FILL
# =============================================================================
section "CHECK 7: Non-carrier EOD/LOD fill is correct"
log "  Verifying: (global EOD) - (carrier EOD) = non-carrier EOD"
log "  Testing sim_001.pheno for every (SNP × ratio) combination..."

SEL_FILE="$SIMDIR/02_selected_snps/selected_snps.tsv"

N_NC_PASS=0
N_NC_FAIL=0

tail -n +2 "$SEL_FILE" | while IFS=$'\t' read -r \
    TARGET VID CHROM POS REF ALT N_CARR ACT DELTA; do

    CARR_FILE="$SIMDIR/tmp/carriers/${VID}.txt"
    [[ -f "$CARR_FILE" ]] || { warn "No carrier file for $VID — skip"; continue; }

    # Build the directory label (e.g. 5 → freq_05pct)
    FREQ_LABEL=$(printf "freq_%02dpct" "${TARGET%.*}")

    for RNAME in EOD1_LOD3 EOD1_LOD2 EOD2_LOD1 EOD3_LOD1; do

        EP="${RATIO_EOD_P[$RNAME]}"
        LP="${RATIO_LOD_P[$RNAME]}"
        DENOM=$(( EP + LP ))

        EXP_EOD_C=$( $AWK "BEGIN{printf \"%d\", int($N_CARR * $EP / $DENOM + 0.5)}")
        EXP_LOD_C=$(( N_CARR - EXP_EOD_C ))
        EXP_EOD_NC=$(( EXP_EOD - EXP_EOD_C ))
        EXP_LOD_NC=$(( EXP_LOD - EXP_LOD_C ))

        TEST_FILE="$SIMDIR/03_phenotype_files/${FREQ_LABEL}/${VID}/ratio_${RNAME}/sim_001.pheno"

        if [[ ! -f "$TEST_FILE" ]]; then
            warn "  Missing: $FREQ_LABEL/$VID/ratio_${RNAME}/sim_001.pheno"
            continue
        fi

        # Count EOD and LOD among non-carriers
        RESULT=$($AWK -F'\t' \
            -v carr_file="$CARR_FILE" \
        '
        BEGIN {
            while ((getline idx < carr_file) > 0)
                carrier_idx[idx+0] = 1
            close(carr_file)
        }
        NR == 1 { next }
        {
            idx = NR - 2
            if (!(idx in carrier_idx)) {
                if ($2 == 0) eod_nc++
                else         lod_nc++
            }
        }
        END { printf "%d\t%d", eod_nc+0, lod_nc+0 }
        ' "$TEST_FILE")

        OBS_EOD_NC=$(echo "$RESULT" | cut -f1)
        OBS_LOD_NC=$(echo "$RESULT" | cut -f2)

        if [[ "$OBS_EOD_NC" -eq "$EXP_EOD_NC" && \
              "$OBS_LOD_NC" -eq "$EXP_LOD_NC" ]]; then
            N_NC_PASS=$(( N_NC_PASS + 1 ))
            log "  ✓  ${TARGET}%  $RNAME"
            log "       non-carrier EOD=$OBS_EOD_NC (expected $EXP_EOD_NC)  \
LOD=$OBS_LOD_NC (expected $EXP_LOD_NC)"
        else
            N_NC_FAIL=$(( N_NC_FAIL + 1 ))
            err "  ${TARGET}%  $RNAME  NON-CARRIER FILL WRONG"
            log "       got      EOD=$OBS_EOD_NC  LOD=$OBS_LOD_NC"
            log "       expected EOD=$EXP_EOD_NC  LOD=$EXP_LOD_NC"
            echo "NONCARRIER_FAIL: ${TARGET}% $RNAME got_eod=$OBS_EOD_NC \
exp_eod=$EXP_EOD_NC got_lod=$OBS_LOD_NC exp_lod=$EXP_LOD_NC" >> "$ERRORS"
        fi

    done
done

if [[ "$N_NC_FAIL" -eq 0 ]]; then
    ok "All non-carrier EOD/LOD fills are correct"
else
    err "$N_NC_FAIL non-carrier fill failures"
fi

# =============================================================================
# CHECK 8 — SIMULATIONS ARE DIFFERENT FROM EACH OTHER
# =============================================================================
section "CHECK 8: Simulations differ from each other"
log "  Checking that the 20 simulations are NOT identical..."
log "  Method 1: md5 checksum — all 20 files in a combo should be unique"
log "  Method 2: diff sim_001 vs sim_002 — should differ at many rows"
log "  Testing 3 frequency bins × 2 ratios = up to 6 combinations"

N_VAR_PASS=0
N_VAR_FAIL=0

# Test first 3 selected SNPs at 2 ratios each
TESTED_BINS=0
while IFS=$'\t' read -r TARGET VID CHROM POS REF ALT N_CARR ACT DELTA; do
    [[ "$TESTED_BINS" -ge 3 ]] && break
    FREQ_LABEL=$(printf "freq_%02dpct" "${TARGET%.*}")
    TESTED_BINS=$(( TESTED_BINS + 1 ))

    for RNAME in EOD1_LOD3 EOD3_LOD1; do

        SIM_DIR="$SIMDIR/03_phenotype_files/${FREQ_LABEL}/${VID}/ratio_${RNAME}"
        [[ -d "$SIM_DIR" ]] || continue

        # ── Test 1: all checksums are unique ─────────────────────────────────
        N_UNIQUE=$(
            for F in "$SIM_DIR"/sim_*.pheno; do
                # md5sum (Linux) or md5 (macOS)
                md5sum "$F" 2>/dev/null \
                || md5 -q "$F"  2>/dev/null \
                || cksum "$F" | awk '{print $1}'
            done | awk '{print $1}' | sort -u | wc -l | tr -d ' '
        )
        N_SIMS_HERE=$(find "$SIM_DIR" -name "*.pheno" | wc -l | tr -d ' ')

        if [[ "$N_UNIQUE" -eq "$N_SIMS_HERE" ]]; then
            N_VAR_PASS=$(( N_VAR_PASS + 1 ))
            log "  ✓  ${TARGET}%  $RNAME: all $N_SIMS_HERE sims have unique checksums"
        else
            N_VAR_FAIL=$(( N_VAR_FAIL + 1 ))
            err "  ${TARGET}%  $RNAME: only $N_UNIQUE unique checksums out of $N_SIMS_HERE sims"
        fi

        # ── Test 2: sim_001 differs from sim_002 at many lines ────────────────
        F1="$SIM_DIR/sim_001.pheno"
        F2="$SIM_DIR/sim_002.pheno"

        if [[ -f "$F1" && -f "$F2" ]]; then
            N_DIFF=$(diff \
                <(tail -n +2 "$F1") \
                <(tail -n +2 "$F2") \
                | grep "^[<>]" | wc -l | tr -d ' ')

            if [[ "$N_DIFF" -gt 0 ]]; then
                N_VAR_PASS=$(( N_VAR_PASS + 1 ))
                log "  ✓  ${TARGET}%  $RNAME: sim_001 vs sim_002 differ at $N_DIFF assignments"
            else
                N_VAR_FAIL=$(( N_VAR_FAIL + 1 ))
                err "  ${TARGET}%  $RNAME: sim_001 and sim_002 are IDENTICAL — randomisation failed"
                echo "VARIABILITY_FAIL: ${TARGET}% $RNAME sim_001 == sim_002" >> "$ERRORS"
            fi
        fi

    done
done < <(tail -n +2 "$SEL_FILE")

if [[ "$N_VAR_FAIL" -eq 0 ]]; then
    ok "All tested simulations are distinct"
else
    err "$N_VAR_FAIL variability check failures"
fi

# =============================================================================
# SUMMARY STATISTICS TABLE
# =============================================================================
section "SUMMARY STATISTICS"

log ""
log "  Mean observed EOD carrier fraction per ratio:"
printf "  %-15s  %-15s  %-15s  %-10s\n" \
    "Ratio" "Mean_Obs_EOD_frac" "Expected_EOD_frac" "N_files" | tee -a "$LOG"
printf "  %-15s  %-15s  %-15s  %-10s\n" \
    "─────" "─────────────────" "─────────────────" "───────" | tee -a "$LOG"

$AWK -F'\t' '
NR == 1 { next }
{
    ratio  = $3
    obs    = $7 + 0
    exp    = $8 + 0
    sum[ratio]    += obs
    count[ratio]  += 1
    exp_val[ratio] = exp
}
END {
    for (r in sum) {
        mean = sum[r] / count[r]
        printf "  %-15s  %-15.4f  %-15.4f  %-10d\n",
            r, mean, exp_val[r], count[r]
    }
}
' "$CARRIER_REPORT" | sort | tee -a "$LOG"

log ""
log "  Global EOD/LOD count distribution across all files:"
$AWK -F'\t' '
NR == 1 { next }
{
    eod[$3]++
    lod[$4]++
}
END {
    print "  EOD values seen:"
    for (v in eod) printf "    EOD=%s  in %d files\n", v, eod[v]
    print "  LOD values seen:"
    for (v in lod) printf "    LOD=%s  in %d files\n", v, lod[v]
}
' "$VAL_DIR/global_totals_check.tsv" | tee -a "$LOG"

# =============================================================================
# FINAL REPORT
# =============================================================================
log ""
log "======================================================================"
log "  VALIDATION COMPLETE"
log "======================================================================"
log "  Total checks passed : $TOTAL_PASS"
log "  Total checks failed : $TOTAL_FAIL"
log ""
log "  Reports written to : $VAL_DIR/"
log "    validation.log              ← full log of all checks"
log "    errors.txt                  ← failures only"
log "    global_totals_check.tsv     ← EOD/LOD counts for every file"
log "    carrier_ratio_check.tsv     ← carrier ratio for every file"
log "    frequency_recheck.tsv       ← allele frequency verification"
log ""

if [[ "$TOTAL_FAIL" -eq 0 ]]; then
    log "  ✅  ALL CHECKS PASSED — simulation output is correct"
else
    log "  ❌  $TOTAL_FAIL CHECKS FAILED"
    log ""
    log "  Error summary:"
    cat "$ERRORS" | sed 's/^/    /' | tee -a "$LOG"
fi

log "======================================================================"