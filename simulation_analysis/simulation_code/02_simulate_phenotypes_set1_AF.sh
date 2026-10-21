#!/bin/bash
# =============================================================================
#  GBS-GWAS Simulation Analysis — Pure Bash + AWK
#  Fixed: Step 3 selection done entirely in AWK (no sort|head SIGPIPE)
#  5%  bin: selects SNPs strictly between 6%  and 10% (above FIRST_BIN_LO)
#  95% bin: selects SNPs strictly between 90% and 94% (below LAST_BIN_HI)
# =============================================================================

#SBATCH --job-name=gwas_sim
#SBATCH --partition=fasse
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/gwas_sim_%j.out
#SBATCH --error=logs/gwas_sim_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rbalasubramanian@g.harvard.edu

set -euo pipefail

# =============================================================================
# ████  USER INPUTS  ████
# =============================================================================
VCF_FILE="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/\
Nextflow_pipeline_EOD_LOD_VLOD_CDC/11.SNIPPY_MULTI/merged.vcf"

OUTDIR="gwas_simulations_$(date +%Y%m%d)"
TOT_EOD=623
TOT_LOD=917
N_SIMS=20
N_PER_BIN=10
SEED=2024
FREQ_MIN=5
FREQ_MAX=95         # ← CHANGED from 50 to 95
FREQ_STEP=5
TOL="2.5"
FMT="tab"
FIRST_BIN_LO=6      # 5%  bin: select SNPs strictly above this  (avoids pyseer --min-af 0.05)
LAST_BIN_HI=94      # 95% bin: select SNPs strictly below this  (avoids pyseer --max-af 0.95)

# =============================================================================
# ████  STOP EDITING  ████
# =============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)              VCF_FILE="$2";      shift 2 ;;
        --outdir)           OUTDIR="$2";         shift 2 ;;
        --eod)              TOT_EOD="$2";        shift 2 ;;
        --lod)              TOT_LOD="$2";        shift 2 ;;
        --nsims)            N_SIMS="$2";         shift 2 ;;
        --snps-per-bin)     N_PER_BIN="$2";      shift 2 ;;
        --seed)             SEED="$2";           shift 2 ;;
        --freq-min)         FREQ_MIN="$2";       shift 2 ;;
        --freq-max)         FREQ_MAX="$2";       shift 2 ;;
        --freq-step)        FREQ_STEP="$2";      shift 2 ;;
        --tol)              TOL="$2";            shift 2 ;;
        --fmt)              FMT="$2";            shift 2 ;;
        --first-bin-lo)     FIRST_BIN_LO="$2";  shift 2 ;;
        --last-bin-hi)      LAST_BIN_HI="$2";   shift 2 ;;  # ← NEW
        -h|--help)          head -20 "$0";       exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

TOTAL=$((TOT_EOD + TOT_LOD))
N_BINS=$(( (FREQ_MAX - FREQ_MIN) / FREQ_STEP + 1 ))
EXPECTED_FILES=$(( N_BINS * N_PER_BIN * 4 * N_SIMS ))
TOL2=$(awk "BEGIN{printf \"%.4f\", $TOL * 2}")

AWK=awk
command -v gawk &>/dev/null && AWK=gawk

mkdir -p "$OUTDIR/logs" \
         "$OUTDIR/tmp/carriers" \
         "$OUTDIR/01_allele_frequencies" \
         "$OUTDIR/02_selected_snps" \
         "$OUTDIR/03_phenotype_files"

LOG="$OUTDIR/logs/run.log"

log()  { printf "[%s]  %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR:   $*"; exit 1; }

log "======================================================================"
log "  GBS-GWAS Simulation Analysis  (Pure Bash + AWK)"
log "======================================================================"
log "  VCF file         : $VCF_FILE"
log "  Output dir       : $OUTDIR"
log "  EOD (0) / LOD(1) : $TOT_EOD / $TOT_LOD  (total n=$TOTAL)"
log "  Frequency bins   : ${FREQ_MIN}% → ${FREQ_MAX}%  step ${FREQ_STEP}%  ($N_BINS bins)"
log "  SNPs per bin     : $N_PER_BIN"
log "  Ratios (EOD:LOD) : 1:3  1:2  2:1  3:1"
log "  Simulations      : $N_SIMS per combination"
log "  Seed             : $SEED"
log "  Tolerance        : ±${TOL}pp  (fallback ±${TOL2}pp)"
log "  AWK binary       : $AWK"
log "  Missing GT       : treated as REF=0  (SNIPPY convention)"
log "  5%  bin special  : (${FIRST_BIN_LO}%%, 10%%)   strictly above ${FIRST_BIN_LO}%% — pyseer --min-af safe"
log "  95% bin special  : (90%%, ${LAST_BIN_HI}%%)  strictly below ${LAST_BIN_HI}%% — pyseer --max-af safe"
log "  Expected files   : $N_BINS × $N_PER_BIN × 4 × $N_SIMS = $EXPECTED_FILES"
log "======================================================================"

[[ -f "$VCF_FILE" ]] || die "VCF not found: $VCF_FILE"

N_VCF_SAMPLES=$(grep -m1 "^#CHROM" "$VCF_FILE" \
    | tr '\t' '\n' | tail -n +10 | wc -l | tr -d ' ')
log "VCF sample count  : $N_VCF_SAMPLES  (expected $TOTAL)"
[[ "$N_VCF_SAMPLES" -eq "$TOTAL" ]] \
    || warn "Sample count mismatch — got $N_VCF_SAMPLES expected $TOTAL"

# =============================================================================
# STEP 1 — EXTRACT SAMPLE NAMES  (unchanged)
# =============================================================================
log ""
log "── Step 1: Extracting sample names ──────────────────────────────────"

SAMPLES_FILE="$OUTDIR/tmp/samples.txt"

grep -m1 "^#CHROM" "$VCF_FILE" \
    | tr '\t' '\n' \
    | tail -n +10 \
    > "$SAMPLES_FILE"

N_SAMP=$(wc -l < "$SAMPLES_FILE" | tr -d ' ')
log "  $N_SAMP sample names written → tmp/samples.txt"
log "  First 3 : $(head -3 "$SAMPLES_FILE" | tr '\n' '  ')"
log "  Last  3 : $(tail -3 "$SAMPLES_FILE" | tr '\n' '  ')"

# =============================================================================
# STEP 2 — COMPUTE ALLELE FREQUENCIES  (unchanged)
# =============================================================================
log ""
log "── Step 2: Computing allele frequencies ─────────────────────────────"
log "  Streaming full VCF — may take several minutes..."

FREQ_FILE="$OUTDIR/01_allele_frequencies/all_variant_frequencies.tsv"
printf "variant_id\tchrom\tpos\tref\talt\talt_count\tfreq_pct\n" > "$FREQ_FILE"

$AWK -F'\t' '
/^##/ { next }
/^#CHROM/ { n_samples = NF - 9; next }
{
    ref = $4; alt = $5
    if (index(alt, ",") > 0)       next
    if (alt == "." || alt == "*")  next
    if (substr(alt, 1, 1) == "<")  next
    if (length(ref) != 1)          next
    if (length(alt) != 1)          next

    gt_pos = 1
    n_fmt  = split($9, fmt_arr, ":")
    for (f = 1; f <= n_fmt; f++) {
        if (fmt_arr[f] == "GT") { gt_pos = f; break }
    }

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

    if (alt_count == 0 || alt_count == n_samples) next

    freq_pct = alt_count / n_samples * 100
    vid      = $1 "_" $2 "_" ref "_" alt
    printf "%s\t%s\t%s\t%s\t%s\t%d\t%.6f\n",
        vid, $1, $2, ref, alt, alt_count, freq_pct
}
' "$VCF_FILE" >> "$FREQ_FILE"

N_VARS=$(( $(wc -l < "$FREQ_FILE") - 1 ))
log "  $N_VARS biallelic SNPs → 01_allele_frequencies/all_variant_frequencies.tsv"
[[ "$N_VARS" -gt 0 ]] || die "No variants found — check VCF format"

log "  Frequency distribution:"
$AWK -F'\t' '
NR == 1 { next }
{
    f = $7 + 0
    if      (f <  1)  bins["1_<1%"]++
    else if (f <  5)  bins["2_1-5%"]++
    else if (f < 10)  bins["3_5-10%"]++
    else if (f < 20)  bins["4_10-20%"]++
    else if (f < 30)  bins["5_20-30%"]++
    else if (f < 40)  bins["6_30-40%"]++
    else if (f < 50)  bins["7_40-50%"]++
    else if (f < 75)  bins["8_50-75%"]++
    else              bins["9_>75%"]++
}
END {
    for (b in bins) {
        label = substr(b, 3)
        printf "    %-10s  %d variants\n", label, bins[b]
    }
}
' "$FREQ_FILE" | sort | tee -a "$LOG"

# =============================================================================
# STEP 3 — SELECT N_PER_BIN SNPS PER FREQUENCY BIN
#
# UNCHANGED from original EXCEPT for the addition of the 95% bin special rule,
# which mirrors the existing 5% bin special rule.
#
# 5%  bin (is_first): range (FIRST_BIN_LO, 10%)  = (6%, 10%)  — above 6%
#                     delta = f - FIRST_BIN_LO    → SNPs closest to  6% rank first
#                     avoids pyseer --min-af 0.05
#
# 95% bin (is_last):  range (90%, LAST_BIN_HI)   = (90%, 94%) — below 94%
#                     delta = LAST_BIN_HI - f      → SNPs closest to 94% rank first
#                     avoids pyseer --max-af 0.95
#
# All other bins: symmetric ±TOL window, fallback ±2×TOL  (unchanged)
# =============================================================================
log ""
log "── Step 3: Selecting $N_PER_BIN SNPs per frequency bin ─────────────"
log "  Method: pure AWK — no sort|head pipe (avoids SIGPIPE under pipefail)"
log "  5%  bin : strictly (${FIRST_BIN_LO}%%, 10%%)   — SNPs > ${FIRST_BIN_LO}%%"
log "  95% bin : strictly (90%%, ${LAST_BIN_HI}%%) — SNPs < ${LAST_BIN_HI}%%"
log "  Other   : symmetric ±${TOL}pp  (fallback ±${TOL2}pp)"

SEL_FILE="$OUTDIR/02_selected_snps/selected_snps.tsv"
printf "target_pct\tsnp_rank\tvariant_id\tchrom\tpos\tref\talt\tn_carriers\tactual_pct\tdelta_pct\n" \
    > "$SEL_FILE"

$AWK -F'\t' \
    -v freq_min="$FREQ_MIN" \
    -v freq_max="$FREQ_MAX" \
    -v freq_step="$FREQ_STEP" \
    -v tol_pri="$TOL" \
    -v tol_ext="$TOL2" \
    -v n_per_bin="$N_PER_BIN" \
    -v first_bin_lo="$FIRST_BIN_LO" \
    -v last_bin_hi="$LAST_BIN_HI" \
'
BEGIN {
    n_targets = 0
    for (t = freq_min; t <= freq_max; t += freq_step) {
        targets[n_targets] = t
        n_targets++
    }
}

NR == 1 { next }

{
    row++
    vid_arr[row]   = $1
    chrom_arr[row] = $2
    pos_arr[row]   = $3
    ref_arr[row]   = $4
    alt_arr[row]   = $5
    ac_arr[row]    = $6 + 0
    freq_arr[row]  = $7 + 0
}

END {
    for (ti = 0; ti < n_targets; ti++) {
        tgt      = targets[ti]
        is_first = (tgt == freq_min)   # 5%  bin special rule
        is_last  = (tgt == freq_max)   # 95% bin special rule  ← NEW

        n_cands = 0

        if (is_first) {

            # ── 5% bin: open interval (first_bin_lo, freq_min + freq_step) ───
            # i.e. (6%, 10%) strictly above 6% — passes pyseer --min-af 0.05
            # delta = f - first_bin_lo → SNPs closest to 6% rank first
            hi = freq_min + freq_step   # = 10

            for (r = 1; r <= row; r++) {
                f = freq_arr[r]
                if (f > first_bin_lo && f < hi) {
                    d = f - first_bin_lo
                    cand_row[n_cands]   = r
                    cand_delta[n_cands] = d
                    n_cands++
                }
            }

            # Fallback: extend upper bound to 15%
            if (n_cands < n_per_bin) {
                hi2 = freq_min + freq_step * 2   # = 15
                n_cands = 0
                for (r = 1; r <= row; r++) {
                    f = freq_arr[r]
                    if (f > first_bin_lo && f < hi2) {
                        d = f - first_bin_lo
                        cand_row[n_cands]   = r
                        cand_delta[n_cands] = d
                        n_cands++
                    }
                }
                printf "# WARN 5%% bin extended to (%s%%, %s%%) (%d candidates)\n",
                    first_bin_lo, hi2, n_cands > "/dev/stderr"
            }

        } else if (is_last) {

            # ── 95% bin: open interval (freq_max - freq_step, last_bin_hi) ──
            # i.e. (90%, 94%) strictly below 94% — passes pyseer --max-af 0.95
            # delta = last_bin_hi - f → SNPs closest to 94% rank first
            lo = freq_max - freq_step   # = 90

            for (r = 1; r <= row; r++) {
                f = freq_arr[r]
                if (f > lo && f < last_bin_hi) {
                    d = last_bin_hi - f
                    cand_row[n_cands]   = r
                    cand_delta[n_cands] = d
                    n_cands++
                }
            }

            # Fallback: extend lower bound to 85%
            if (n_cands < n_per_bin) {
                lo2 = freq_max - freq_step * 2   # = 85
                n_cands = 0
                for (r = 1; r <= row; r++) {
                    f = freq_arr[r]
                    if (f > lo2 && f < last_bin_hi) {
                        d = last_bin_hi - f
                        cand_row[n_cands]   = r
                        cand_delta[n_cands] = d
                        n_cands++
                    }
                }
                printf "# WARN 95%% bin extended to (%s%%, %s%%) (%d candidates)\n",
                    lo2, last_bin_hi, n_cands > "/dev/stderr"
            }

        } else {

            # ── All other bins: symmetric ±tol_pri  (unchanged) ──────────────
            for (r = 1; r <= row; r++) {
                d = freq_arr[r] - tgt
                if (d < 0) d = -d
                if (d <= tol_pri) {
                    cand_row[n_cands]   = r
                    cand_delta[n_cands] = d
                    n_cands++
                }
            }

            # Fallback: widen to ±tol_ext
            if (n_cands < n_per_bin) {
                n_cands = 0
                for (r = 1; r <= row; r++) {
                    d = freq_arr[r] - tgt
                    if (d < 0) d = -d
                    if (d <= tol_ext) {
                        cand_row[n_cands]   = r
                        cand_delta[n_cands] = d
                        n_cands++
                    }
                }
                printf "# WARN %s%% bin extended to +-%.1f%% (%d candidates)\n",
                    tgt, tol_ext, n_cands > "/dev/stderr"
            }
        }

        if (n_cands == 0) {
            printf "# WARN target=%s%% NO candidates found\n",
                tgt > "/dev/stderr"
            continue
        }

        # ── Partial insertion sort — only sort as far as needed  (unchanged) ──
        n_take = (n_cands < n_per_bin) ? n_cands : n_per_bin

        for (i = 0; i < n_take; i++) {
            best = i
            for (j = i + 1; j < n_cands; j++) {
                if (cand_delta[j] < cand_delta[best]) {
                    best = j
                } else if (cand_delta[j] == cand_delta[best] &&
                           ac_arr[cand_row[j]] > ac_arr[cand_row[best]]) {
                    best = j
                }
            }
            tmp_r            = cand_row[i]
            tmp_d            = cand_delta[i]
            cand_row[i]      = cand_row[best]
            cand_delta[i]    = cand_delta[best]
            cand_row[best]   = tmp_r
            cand_delta[best] = tmp_d
        }

        # ── Emit top n_take rows  (unchanged) ────────────────────────────────
        for (i = 0; i < n_take; i++) {
            r     = cand_row[i]
            delta = cand_delta[i]
            rank  = i + 1
            printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%.4f\n",
                tgt, rank,
                vid_arr[r], chrom_arr[r], pos_arr[r],
                ref_arr[r], alt_arr[r],
                ac_arr[r], freq_arr[r], delta
        }
    }
}
' "$FREQ_FILE" >> "$SEL_FILE"

N_SEL_ROWS=$(( $(wc -l < "$SEL_FILE") - 1 ))
log "  Total SNPs selected : $N_SEL_ROWS  (expected up to $(( N_BINS * N_PER_BIN )))"
[[ "$N_SEL_ROWS" -gt 0 ]] || die "No SNPs selected — check frequency file and AWK step"

# Per-bin summary
log ""
log "  Per-bin selection summary:"
$AWK -F'\t' '
NR == 1 { next }
{
    bin  = $1
    freq = $9 + 0
    count[bin]++
    if (count[bin] == 1 || freq < min_f[bin]) min_f[bin] = freq
    if (count[bin] == 1 || freq > max_f[bin]) max_f[bin] = freq
}
END {
    printf "  %-14s %-12s %-14s %-14s\n",
        "target_pct","n_selected","min_actual_%","max_actual_%"
    printf "  %-14s %-12s %-14s %-14s\n",
        "──────────","──────────","────────────","────────────"
    for (b in count)
        printf "  %-14s %-12d %-14.4f %-14.4f\n",
            b, count[b], min_f[b], max_f[b]
}
' "$SEL_FILE" | sort -k1,1n | tee -a "$LOG"

# Verify 5% bin: all SNPs must be strictly above FIRST_BIN_LO
log ""
log "  Verifying 5%% bin SNPs are strictly above ${FIRST_BIN_LO}%%:"
$AWK -F'\t' -v lo="$FREQ_MIN" -v lo_actual="$FIRST_BIN_LO" '
NR == 1 { next }
$1 + 0 == lo {
    if ($9 + 0 <= lo_actual + 0)
        printf "  ✗  FAIL: %s actual=%.4f%% is NOT above %s%%\n", $3, $9+0, lo_actual
    else
        printf "  ✓  OK:   %s actual=%.4f%% > %s%%\n", $3, $9+0, lo_actual
}
' "$SEL_FILE" | tee -a "$LOG"

# Verify 95% bin: all SNPs must be strictly below LAST_BIN_HI
log ""
log "  Verifying 95%% bin SNPs are strictly below ${LAST_BIN_HI}%%:"
$AWK -F'\t' -v hi="$FREQ_MAX" -v hi_actual="$LAST_BIN_HI" '
NR == 1 { next }
$1 + 0 == hi {
    if ($9 + 0 >= hi_actual + 0)
        printf "  ✗  FAIL: %s actual=%.4f%% is NOT below %s%%\n", $3, $9+0, hi_actual
    else
        printf "  ✓  OK:   %s actual=%.4f%% < %s%%\n", $3, $9+0, hi_actual
}
' "$SEL_FILE" | tee -a "$LOG"

log "  Full list → 02_selected_snps/selected_snps.tsv"

# =============================================================================
# STEP 4 — EXTRACT CARRIER INDICES  (unchanged)
# =============================================================================
log ""
log "── Step 4: Extracting carrier indices ───────────────────────────────"

VIDS_FILE="$OUTDIR/tmp/selected_vids.txt"
cut -f3 "$SEL_FILE" | tail -n +2 | sort -u > "$VIDS_FILE"
N_UNIQUE_VIDS=$(wc -l < "$VIDS_FILE" | tr -d ' ')

log "  Unique SNPs to extract: $N_UNIQUE_VIDS"
log "  Running single-pass VCF scan..."

$AWK -F'\t' \
    -v vids_file="$VIDS_FILE" \
    -v out_dir="$OUTDIR/tmp/carriers" \
    -v n_targets="$N_UNIQUE_VIDS" \
'
BEGIN {
    while ((getline line < vids_file) > 0)
        target_vids[line] = 1
    close(vids_file)
    found = 0
}
/^##/ { next }
/^#CHROM/ { n_samples = NF - 9; next }
{
    ref = $4; alt = $5
    if (index(alt, ",") > 0)       next
    if (alt == "." || alt == "*")  next
    if (substr(alt, 1, 1) == "<")  next
    if (length(ref) != 1)          next
    if (length(alt) != 1)          next

    vid = $1 "_" $2 "_" ref "_" alt
    if (!(vid in target_vids)) next

    gt_pos = 1
    n_fmt  = split($9, fmt_arr, ":")
    for (f = 1; f <= n_fmt; f++) {
        if (fmt_arr[f] == "GT") { gt_pos = f; break }
    }

    outfile = out_dir "/" vid ".txt"
    for (i = 10; i <= NF; i++) {
        n_f = split($i, sf, ":")
        gt  = (n_f >= gt_pos) ? sf[gt_pos] : sf[1]
        gsub(/\|/, "/", gt)
        n_al = split(gt, al, "/")
        for (a = 1; a <= n_al; a++) {
            if (al[a] == "1") {
                print (i - 10) >> outfile
                break
            }
        }
    }
    close(outfile)

    found++
    if (found >= n_targets) exit
}
' "$VCF_FILE"

N_CF_OK=0
N_CF_MISS=0
while IFS= read -r VID; do
    CF="$OUTDIR/tmp/carriers/${VID}.txt"
    if [[ -f "$CF" ]]; then
        NC=$(wc -l < "$CF" | tr -d ' ')
        log "  ✓  $VID  →  $NC carriers"
        N_CF_OK=$(( N_CF_OK + 1 ))
    else
        warn "  ✗  No carrier file for $VID"
        N_CF_MISS=$(( N_CF_MISS + 1 ))
    fi
done < "$VIDS_FILE"

log "  Carrier files created : $N_CF_OK / $N_UNIQUE_VIDS"
[[ "$N_CF_MISS" -gt 0 ]] \
    && warn "$N_CF_MISS missing carrier files — those SNPs will be skipped"

# =============================================================================
# STEP 5 — GENERATE PHENOTYPE FILES  (unchanged)
# =============================================================================
log ""
log "── Step 5: Generating phenotype files ───────────────────────────────"
log "  $N_SEL_ROWS SNPs × 4 ratios × $N_SIMS sims = up to $(( N_SEL_ROWS * 4 * N_SIMS )) files"

N_WRITTEN=0
N_SKIP=0

declare -A RATIOS=(
    ["EOD1_LOD3"]="1 3"
    ["EOD1_LOD2"]="1 2"
    ["EOD2_LOD1"]="2 1"
    ["EOD3_LOD1"]="3 1"
)
RATIO_ORDER=("EOD1_LOD3" "EOD1_LOD2" "EOD2_LOD1" "EOD3_LOD1")

while IFS=$'\t' read -r TARGET RANK VID CHROM POS REF ALT N_CARR ACT DELTA; do

    FREQ_LABEL=$(printf "freq_%02dpct" "$TARGET")
    CF="$OUTDIR/tmp/carriers/${VID}.txt"

    if [[ ! -f "$CF" ]]; then
        warn "  No carrier file for $VID — skipping"
        N_SKIP=$(( N_SKIP + 4 * N_SIMS ))
        continue
    fi

    N_CARR_ACTUAL=$(wc -l < "$CF" | tr -d ' ')

    log ""
    log "  [$FREQ_LABEL rank${RANK}]  $VID"
    log "    actual_freq=${ACT}%  n_carriers=${N_CARR_ACTUAL}"

    for RNAME in "${RATIO_ORDER[@]}"; do

        IFS=' ' read -r EOD_P LOD_P <<< "${RATIOS[$RNAME]}"
        DENOM=$(( EOD_P + LOD_P ))

        N_EOD_C=$(  $AWK "BEGIN{printf \"%d\", int($N_CARR_ACTUAL * $EOD_P / $DENOM + 0.5)}")
        N_LOD_C=$(( N_CARR_ACTUAL - N_EOD_C ))
        N_EOD_NC=$(( TOT_EOD - N_EOD_C ))
        N_LOD_NC=$(( TOT_LOD - N_LOD_C ))

        if [[ "$N_EOD_NC" -lt 0 || "$N_LOD_NC" -lt 0 ]]; then
            warn "    SKIP $RNAME — infeasible"
            N_SKIP=$(( N_SKIP + N_SIMS ))
            continue
        fi

        log "    $RNAME  carriers → EOD=$N_EOD_C LOD=$N_LOD_C | non-carriers → EOD=$N_EOD_NC LOD=$N_LOD_NC"

        RATIO_DIR="$OUTDIR/03_phenotype_files/$FREQ_LABEL/$VID/ratio_${RNAME}"
        mkdir -p "$RATIO_DIR"

        COMBO_HASH=$(printf "%s|%s" "$VID" "$RNAME" | cksum | cut -d' ' -f1)
        COMBO_SEED=$(( (SEED + COMBO_HASH) % 16777216 ))

        $AWK \
            -v carriers_file="$CF" \
            -v samples_file="$SAMPLES_FILE" \
            -v n_eod_c="$N_EOD_C" \
            -v n_eod_nc="$N_EOD_NC" \
            -v tot_eod="$TOT_EOD" \
            -v tot_lod="$TOT_LOD" \
            -v n_sims="$N_SIMS" \
            -v combo_seed="$COMBO_SEED" \
            -v ratio_dir="$RATIO_DIR" \
            -v fmt="$FMT" \
        '
        BEGIN {
            ns = 0
            while ((getline sname < samples_file) > 0) {
                samples[ns] = sname; ns++
            }
            close(samples_file)

            nc = 0
            while ((getline idx < carriers_file) > 0) {
                carriers[nc] = idx + 0; is_carrier[idx+0] = 1; nc++
            }
            close(carriers_file)

            nnc = 0
            for (j = 0; j < ns; j++) {
                if (!is_carrier[j]) { noncarriers[nnc] = j; nnc++ }
            }

            for (sim = 1; sim <= n_sims; sim++) {

                carrier_seed    = (combo_seed + sim * 999983)  % 16777216
                noncarrier_seed = (combo_seed + sim * 1000003) % 16777216

                srand(carrier_seed)
                for (j = 0; j < nc; j++) ct[j] = carriers[j]
                for (j = nc - 1; j > 0; j--) {
                    k = int(rand() * (j + 1))
                    tmp = ct[j]; ct[j] = ct[k]; ct[k] = tmp
                }

                srand(noncarrier_seed)
                for (j = 0; j < nnc; j++) nct[j] = noncarriers[j]
                for (j = nnc - 1; j > 0; j--) {
                    k = int(rand() * (j + 1))
                    tmp = nct[j]; nct[j] = nct[k]; nct[k] = tmp
                }

                for (k in pheno) delete pheno[k]
                for (j = 0;        j < n_eod_c;  j++) pheno[ct[j]]  = 0
                for (j = n_eod_c;  j < nc;       j++) pheno[ct[j]]  = 1
                for (j = 0;        j < n_eod_nc; j++) pheno[nct[j]] = 0
                for (j = n_eod_nc; j < nnc;      j++) pheno[nct[j]] = 1

                eod_n = 0; lod_n = 0
                for (j = 0; j < ns; j++) {
                    if (pheno[j] == 0) eod_n++; else lod_n++
                }
                if (eod_n != tot_eod || lod_n != tot_lod) {
                    printf "ASSERTION FAILED sim=%d EOD=%d LOD=%d\n",
                        sim, eod_n, lod_n > "/dev/stderr"
                    exit 1
                }

                outfile = ratio_dir "/" sprintf("sim_%03d.pheno", sim)
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

        N_WRITTEN=$(( N_WRITTEN + N_SIMS ))
        log "      → $N_SIMS files written"

    done
done < <(tail -n +2 "$SEL_FILE")

# =============================================================================
# STEP 6 — SPOT CHECK  (unchanged)
# =============================================================================
log ""
log "── Step 6: Spot-checking 10 random phenotype files ─────────────────"

PASS=0; FAIL=0
while IFS= read -r PFILE; do
    [[ -f "$PFILE" ]] || { log "  ✗  MISSING: $PFILE"; FAIL=$((FAIL+1)); continue; }
    READ_EOD=$(awk -F'\t' 'NR>1 && $2==0{c++} END{print c+0}' "$PFILE")
    READ_LOD=$(awk -F'\t' 'NR>1 && $2==1{c++} END{print c+0}' "$PFILE")
    LABEL="$(basename "$(dirname "$PFILE")")/$(basename "$PFILE")"
    if [[ "$READ_EOD" -eq "$TOT_EOD" && "$READ_LOD" -eq "$TOT_LOD" ]]; then
        log "  ✓  $LABEL  EOD=$READ_EOD LOD=$READ_LOD"
        PASS=$((PASS+1))
    else
        log "  ✗  $LABEL  GOT EOD=$READ_EOD LOD=$READ_LOD  WANT $TOT_EOD/$TOT_LOD"
        FAIL=$((FAIL+1))
    fi
done < <(find "$OUTDIR/03_phenotype_files" -name "*.pheno" 2>/dev/null \
          | sort -R | head -10)

# =============================================================================
# FINAL SUMMARY
# =============================================================================
FINAL_N=$(find "$OUTDIR/03_phenotype_files" -name "*.pheno" 2>/dev/null \
           | wc -l | tr -d ' ')

log ""
log "======================================================================"
log "  FINAL SUMMARY"
log "======================================================================"
log "  SNPs selected    : $N_SEL_ROWS"
log "  Phenotype files  : $FINAL_N  (expected $EXPECTED_FILES)"
log "  Files skipped    : $N_SKIP"
log "  Spot check       : $PASS passed  /  $FAIL failed"
log "  5%%  bin         : SNPs from (${FIRST_BIN_LO}%%, 10%%) — safe for pyseer --min-af"
log "  95%% bin         : SNPs from (90%%, ${LAST_BIN_HI}%%) — safe for pyseer --max-af"

if [[ "$FAIL" -eq 0 && "$FINAL_N" -gt 0 ]]; then
    log ""
    log "  ✅  All checks passed — output is ready"
else
    [[ "$FINAL_N" -eq 0 ]] && log "  ❌  No files written — check $LOG"
    [[ "$FAIL"    -gt 0 ]] && log "  ❌  Spot-check failures — check $LOG"
fi

log ""
log "======================================================================"
log "  DONE  $(date)"
log "======================================================================"
