#!/bin/bash
# =============================================================================
#  run_simulation_homoplasy.sh
#  GBS-GWAS Simulation — Homoplasy-Stratified
#
#  CHANGE vs previous working version:
#    Step 3 now RANDOMLY SAMPLES from AF ranges instead of selecting
#    the SNPs closest to a target frequency.
#
#    AF ranges:
#      5%  bin  →  [6.0,  7.5]%   (lower bound 6 to pass pyseer --min-af)
#      10% bin  →  [7.5, 12.5]%
#      15% bin  →  [12.5,17.5]%
#      20% bin  →  [17.5,22.5]%
#      25% bin  →  [22.5,27.5]%
#      40% bin  →  [37.5,42.5]%
#
#    Format of AF_BINS: "lo:hi:target_pct" separated by |
#      target_pct = the label used for directory naming (5,10,15...)
#      lo/hi      = inclusive AF range bounds in %
#
#  All other steps unchanged from working version.
# =============================================================================

#SBATCH --job-name=gwas_sim_homo
#SBATCH --partition=fasse
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/gwas_sim_homo_%j.out
#SBATCH --error=logs/gwas_sim_homo_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rbalasubramanian@g.harvard.edu

set -euo pipefail

# =============================================================================
# ████  USER INPUTS  ████
# =============================================================================

VCF_FILE="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/\
Nextflow_pipeline_EOD_LOD_VLOD_CDC/11.SNIPPY_MULTI/merged.vcf"

HOMOPLASY_FILE="/n/hanage_gbs_gwas_l3/Lab/Variant_Simulation/homoplasy_results.tsv"

OUTDIR="gwas_simulations_hi_$(date +%Y%m%d)"
TOT_EOD=623
TOT_LOD=917
N_SIMS=20
N_PER_COMBO=3
SEED=2024
FMT="tab"

# AF bins: "lo:hi:target_pct" separated by |
# lo and hi are INCLUSIVE bounds in percentage points
# target_pct is used for directory naming (freq_05pct, freq_10pct, ...)
AF_BINS="6:7.5:5|7.5:12.5:10|12.5:17.5:15|17.5:22.5:20|22.5:27.5:25|37.5:42.5:40"

# HI bins: "lo:hi:label" separated by |  (bounds inclusive)
HI_BINS="0:0.1:hi_00_10|0.45:0.55:hi_45_55|0.60:0.70:hi_60_70|0.70:0.80:hi_70_80|0.80:0.90:hi_80_90|0.90:1.0:hi_90_100"

# =============================================================================
# ████  STOP EDITING  ████
# =============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)           VCF_FILE="$2";       shift 2 ;;
        --homoplasy)     HOMOPLASY_FILE="$2";  shift 2 ;;
        --outdir)        OUTDIR="$2";          shift 2 ;;
        --eod)           TOT_EOD="$2";         shift 2 ;;
        --lod)           TOT_LOD="$2";         shift 2 ;;
        --nsims)         N_SIMS="$2";          shift 2 ;;
        --per-combo)     N_PER_COMBO="$2";     shift 2 ;;
        --seed)          SEED="$2";            shift 2 ;;
        --fmt)           FMT="$2";             shift 2 ;;
        -h|--help)       head -40 "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

TOTAL=$(( TOT_EOD + TOT_LOD ))
N_AF_BINS=$(echo "$AF_BINS" | tr '|' '\n' | wc -l | tr -d ' ')
N_HI_BINS=$(echo "$HI_BINS" | tr '|' '\n' | wc -l | tr -d ' ')
MAX_SNPS=$(( N_AF_BINS * N_HI_BINS * N_PER_COMBO ))
MAX_FILES=$(( MAX_SNPS * 4 * N_SIMS ))

# Pre-filter bounds for Step 2.5
# LO: min(AF_LO) - 0.01 so that inclusive lower bounds are captured
# HI: max(AF_HI) + 1.0  as a buffer
PREFILTER_LO=$(echo "$AF_BINS" | tr '|' '\n' | \
    awk -F: 'BEGIN{m=99999}{v=$1+0; if(NR==1||v<m)m=v}END{printf "%.4f", m-0.01}')
PREFILTER_HI=$(echo "$AF_BINS" | tr '|' '\n' | \
    awk -F: 'BEGIN{m=0}{v=$2+0; if(v>m)m=v}END{printf "%.4f", m+1.0}')

AWK=awk
command -v gawk &>/dev/null && AWK=gawk

mkdir -p "$OUTDIR/logs" \
         "$OUTDIR/tmp/carriers" \
         "$OUTDIR/tmp/cands" \
         "$OUTDIR/01_allele_frequencies" \
         "$OUTDIR/02_selected_snps" \
         "$OUTDIR/03_phenotype_files"

LOG="$OUTDIR/logs/run.log"

log()  { printf "[%s]  %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR:   $*"; exit 1; }

log "======================================================================"
log "  GBS-GWAS Simulation — Homoplasy-Stratified (random AF sampling)"
log "  (run_simulation_homoplasy.sh)"
log "======================================================================"
log "  VCF file         : $VCF_FILE"
log "  Homoplasy file   : $HOMOPLASY_FILE"
log "  Output dir       : $OUTDIR"
log "  EOD / LOD        : $TOT_EOD / $TOT_LOD  (n=$TOTAL)"
log "  AF bins          : $N_AF_BINS  (random sampling from range)"
log "  AF ranges        : $AF_BINS"
log "  HI bins          : $N_HI_BINS"
log "  SNPs per combo   : up to $N_PER_COMBO  (fewer OK)"
log "  Max SNPs         : $N_AF_BINS × $N_HI_BINS × $N_PER_COMBO = $MAX_SNPS"
log "  Ratios (EOD:LOD) : 1:3  1:2  2:1  3:1"
log "  Simulations      : $N_SIMS per combination"
log "  Seed             : $SEED"
log "  Pre-filter AF    : ($PREFILTER_LO, $PREFILTER_HI)%%"
log "  Max files        : $MAX_SNPS × 4 × $N_SIMS = $MAX_FILES"
log "======================================================================"

[[ -f "$VCF_FILE" ]]       || die "VCF not found: $VCF_FILE"
[[ -f "$HOMOPLASY_FILE" ]] || die "Homoplasy file not found: $HOMOPLASY_FILE"

N_VCF_SAMPLES=$(grep -m1 "^#CHROM" "$VCF_FILE" \
    | tr '\t' '\n' | tail -n +10 | wc -l | tr -d ' ')
log "VCF sample count  : $N_VCF_SAMPLES  (expected $TOTAL)"
[[ "$N_VCF_SAMPLES" -eq "$TOTAL" ]] \
    || warn "Sample count mismatch — got $N_VCF_SAMPLES expected $TOTAL"

# =============================================================================
# STEP 0 — BUILD HI LOOKUP FILE  (unchanged)
# =============================================================================
log ""
log "── Step 0: Building HI lookup file ──────────────────────────────────"

HI_LOOKUP="$OUTDIR/tmp/hi_lookup.tsv"

$AWK -F'\t' '
NR == 1 { next }
NF < 10 { next }
{
    for (i = 1; i <= NF; i++) gsub(/\r/, "", $i)
    if ($10 == "" || $10 == "HI") next
    printf "%s_%s_%s_%s\t%s\n", $1, $2, $4, $5, $10
}
' "$HOMOPLASY_FILE" > "$HI_LOOKUP"

N_HI=$(wc -l < "$HI_LOOKUP" | tr -d ' ')
log "  $N_HI entries → tmp/hi_lookup.tsv"
log "  Example: $(head -1 "$HI_LOOKUP")"
[[ "$N_HI" -gt 0 ]] \
    || die "HI lookup is empty — check cols: CHROM(1) POS(2) REF(4) ALT(5) HI(10)"

# =============================================================================
# STEP 1 — EXTRACT SAMPLE NAMES  (unchanged)
# =============================================================================
log ""
log "── Step 1: Extracting sample names ──────────────────────────────────"

SAMPLES_FILE="$OUTDIR/tmp/samples.txt"
grep -m1 "^#CHROM" "$VCF_FILE" | tr '\t' '\n' | tail -n +10 > "$SAMPLES_FILE"

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
# STEP 2.5 — PRE-JOIN: FREQ + HI → SMALL COMBINED FILE  (bounds updated)
#
# Uses NR==FNR two-file join — no getline in BEGIN.
# Pre-filter bounds now come from AF_BINS (PREFILTER_LO / PREFILTER_HI).
# Output columns: vid(1) chrom(2) pos(3) ref(4) alt(5) ac(6) freq(7) hi(8)
# =============================================================================
log ""
log "── Step 2.5: Pre-joining AF data with HI data ───────────────────────"
log "  AF pre-filter : f > $PREFILTER_LO%%  AND  f < $PREFILTER_HI%%"

FREQ_HI_FILE="$OUTDIR/tmp/freq_hi_filtered.tsv"

$AWK -F'\t' \
    -v af_lo="$PREFILTER_LO" \
    -v af_hi="$PREFILTER_HI" \
'
NR == FNR {
    if ($1 != "" && $2 != "") hi_map[$1] = $2 + 0
    next
}
FNR == 1 { next }
{
    f = $7 + 0
    if (f <= af_lo) next
    if (f >= af_hi) next
    if (!($1 in hi_map)) next
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.6f\n",
        $1, $2, $3, $4, $5, $6, $7, hi_map[$1]
}
' "$HI_LOOKUP" "$FREQ_FILE" > "$FREQ_HI_FILE"

N_FREQ_HI=$(wc -l < "$FREQ_HI_FILE" | tr -d ' ')
log "  $N_FREQ_HI variants in combined file → tmp/freq_hi_filtered.tsv"
[[ "$N_FREQ_HI" -gt 0 ]] \
    || die "Combined freq+HI file is empty — check HOMOPLASY_FILE format"
log "  Example: $(head -1 "$FREQ_HI_FILE")"

# =============================================================================
# STEP 3 — RANDOMLY SAMPLE N_PER_COMBO SNPS PER (AF BIN × HI BIN)
#
# KEY CHANGE from previous version:
#   OLD: filter ±tol around target → sort by delta → take closest N
#   NEW: filter within [AF_LO, AF_HI] range → RANDOMLY SAMPLE N
#
# Random sampling method:
#   1. Filter FREQ_HI_FILE to candidates in [AF_LO, AF_HI] AND [HI_LO, HI_HI]
#   2. Add a seeded random key per row using awk (reproducible via SEED)
#   3. Sort by that random key (random shuffle)
#   4. awk NR<=N takes first N rows — reads ALL of sort output (no SIGPIPE)
#
# Seed per combination: SEED + af_idx * 1000 + hi_idx * 13
# This ensures different (AF×HI) combinations get different random selections
# while remaining fully reproducible.
#
# delta_pct in SEL_FILE = |actual_freq - bin_center| (informational only)
#
# FREQ_HI_FILE cols: vid(1) chrom(2) pos(3) ref(4) alt(5) ac(6) freq(7) hi(8)
# SEL_FILE cols:     target_pct hi_bin rank vid chrom pos ref alt ac freq delta hi
# =============================================================================
log ""
log "── Step 3: Randomly sampling $N_PER_COMBO SNPs per (AF bin × HI bin) ─"
log "  Method  : random sampling within AF range (seeded, reproducible)"
log "  Input   : $FREQ_HI_FILE  ($N_FREQ_HI rows)"
log "  Combos  : $N_AF_BINS AF × $N_HI_BINS HI = $(( N_AF_BINS * N_HI_BINS ))"
log ""
log "  AF bin ranges:"
echo "$AF_BINS" | tr '|' '\n' | \
    awk -F: '{printf "    target %s%%  →  [%s%%, %s%%]\n", $3, $1, $2}' | \
    tee -a "$LOG"

SEL_FILE="$OUTDIR/02_selected_snps/selected_snps.tsv"
printf "target_pct\thi_bin\tsnp_rank\tvariant_id\tchrom\tpos\tref\talt\tn_carriers\tactual_pct\tdelta_pct\tactual_hi\n" \
    > "$SEL_FILE"

# Parse AF bins into bash arrays
declare -a AF_LO_ARR AF_HI_ARR AF_TGT_ARR
while IFS=':' read -r AL AH AT; do
    AF_LO_ARR+=("$AL")
    AF_HI_ARR+=("$AH")
    AF_TGT_ARR+=("$AT")
done < <(echo "$AF_BINS" | tr '|' '\n')
N_AF_ARR=${#AF_LO_ARR[@]}

# Parse HI bins into bash arrays (same as before)
declare -a HI_LO_ARR HI_HI_ARR HI_LBL_ARR
while IFS=':' read -r HL HH LBL; do
    HI_LO_ARR+=("$HL")
    HI_HI_ARR+=("$HH")
    HI_LBL_ARR+=("$LBL")
done < <(echo "$HI_BINS" | tr '|' '\n')
N_HI_ARR=${#HI_LO_ARR[@]}

TOTAL_SELECTED=0

for (( af_idx=0; af_idx < N_AF_ARR; af_idx++ )); do

    AF_LO="${AF_LO_ARR[$af_idx]}"
    AF_HI="${AF_HI_ARR[$af_idx]}"
    AF_TGT="${AF_TGT_ARR[$af_idx]}"

    # Bin center — used to compute delta_pct (informational)
    BIN_CENTER=$(awk "BEGIN{printf \"%.6f\", ($AF_LO + $AF_HI) / 2.0}")

    for (( hi_idx=0; hi_idx < N_HI_ARR; hi_idx++ )); do

        HL="${HI_LO_ARR[$hi_idx]}"
        HH="${HI_HI_ARR[$hi_idx]}"
        LBL="${HI_LBL_ARR[$hi_idx]}"

        CAND_TMP="$OUTDIR/tmp/cands/cand_af${AF_TGT}_${LBL}.tsv"

        # ── Filter: all candidates in range [AF_LO, AF_HI] AND [HL, HH] ──────
        # Inclusive bounds on both ends.
        # No delta prepended — we random-sample, not delta-sort.
        # Output format same as FREQ_HI_FILE:
        #   vid(1) chrom(2) pos(3) ref(4) alt(5) ac(6) freq(7) hi(8)
        $AWK -F'\t' \
            -v alo="$AF_LO" -v ahi="$AF_HI" \
            -v hlo="$HL"    -v hhi="$HH" \
        '{
            f = $7+0; h = $8+0
            if (f < alo || f > ahi) next    # inclusive AF bounds
            if (h < hlo || h > hhi) next    # inclusive HI bounds
            print
        }' "$FREQ_HI_FILE" > "$CAND_TMP"

        N_CANDS=$(wc -l < "$CAND_TMP" | tr -d ' ')

        if [[ "$N_CANDS" -eq 0 ]]; then
            warn "  No candidates for AF=${AF_TGT}% (${AF_LO}–${AF_HI}%) HI=${LBL} — skipping"
            continue
        fi

        # ── Random sample ─────────────────────────────────────────────────────
        # Unique seed per (AF bin × HI bin) combination
        SAMPLE_SEED=$(( (SEED + af_idx * 1000 + hi_idx * 13) % 16777216 ))

        # Step 1: awk assigns a seeded random float key to each candidate row
        # Step 2: sort by that key → random shuffle of candidates
        # Step 3: awk NR<=nmax reads ALL of sort's output (no SIGPIPE),
        #         writes first nmax rows to SEL_FILE, returns count

        N_TAKEN=$(
            $AWK -F'\t' -v seed="$SAMPLE_SEED" '
            BEGIN { srand(seed) }
            { printf "%.10f\t%s\n", rand(), $0 }
            ' "$CAND_TMP" | \
            sort -t$'\t' -k1,1n | \
            $AWK -F'\t' \
                -v tgt="$AF_TGT" \
                -v lbl="$LBL" \
                -v nmax="$N_PER_COMBO" \
                -v bcenter="$BIN_CENTER" \
                -v sel="$SEL_FILE" \
            '
            # After prepending random key, columns shift by 1:
            # rnd(1) vid(2) chrom(3) pos(4) ref(5) alt(6) ac(7) freq(8) hi(9)
            NR <= nmax {
                freq  = $8 + 0
                delta = freq - bcenter
                if (delta < 0) delta = -delta
                printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%.6f\t%.4f\t%.6f\n",
                    tgt, lbl, NR,
                    $2, $3, $4, $5, $6, $7, freq, delta, $9 >> sel
                n++
            }
            END { print n+0 }
            '
        )

        log "    AF=${AF_TGT}% [${AF_LO}–${AF_HI}%%]  ${LBL}  →  ${N_CANDS} candidates → sampled ${N_TAKEN}  (seed=${SAMPLE_SEED})"
        TOTAL_SELECTED=$(( TOTAL_SELECTED + N_TAKEN ))

    done   # HI bin loop
done   # AF bin loop

N_SEL_ROWS=$(( $(wc -l < "$SEL_FILE") - 1 ))
log ""
log "  Total SNPs selected : $N_SEL_ROWS  (max possible $MAX_SNPS)"
[[ "$N_SEL_ROWS" -gt 0 ]] || die "No SNPs selected — check HI bins and AF ranges"
log "  Full list → 02_selected_snps/selected_snps.tsv"

# Quick summary of actual AF values selected per bin
log ""
log "  Actual AF range achieved per target bin:"
$AWK -F'\t' '
NR == 1 { next }
{
    bin = $1
    f   = $10 + 0
    cnt[bin]++
    if (cnt[bin] == 1 || f < mn[bin]) mn[bin] = f
    if (cnt[bin] == 1 || f > mx[bin]) mx[bin] = f
    sm[bin] += f
}
END {
    printf "  %-12s  %-10s  %-10s  %-10s  %-10s\n",
        "target_%","n","min_AF","mean_AF","max_AF"
    printf "  %-12s  %-10s  %-10s  %-10s  %-10s\n",
        "────────","─","──────","───────","──────"
    for (b in cnt)
        printf "  %-12s  %-10d  %-10.3f  %-10.3f  %-10.3f\n",
            b, cnt[b], mn[b], sm[b]/cnt[b], mx[b]
}
' "$SEL_FILE" | sort -k1,1n | tee -a "$LOG"

# =============================================================================
# STEP 4 — EXTRACT CARRIER INDICES  (unchanged)
# NOTE: variant_id is column 4 in selected_snps.tsv
# =============================================================================
log ""
log "── Step 4: Extracting carrier indices ───────────────────────────────"

VIDS_FILE="$OUTDIR/tmp/selected_vids.txt"
cut -f4 "$SEL_FILE" | tail -n +2 | sort -u > "$VIDS_FILE"
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

while IFS=$'\t' read -r TARGET HI_BIN RANK VID CHROM POS REF ALT \
                              N_CARR ACT DELTA ACTUAL_HI; do

    FREQ_LABEL=$(printf "freq_%02dpct" "$TARGET")
    CF="$OUTDIR/tmp/carriers/${VID}.txt"

    if [[ ! -f "$CF" ]]; then
        warn "  No carrier file for $VID — skipping"
        N_SKIP=$(( N_SKIP + 4 * N_SIMS ))
        continue
    fi

    N_CARR_ACTUAL=$(wc -l < "$CF" | tr -d ' ')

    log ""
    log "  [$FREQ_LABEL / $HI_BIN  rank${RANK}]  $VID"
    log "    actual_freq=${ACT}%   n_carriers=${N_CARR_ACTUAL}   HI=${ACTUAL_HI}"

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

        RATIO_DIR="$OUTDIR/03_phenotype_files/$FREQ_LABEL/$HI_BIN/$VID/ratio_${RNAME}"
        mkdir -p "$RATIO_DIR"

        COMBO_HASH=$(printf "%s|%s|%s" "$VID" "$HI_BIN" "$RNAME" | cksum | cut -d' ' -f1)
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
log "  AF bins (random)  : $AF_BINS"
log "  HI bins           : $N_HI_BINS"
log "  SNPs per combo    : up to $N_PER_COMBO  (fewer OK)"
log "  SNPs selected     : $N_SEL_ROWS  (max possible $MAX_SNPS)"
log "  Phenotype files   : $FINAL_N  (max possible $MAX_FILES)"
log "  Files skipped     : $N_SKIP"
log "  Spot check        : $PASS passed  /  $FAIL failed"
log ""
log "  Output layout:"
log "    $OUTDIR/03_phenotype_files/"
log "    ├── freq_05pct/hi_00_10/<SNP>/ratio_EOD1_LOD3/sim_001.pheno"
log "    └── freq_40pct/hi_90_100/<SNP>/..."

if [[ "$FAIL" -eq 0 && "$FINAL_N" -gt 0 ]]; then
    log "  ✅  All checks passed — output is ready"
else
    [[ "$FINAL_N" -eq 0 ]] && log "  ❌  No files written — check $LOG"
    [[ "$FAIL"    -gt 0 ]] && log "  ❌  Spot-check failures — check $LOG"
fi

log "======================================================================"
log "  DONE  $(date)"
log "======================================================================"
