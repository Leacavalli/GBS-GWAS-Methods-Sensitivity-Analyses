#!/bin/bash
# =============================================================================
#  GBS-GWAS Simulation Analysis — Pure Bash + AWK
#  No Python, no R, no conda required
#  Tools used: bash awk grep sort cksum find (all standard on O2/Rocky8)
# =============================================================================
#
#  WHAT THIS DOES
#  ──────────────
#  1. Parses merged.vcf → sample names + allele frequencies
#     Missing genotypes (./.) treated as REF=0  (SNIPPY convention)
#  2. Selects 1 SNP per bin: 5% 10% 15% 20% 25% 30% 35% 40% 45% 50%
#  3. For each SNP × 4 ratios × 20 simulations:
#       Ratios (EOD:LOD among SNP carriers): 1:3  1:2  2:1  3:1
#       Global totals FIXED: EOD=623  LOD=917  in every file
#  → 10 × 4 × 20 = 800 phenotype files
#
#  USAGE
#  ─────
#  bash run_simulation.sh                        ← use defaults below
#  bash run_simulation.sh --vcf /path/to/v.vcf   ← override VCF path
#  sbatch run_simulation.sh                      ← cluster submission
#
# =============================================================================

## ── SLURM  (remove one '#' per line to activate) ─────────────────────────────
##SBATCH --job-name=gwas_sim
##SBATCH --partition=short
##SBATCH --time=12:00:00
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=2
##SBATCH --mem=32G
##SBATCH --output=logs/gwas_sim_%j.out
##SBATCH --error=logs/gwas_sim_%j.err
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=YOUR_EMAIL@partners.org

set -euo pipefail

# =============================================================================
# ████  USER INPUTS — edit these  ████
# =============================================================================

VCF_FILE="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/\
Nextflow_pipeline_EOD_LOD_VLOD_CDC/11.SNIPPY_MULTI/merged.vcf"

OUTDIR="gwas_simulations_$(date +%Y%m%d)"

TOT_EOD=623      # total EOD isolates  → phenotype value 0
TOT_LOD=917      # total LOD isolates  → phenotype value 1
N_SIMS=20        # random draws per (SNP × ratio) combination
SEED=2024        # base random seed
FREQ_MIN=5       # lowest  target allele frequency (%)
FREQ_MAX=50      # highest target allele frequency (%) → 10 SNPs
FREQ_STEP=5      # step between bins (%)
TOL="2.5"        # ±pp window for SNP frequency matching
FMT="tab"        # tab  → sample / phenotype
                 # plink → FID / IID / PHENOTYPE

# =============================================================================
# ████  STOP EDITING  ████
# =============================================================================

# ── Optional CLI overrides ────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)       VCF_FILE="$2";  shift 2 ;;
        --outdir)    OUTDIR="$2";    shift 2 ;;
        --eod)       TOT_EOD="$2";   shift 2 ;;
        --lod)       TOT_LOD="$2";   shift 2 ;;
        --nsims)     N_SIMS="$2";    shift 2 ;;
        --seed)      SEED="$2";      shift 2 ;;
        --freq-min)  FREQ_MIN="$2";  shift 2 ;;
        --freq-max)  FREQ_MAX="$2";  shift 2 ;;
        --freq-step) FREQ_STEP="$2"; shift 2 ;;
        --tol)       TOL="$2";       shift 2 ;;
        --fmt)       FMT="$2";       shift 2 ;;
        -h|--help)   head -35 "$0";  exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── Derived values ─────────────────────────────────────────────────────────
TOTAL=$((TOT_EOD + TOT_LOD))
N_BINS=$(( (FREQ_MAX - FREQ_MIN) / FREQ_STEP + 1 ))
EXPECTED_FILES=$(( N_BINS * 4 * N_SIMS ))    # 10 × 4 × 20 = 800

# ── AWK: prefer gawk for better srand() quality ────────────────────────────
AWK=awk
command -v gawk &>/dev/null && AWK=gawk

# ── Directories ────────────────────────────────────────────────────────────
mkdir -p "$OUTDIR"/{logs,tmp/carriers,\
01_allele_frequencies,02_selected_snps,03_phenotype_files}
LOG="$OUTDIR/logs/run.log"

# ── Logging helpers ────────────────────────────────────────────────────────
log()  { printf "[%s]  %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR:   $*"; exit 1; }

log "======================================================================"
log "  GBS-GWAS Simulation Analysis  (Pure Bash + AWK)"
log "======================================================================"
log "  VCF              : $VCF_FILE"
log "  Output dir       : $OUTDIR"
log "  EOD (0) / LOD(1) : $TOT_EOD / $TOT_LOD  (n=$TOTAL)"
log "  Frequency bins   : ${FREQ_MIN}% → ${FREQ_MAX}%  step ${FREQ_STEP}%  ($N_BINS bins)"
log "  Ratios (EOD:LOD) : 1:3  1:2  2:1  3:1  (among SNP carriers)"
log "  Simulations      : $N_SIMS per combination"
log "  Seed             : $SEED"
log "  Tolerance        : ±${TOL} pp"
log "  AWK              : $AWK"
log "  Missing GT       : treated as REF=0  (SNIPPY convention)"
log "  Expected files   : $N_BINS × 4 × $N_SIMS = $EXPECTED_FILES"
log "======================================================================"

# ── Input validation ───────────────────────────────────────────────────────
[[ -f "$VCF_FILE" ]] || die "VCF not found: $VCF_FILE"

N_VCF_SAMPLES=$(grep -m1 "^#CHROM" "$VCF_FILE" \
                | tr '\t' '\n' | tail -n +10 | wc -l | tr -d ' ')
log "VCF sample count  : $N_VCF_SAMPLES  (expected $TOTAL)"
[[ "$N_VCF_SAMPLES" -eq "$TOTAL" ]] \
    || warn "Sample count mismatch — got $N_VCF_SAMPLES expected $TOTAL"

# =============================================================================
# STEP 1 — EXTRACT SAMPLE NAMES
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
# STEP 2 — COMPUTE ALLELE FREQUENCIES  (single VCF pass)
# =============================================================================
log ""
log "── Step 2: Computing allele frequencies ─────────────────────────────"
log "  Streaming VCF — may take several minutes for large files..."

FREQ_FILE="$OUTDIR/01_allele_frequencies/all_variant_frequencies.tsv"

# Write header first
printf "variant_id\tchrom\tpos\tref\talt\talt_count\tfreq_pct\n" > "$FREQ_FILE"

$AWK -F'\t' '
# ── skip meta-information lines ────────────────────────────────────────
/^##/ { next }

# ── capture sample count from header ───────────────────────────────────
/^#CHROM/ { n_samples = NF - 9; next }

# ── process variant lines ───────────────────────────────────────────────
{
    ref = $4
    alt = $5

    # keep only biallelic SNPs
    if (index(alt, ",") > 0)          next   # multi-allelic
    if (alt == "." || alt == "*")     next   # missing / spanning del
    if (substr(alt, 1, 1) == "<")    next   # symbolic allele
    if (length(ref) != 1)            next   # indel
    if (length(alt) != 1)            next   # indel

    # locate GT within FORMAT  (default position = 1 if GT not labelled)
    gt_pos = 1
    n_fmt  = split($9, fmt_arr, ":")
    for (f = 1; f <= n_fmt; f++) {
        if (fmt_arr[f] == "GT") { gt_pos = f; break }
    }

    # count samples carrying ≥1 ALT allele
    # missing (. or ./.) → REF → count = 0  (SNIPPY convention)
    alt_count = 0
    for (i = 10; i <= NF; i++) {
        n_f = split($i, sf, ":")
        gt  = (n_f >= gt_pos) ? sf[gt_pos] : sf[1]
        gsub(/\|/, "/", gt)              # phased → unphased
        n_al = split(gt, al, "/")
        for (a = 1; a <= n_al; a++) {
            if (al[a] == "1") { alt_count++; break }
        }
    }

    # skip monomorphic sites
    if (alt_count == 0 || alt_count == n_samples) next

    freq_pct = alt_count / n_samples * 100
    vid      = $1 "_" $2 "_" ref "_" alt

    printf "%s\t%s\t%s\t%s\t%s\t%d\t%.6f\n",
        vid, $1, $2, ref, alt, alt_count, freq_pct
}
' "$VCF_FILE" >> "$FREQ_FILE"

N_VARS=$(( $(wc -l < "$FREQ_FILE") - 1 ))
log "  $N_VARS biallelic SNPs → 01_allele_frequencies/all_variant_frequencies.tsv"
[[ "$N_VARS" -gt 0 ]] || die "No variants found — check VCF format and content"

# Frequency distribution summary for QC
log "  Frequency distribution:"
$AWK -F'\t' '
NR == 1 { next }
{
    f = $7 + 0
    if      (f <  1)  bins["<1%"]++
    else if (f <  5)  bins["1-5%"]++
    else if (f < 10)  bins["5-10%"]++
    else if (f < 20)  bins["10-20%"]++
    else if (f < 30)  bins["20-30%"]++
    else if (f < 40)  bins["30-40%"]++
    else if (f < 50)  bins["40-50%"]++
    else if (f < 75)  bins["50-75%"]++
    else              bins[">75%"]++
}
END {
    for (b in bins) printf "    %-10s  %d\n", b, bins[b]
}
' "$FREQ_FILE" | sort | tee -a "$LOG"

# =============================================================================
# STEP 3 — SELECT ONE SNP PER FREQUENCY BIN
# =============================================================================
log ""
log "── Step 3: Selecting one SNP per frequency bin ──────────────────────"

SEL_FILE="$OUTDIR/02_selected_snps/selected_snps.tsv"
printf "target_pct\tvariant_id\tchrom\tpos\tref\talt\tn_carriers\tactual_pct\tdelta_pct\n" \
    > "$SEL_FILE"

# Arrays keyed by target frequency (integer %)
declare -a SEL_TARGETS      # ordered list: 5 10 15 ... 50
declare -A SEL_VID          # target → variant_id
declare -A SEL_NCARR        # target → n_carriers

for (( T=FREQ_MIN; T<=FREQ_MAX; T+=FREQ_STEP )); do

    RESULT=""

    # Try primary tolerance, then doubled tolerance
    for PASS in 1 2; do
        [[ -n "$RESULT" ]] && break
        PASS_TOL=$(awk "BEGIN{printf \"%.4f\", $TOL * $PASS}")

        RESULT=$(
            $AWK -F'\t' -v target="$T" -v tol="$PASS_TOL" '
            NR == 1 { next }
            {
                diff = $7 - target
                if (diff < 0) diff = -diff
                if (diff <= tol) {
                    if (best_vid == "" ||
                        diff < best_diff ||
                        (diff == best_diff && $6 + 0 > best_count)) {
                        best_vid   = $1
                        best_diff  = diff
                        best_count = $6 + 0
                        best_line  = $0
                    }
                }
            }
            END {
                if (best_vid == "") print "NOTFOUND"
                else               print best_line
            }
            ' "$FREQ_FILE"
        )

        if [[ "$RESULT" == "NOTFOUND" ]]; then
            RESULT=""
            [[ $PASS -eq 1 ]] && \
                warn "  ${T}%: no SNP within ±${TOL}pp — widening to ±$(awk "BEGIN{print $TOL*2}")pp"
        fi
    done

    if [[ -z "$RESULT" ]]; then
        warn "  ${T}%: NO SNP FOUND — skipping this bin"
        continue
    fi

    VID=$(    echo "$RESULT" | cut -f1 | tr -d ' ')
    N_CARR=$( echo "$RESULT" | cut -f6 | tr -d ' ')
    ACT=$(    echo "$RESULT" | cut -f7 | tr -d ' ')
    DELTA=$(  awk "BEGIN{d=$ACT-$T; if(d<0)d=-d; printf \"%.4f\",d}")

    SEL_TARGETS+=("$T")
    SEL_VID[$T]="$VID"
    SEL_NCARR[$T]="$N_CARR"

    # Append row to selected_snps.tsv
    printf "%s\t%s\t%s\n" "$T" "$RESULT" "$DELTA" \
        | $AWK -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                         $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
        >> "$SEL_FILE"

    log "  ${T}%  →  $VID  actual=${ACT}%  n_carriers=${N_CARR}"
done

N_SEL=${#SEL_TARGETS[@]}
log "  $N_SEL / $N_BINS SNPs selected → 02_selected_snps/selected_snps.tsv"
[[ "$N_SEL" -gt 0 ]] || die "No SNPs selected — check VCF content and parameters"

# =============================================================================
# STEP 4 — EXTRACT CARRIER INDICES  (single second VCF pass)
# =============================================================================
log ""
log "── Step 4: Extracting carrier indices (one VCF pass for all SNPs) ───"

# Write the selected variant IDs to a file for awk lookup
VIDS_FILE="$OUTDIR/tmp/selected_vids.txt"
for T in "${SEL_TARGETS[@]}"; do
    echo "${SEL_VID[$T]}"
done > "$VIDS_FILE"

N_TARGETS=$(wc -l < "$VIDS_FILE" | tr -d ' ')

$AWK -F'\t' \
    -v vids_file="$VIDS_FILE" \
    -v out_dir="$OUTDIR/tmp/carriers" \
    -v n_targets="$N_TARGETS" \
'
BEGIN {
    while ((getline line < vids_file) > 0) {
        target_vids[line] = 1
    }
    close(vids_file)
    found = 0
}
/^##/ { next }
/^#CHROM/ { n_samples = NF - 9; next }
{
    ref = $4; alt = $5

    # same filters as step 2
    if (index(alt, ",") > 0)       next
    if (alt == "." || alt == "*")  next
    if (substr(alt,1,1) == "<")    next
    if (length(ref) != 1)         next
    if (length(alt) != 1)         next

    vid = $1 "_" $2 "_" ref "_" alt
    if (!(vid in target_vids)) next

    # locate GT field
    gt_pos = 1
    n_fmt  = split($9, fmt_arr, ":")
    for (f = 1; f <= n_fmt; f++) {
        if (fmt_arr[f] == "GT") { gt_pos = f; break }
    }

    # write 0-based sample indices of ALT carriers to per-SNP file
    outfile = out_dir "/" vid ".txt"
    for (i = 10; i <= NF; i++) {
        n_f = split($i, sf, ":")
        gt  = (n_f >= gt_pos) ? sf[gt_pos] : sf[1]
        gsub(/\|/, "/", gt)
        n_al = split(gt, al, "/")
        for (a = 1; a <= n_al; a++) {
            if (al[a] == "1") {
                print (i - 10) >> outfile   # 0-based index
                break
            }
        }
    }
    close(outfile)

    found++
    if (found >= n_targets) exit    # all SNPs found — stop reading VCF
}
' "$VCF_FILE"

# Verify carrier files were created and match expected counts
for T in "${SEL_TARGETS[@]}"; do
    VID="${SEL_VID[$T]}"
    CF="$OUTDIR/tmp/carriers/${VID}.txt"
    if [[ -f "$CF" ]]; then
        NC=$(wc -l < "$CF" | tr -d ' ')
        EXP="${SEL_NCARR[$T]}"
        STATUS="✓"
        [[ "$NC" -ne "$EXP" ]] && STATUS="⚠ count mismatch: got $NC expected $EXP"
        log "  ${T}%  $VID  →  $NC carriers  $STATUS"
    else
        warn "  Carrier file NOT created for $VID — this SNP will be skipped"
    fi
done

# =============================================================================
# STEP 5 — GENERATE PHENOTYPE FILES
# =============================================================================
log ""
log "── Step 5: Generating phenotype files ───────────────────────────────"
log "  $N_SEL SNPs × 4 ratios × $N_SIMS sims = up to $EXPECTED_FILES files"

N_WRITTEN=0
N_SKIP=0

# Ratio definitions  (EOD_parts LOD_parts)
declare -A RATIOS=(
    ["EOD1_LOD3"]="1 3"
    ["EOD1_LOD2"]="1 2"
    ["EOD2_LOD1"]="2 1"
    ["EOD3_LOD1"]="3 1"
)
RATIO_ORDER=("EOD1_LOD3" "EOD1_LOD2" "EOD2_LOD1" "EOD3_LOD1")

for T in "${SEL_TARGETS[@]}"; do

    VID="${SEL_VID[$T]}"
    FREQ_LABEL=$(printf "freq_%02dpct" "$T")
    CF="$OUTDIR/tmp/carriers/${VID}.txt"

    if [[ ! -f "$CF" ]]; then
        warn "No carrier file for $VID — skipping"
        continue
    fi

    N_CARR=$(wc -l < "$CF" | tr -d ' ')

    log ""
    log "  [$FREQ_LABEL]  $VID  n_carriers=$N_CARR"

    for RNAME in "${RATIO_ORDER[@]}"; do

        IFS=' ' read -r EOD_P LOD_P <<< "${RATIOS[$RNAME]}"
        DENOM=$(( EOD_P + LOD_P ))

        # Carrier splits — use awk for correct rounding
        N_EOD_C=$( awk "BEGIN{printf \"%d\", int($N_CARR * $EOD_P / $DENOM + 0.5)}")
        N_LOD_C=$(( N_CARR - N_EOD_C ))
        N_EOD_NC=$(( TOT_EOD - N_EOD_C ))
        N_LOD_NC=$(( TOT_LOD - N_LOD_C ))

        # Feasibility check (should never fail with 1:3 1:2 2:1 3:1 at ≤50%)
        if [[ "$N_EOD_NC" -lt 0 || "$N_LOD_NC" -lt 0 ]]; then
            warn "    SKIP $RNAME: infeasible (n_eod_nc=$N_EOD_NC n_lod_nc=$N_LOD_NC)"
            N_SKIP=$(( N_SKIP + 1 ))
            continue
        fi

        log "    $RNAME  carriers   → EOD=$N_EOD_C  LOD=$N_LOD_C"
        log "           non-carriers → EOD=$N_EOD_NC LOD=$N_LOD_NC"

        RATIO_DIR="$OUTDIR/03_phenotype_files/$FREQ_LABEL/$VID/ratio_${RNAME}"
        mkdir -p "$RATIO_DIR"

        # Deterministic base seed for this (SNP, ratio) combo
        COMBO_HASH=$(printf "%s|%s" "$VID" "$RNAME" | cksum | cut -d' ' -f1)
        COMBO_SEED=$(( (SEED + COMBO_HASH) % 16777216 ))

        # ── Generate ALL N_SIMS files in ONE awk call ────────────────────
        # All data is loaded once; each simulation uses its own shuffled copy
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

            # ── Load sample names (0-indexed) ─────────────────────────────
            ns = 0
            while ((getline sname < samples_file) > 0) {
                samples[ns] = sname
                ns++
            }
            close(samples_file)

            # ── Load carrier indices ──────────────────────────────────────
            nc = 0
            while ((getline idx < carriers_file) > 0) {
                carriers[nc]      = idx + 0
                is_carrier[idx+0] = 1
                nc++
            }
            close(carriers_file)

            # ── Build non-carrier list ────────────────────────────────────
            nnc = 0
            for (j = 0; j < ns; j++) {
                if (!is_carrier[j]) {
                    noncarriers[nnc] = j
                    nnc++
                }
            }

            # ── Generate one file per simulation ─────────────────────────
            for (sim = 1; sim <= n_sims; sim++) {

                # Two independent seeds per simulation:
                #   carrier_seed    → shuffles carrier list
                #   noncarrier_seed → shuffles non-carrier list
                # Using different large-prime multipliers ensures independence
                carrier_seed    = (combo_seed + sim * 999983)  % 16777216
                noncarrier_seed = (combo_seed + sim * 1000003) % 16777216

                # ── Fisher-Yates shuffle: carriers ───────────────────────
                srand(carrier_seed)
                for (j = 0; j < nc; j++) ct[j] = carriers[j]
                for (j = nc - 1; j > 0; j--) {
                    k = int(rand() * (j + 1))
                    tmp = ct[j]; ct[j] = ct[k]; ct[k] = tmp
                }

                # ── Fisher-Yates shuffle: non-carriers ───────────────────
                srand(noncarrier_seed)
                for (j = 0; j < nnc; j++) nct[j] = noncarriers[j]
                for (j = nnc - 1; j > 0; j--) {
                    k = int(rand() * (j + 1))
                    tmp = nct[j]; nct[j] = nct[k]; nct[k] = tmp
                }

                # ── Assign phenotypes ─────────────────────────────────────
                # first n_eod_c  shuffled carriers     → EOD (0)
                # rest of carriers                     → LOD (1)
                # first n_eod_nc shuffled non-carriers → EOD (0)
                # rest of non-carriers                 → LOD (1)
                for (k in pheno) delete pheno[k]   # reset between sims
                for (j = 0; j < n_eod_c;  j++) pheno[ct[j]]  = 0
                for (j = n_eod_c; j < nc; j++) pheno[ct[j]]  = 1
                for (j = 0; j < n_eod_nc;  j++) pheno[nct[j]] = 0
                for (j = n_eod_nc; j < nnc; j++) pheno[nct[j]] = 1

                # ── Hard verification ─────────────────────────────────────
                eod_n = 0; lod_n = 0
                for (j = 0; j < ns; j++) {
                    if (pheno[j] == 0) eod_n++
                    else               lod_n++
                }
                if (eod_n != tot_eod || lod_n != tot_lod) {
                    printf "ASSERTION FAILED: sim=%d EOD=%d LOD=%d (expected %d %d)\n",
                        sim, eod_n, lod_n, tot_eod, tot_lod > "/dev/stderr"
                    exit 1
                }

                # ── Write phenotype file ──────────────────────────────────
                outfile = ratio_dir "/" sprintf("sim_%03d.pheno", sim)

                if (fmt == "tab") {
                    print "sample\tphenotype" > outfile
                    for (j = 0; j < ns; j++) {
                        print samples[j] "\t" pheno[j] > outfile
                    }
                } else {
                    print "FID\tIID\tPHENOTYPE" > outfile
                    for (j = 0; j < ns; j++) {
                        print samples[j] "\t" samples[j] "\t" (pheno[j]+1) > outfile
                    }
                }
                close(outfile)

            }  # end sim loop
        }
        ' /dev/null

        N_WRITTEN=$(( N_WRITTEN + N_SIMS ))
        log "    → $N_SIMS files written to ratio_${RNAME}/"

    done  # ratio loop
done  # target loop

# =============================================================================
# STEP 6 — SPOT CHECK  (10 random files)
# =============================================================================
log ""
log "── Step 6: Spot-checking 10 random phenotype files ─────────────────"

PASS=0
FAIL=0

while IFS= read -r PFILE; do
    if [[ ! -f "$PFILE" ]]; then
        log "  ✗  MISSING: $PFILE"
        FAIL=$(( FAIL + 1 ))
        continue
    fi

    READ_EOD=$(awk -F'\t' 'NR>1 && $2==0{c++} END{print c+0}' "$PFILE")
    READ_LOD=$(awk -F'\t' 'NR>1 && $2==1{c++} END{print c+0}' "$PFILE")
    LABEL=$(basename "$(dirname "$PFILE")")/$(basename "$PFILE")

    if [[ "$READ_EOD" -eq "$TOT_EOD" && "$READ_LOD" -eq "$TOT_LOD" ]]; then
        log "  ✓  $LABEL   EOD=$READ_EOD  LOD=$READ_LOD"
        PASS=$(( PASS + 1 ))
    else
        log "  ✗  $LABEL   GOT EOD=$READ_EOD LOD=$READ_LOD   WANT $TOT_EOD / $TOT_LOD"
        FAIL=$(( FAIL + 1 ))
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
log "  Phenotype files  : $FINAL_N  (expected $EXPECTED_FILES)"
log "  Skipped combos   : $N_SKIP"
log "  Spot check       : $PASS passed  /  $FAIL failed"
log ""
log "  Output layout:"
log "    $OUTDIR/"
log "    ├── 01_allele_frequencies/all_variant_frequencies.tsv"
log "    ├── 02_selected_snps/selected_snps.tsv"
log "    └── 03_phenotype_files/"
log "        ├── freq_05pct/<SNP_ID>/ratio_EOD1_LOD3/sim_001.pheno"
log "        │                                      /sim_020.pheno"
log "        │                      /ratio_EOD1_LOD2/"
log "        │                      /ratio_EOD2_LOD1/"
log "        │                      /ratio_EOD3_LOD1/"
log "        └── freq_50pct/..."
log ""
log "  Each .pheno file contains $TOTAL lines + 1 header:"
log "    sample      phenotype"
log "    isolate_01  0          ← EOD"
log "    isolate_02  1          ← LOD"

if [[ "$FAIL" -eq 0 ]]; then
    log ""
    log "  ✅  All spot-checks PASSED — ready to use"
else
    log ""
    log "  ❌  $FAIL failures — review $LOG"
fi

log ""
log "======================================================================"
log "  DONE  $(date)"
log "======================================================================"