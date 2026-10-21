#!/bin/bash
#SBATCH -J homoplasy
#SBATCH -p fasse
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o homoplasy_%j.out
#SBATCH -e homoplasy_%j.err

set -euo pipefail

# ═══════════════════════════════════════════════════════════════
#  PATHS  —  edit if needed
# ═══════════════════════════════════════════════════════════════
VCF="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/Nextflow_pipeline_EOD_LOD_VLOD_CDC/11.SNIPPY_MULTI/merged.vcf"
TREE="/n/hanage_gbs_gwas_l3/Lab/GBS-GWAS-pipeline/Nextflow_pipeline_EOD_LOD_VLOD_CDC/12.Phylogeny/12.1.FastTree_DP/EOD_LOD_VLOD_CDC_FastTree_DP.tre"
OUT="${SLURM_SUBMIT_DIR}/homoplasy_results.tsv"
ENV_NAME="homoplasy_env"

# ═══════════════════════════════════════════════════════════════
#  CONDA SETUP
# ═══════════════════════════════════════════════════════════════
echo "[CONDA] Initialising conda..."

# Try common conda locations on FASRC / generic clusters
CONDA_PATHS=(
    "/n/sw/Miniforge3-26.1.0-0/etc/profile.d/conda.sh"
    "${HOME}/miniconda3/etc/profile.d/conda.sh"
    "${HOME}/anaconda3/etc/profile.d/conda.sh"
    "/opt/conda/etc/profile.d/conda.sh"
)

CONDA_FOUND=0
for p in "${CONDA_PATHS[@]}"; do
    if [[ -f "$p" ]]; then
        source "$p"
        CONDA_FOUND=1
        echo "       Sourced: $p"
        break
    fi
done

if [[ $CONDA_FOUND -eq 0 ]]; then
    echo "ERROR: conda not found. Install miniconda first:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3"
    exit 1
fi

# ── Create environment (only if it does not already exist) ──────
if conda env list | grep -q "^${ENV_NAME}[[:space:]]"; then
    echo "[CONDA] Environment '${ENV_NAME}' already exists — skipping creation"
else
    echo "[CONDA] Creating environment '${ENV_NAME}'..."
    conda create -y -n "${ENV_NAME}" \
        -c etetoolkit \
        -c conda-forge \
        python=3.10 \
        ete3 \
        numpy
    echo "[CONDA] Environment created."
fi

conda activate "${ENV_NAME}"
echo "[CONDA] Active env: $(conda info --envs | grep '*' | awk '{print $1}')"
python --version

# ═══════════════════════════════════════════════════════════════
#  WRITE PYTHON SCRIPT TO A TEMP FILE VIA HEREDOC
# ═══════════════════════════════════════════════════════════════
PYTHON_SCRIPT=$(mktemp /tmp/homoplasy_XXXXXX.py)
echo "[SCRIPT] Temp script written to: ${PYTHON_SCRIPT}"

cat > "${PYTHON_SCRIPT}" << 'PYEOF'
#!/usr/bin/env python3
"""
Per-SNP Homoplasy Index via vectorised Fitch parsimony.

Metrics
-------
  parsimony score (l) : min # state changes inferred on tree
  min changes    (s)  : (# distinct alleles) - 1
  CI = s / l          : Consistency Index  (1.0 = no homoplasy)
  HI = 1 - CI         : Homoplasy Index    (0.0 = no homoplasy)
"""

import sys, os, argparse
import numpy as np


# ───────────────────────────────────────────────────────────────
#  VCF READER
# ───────────────────────────────────────────────────────────────
def read_vcf(path):
    samples, meta, rows = [], [], []
    print(f"[VCF]  Reading: {path}")
    with open(path) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]
                n = len(samples)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            fmt = fields[8].split(':')
            if 'GT' not in fmt:
                continue
            gi  = fmt.index('GT')
            row = np.full(n, -1, np.int8)
            for i in range(n):
                if 9 + i >= len(fields):
                    continue
                gt_raw = fields[9 + i].split(':')[gi]
                sep    = '/' if '/' in gt_raw else ('|' if '|' in gt_raw else None)
                allele = gt_raw.split(sep)[0] if sep else gt_raw
                if allele not in ('.', ''):
                    try:
                        row[i] = int(allele)
                    except ValueError:
                        pass
            meta.append((fields[0], int(fields[1]),
                          fields[2], fields[3], fields[4]))
            rows.append(row)

    gt = np.vstack(rows) if rows else np.empty((0, n), np.int8)
    print(f"       {len(samples)} samples  |  {len(meta)} SNPs")
    return samples, meta, gt


# ───────────────────────────────────────────────────────────────
#  TREE LOADER
# ───────────────────────────────────────────────────────────────
def load_tree(path):
    try:
        from ete3 import Tree
    except ImportError:
        sys.exit("ERROR: ete3 not installed in this environment.")

    print(f"[TREE] Reading: {path}")
    for fmt in (1, 0, 5):
        try:
            t = Tree(path, format=fmt)
            print(f"       {len(t.get_leaves())} tips  (ete3 format={fmt})")
            return t
        except Exception:
            pass
    sys.exit(f"ERROR: Could not parse tree: {path}")


# ───────────────────────────────────────────────────────────────
#  FITCH PARSIMONY  (vectorised across all SNP sites in a chunk)
# ───────────────────────────────────────────────────────────────
def build_traversal(tree, leaf_order):
    """One-time indexing of tree nodes for fast array ops."""
    leaf_col = {name: i for i, name in enumerate(leaf_order)}
    nid, ctr = {}, [0]

    def _label(nd):
        nid[id(nd)] = ctr[0]; ctr[0] += 1
        for ch in nd.children:
            _label(ch)
    _label(tree)

    traversal = []
    for nd in tree.traverse('postorder'):
        i = nid[id(nd)]
        if nd.is_leaf():
            traversal.append((i, True, leaf_col.get(nd.name, -1), []))
        else:
            traversal.append((i, False, -1,
                               [nid[id(c)] for c in nd.children]))
    return traversal, ctr[0]


def fitch_chunk(traversal, n_nodes, gts):
    """
    Vectorised Fitch over a batch of SNP sites.
    State sets encoded as uint32 bitmasks:
        bit k set  ⟺  allele k in the set
        0xFFFFFFFF  =  missing (all alleles allowed)
    """
    AMBIG  = np.uint32(0xFFFFFFFF)
    n_sites = gts.shape[0]
    S      = np.zeros((n_nodes, n_sites), dtype=np.uint32)
    scores = np.zeros(n_sites,            dtype=np.int32)

    for nid, is_leaf, col, ch_ids in traversal:
        if is_leaf:
            if col < 0:
                S[nid] = AMBIG
            else:
                a      = gts[:, col].astype(np.int32)
                a_safe = np.clip(a, 0, 30).astype(np.uint32)
                S[nid] = np.where(a < 0, AMBIG,
                                  (np.uint32(1) << a_safe).astype(np.uint32))
        else:
            if not ch_ids:
                S[nid] = AMBIG
                continue
            inter = S[ch_ids[0]].copy()
            union = S[ch_ids[0]].copy()
            for cid in ch_ids[1:]:
                inter &= S[cid]
                union |= S[cid]
            no_inter  = (inter == 0)
            scores   += no_inter.astype(np.int32)
            S[nid]    = np.where(no_inter, union, inter).astype(np.uint32)

    return scores


# ───────────────────────────────────────────────────────────────
#  MAIN
# ───────────────────────────────────────────────────────────────
def run(vcf_path, tree_path, out_path, chunk_size):

    samples, snp_meta, gt = read_vcf(vcf_path)
    tree                   = load_tree(tree_path)

    # ── sample name matching ────────────────────────────────────
    tips    = [lf.name for lf in tree.get_leaves()]
    common  = [s for s in tips if s in set(samples)]

    print(f"\n[MATCH] VCF={len(samples)}  tree={len(tips)}  common={len(common)}")
    if not common:
        print("  ERROR: no overlapping sample names!")
        print(f"  Tree (first 5) : {tips[:5]}")
        print(f"  VCF  (first 5) : {samples[:5]}")
        sys.exit(1)

    if len(common) < len(tips):
        print(f"  WARNING: {len(tips)-len(common)} tree tips absent from VCF — pruning")
        tree.prune(common, preserve_branch_length=True)

    leaf_order = [lf.name for lf in tree.get_leaves()]
    sidx       = {s: i for i, s in enumerate(samples)}
    gt_sub     = gt[:, [sidx[s] for s in leaf_order]]

    # ── build traversal index once ──────────────────────────────
    print("\n[INDEX] Building tree traversal...")
    traversal, n_nodes = build_traversal(tree, leaf_order)

    # ── chunked Fitch ───────────────────────────────────────────
    n_snps = gt_sub.shape[0]
    ps     = np.zeros(n_snps, np.int32)
    w      = len(str(n_snps))

    print(f"[FITCH] {n_snps} SNPs  |  chunk={chunk_size}")
    for s0 in range(0, n_snps, chunk_size):
        s1       = min(s0 + chunk_size, n_snps)
        ps[s0:s1] = fitch_chunk(traversal, n_nodes, gt_sub[s0:s1])
        print(f"        {s1:{w}d} / {n_snps}  ({100*s1/n_snps:.0f}%)", end='\r')
    print()

    # ── allele counts and metrics ───────────────────────────────
    print("[STATS] Computing CI / HI...")
    n_al  = np.array([len(np.unique(r[r >= 0])) for r in gt_sub], dtype=np.int32)
    s_min = np.maximum(n_al - 1, 0)

    with np.errstate(divide='ignore', invalid='ignore'):
        ci = np.where((ps > 0) & (n_al >= 2),
                      s_min.astype(float) / ps.astype(float), np.nan)
    hi = np.where(np.isfinite(ci), 1.0 - ci, np.nan)

    # ── write TSV ───────────────────────────────────────────────
    print(f"[OUT]   Writing: {out_path}")
    with open(out_path, 'w') as fh:
        fh.write('CHROM\tPOS\tID\tREF\tALT\t'
                 'n_alleles\tparsimony_score\tmin_changes\tCI\tHI\n')
        for i, (chrom, pos, sid, ref, alt) in enumerate(snp_meta):
            ci_s = f'{ci[i]:.6f}' if np.isfinite(ci[i]) else 'NA'
            hi_s = f'{hi[i]:.6f}' if np.isfinite(hi[i]) else 'NA'
            fh.write(f'{chrom}\t{pos}\t{sid}\t{ref}\t{alt}\t'
                     f'{n_al[i]}\t{ps[i]}\t{s_min[i]}\t{ci_s}\t{hi_s}\n')

    # ── summary stats ───────────────────────────────────────────
    poly = n_al >= 2
    vh   = hi[np.isfinite(hi) & poly]
    print(f"\n{'─'*50}")
    print(f"Total SNPs              : {n_snps}")
    print(f"Polymorphic SNPs        : {int(poly.sum())}")
    print(f"Monomorphic / no data   : {int((~poly).sum())}")
    if len(vh):
        no_h  = int((vh == 0).sum())
        has_h = int((vh  > 0).sum())
        print(f"HI = 0  (no homoplasy)  : {no_h}  ({100*no_h/len(vh):.1f}%)")
        print(f"HI > 0  (homoplastic)   : {has_h}  ({100*has_h/len(vh):.1f}%)")
        print(f"Mean HI                 : {np.mean(vh):.4f}")
        print(f"Median HI               : {np.median(vh):.4f}")
        print(f"Max HI                  : {np.max(vh):.4f}")
    print('─'*50)
    print("Done.")


def main():
    ap = argparse.ArgumentParser(description='Per-SNP Homoplasy Index (Fitch parsimony)')
    ap.add_argument('--vcf',   required=True)
    ap.add_argument('--tree',  required=True)
    ap.add_argument('--out',   required=True)
    ap.add_argument('--chunk', type=int, default=5000)
    args = ap.parse_args()
    for p in (args.vcf, args.tree):
        if not os.path.exists(p):
            sys.exit(f"ERROR: file not found: {p}")
    run(args.vcf, args.tree, args.out, args.chunk)

if __name__ == '__main__':
    main()
PYEOF

# ═══════════════════════════════════════════════════════════════
#  RUN
# ═══════════════════════════════════════════════════════════════
echo ""
echo "[RUN]  Starting homoplasy calculation..."
echo "       VCF  : ${VCF}"
echo "       TREE : ${TREE}"
echo "       OUT  : ${OUT}"
echo ""

python "${PYTHON_SCRIPT}" \
    --vcf   "${VCF}"  \
    --tree  "${TREE}" \
    --out   "${OUT}"  \
    --chunk 5000

# ── cleanup ─────────────────────────────────────────────────────
rm -f "${PYTHON_SCRIPT}"

echo ""
echo "Results written to: ${OUT}"
echo "Elapsed time: ${SECONDS}s"
