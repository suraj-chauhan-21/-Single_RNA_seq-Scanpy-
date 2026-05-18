"""
scripts/02_normalization.py
────────────────────────────
Step 2: Normalize counts, log-transform, select highly variable genes (HVGs),
and scale data for downstream analysis.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scanpy as sc
import config as cfg
from utils import load_adata, save_adata, ensure_dirs, plot_hvg


def run_normalization():
    sc.settings.verbosity = cfg.VERBOSITY
    sc.settings.set_figure_params(dpi=cfg.FIGURE["dpi"], facecolor="white")

    fig_dir = os.path.join(cfg.FIGURES_DIR, "02_normalization")
    ensure_dirs(fig_dir)

    # ── Load ──────────────────────────────────────────────────
    print("\n=== Step 2: Normalization & Feature Selection ===\n")
    adata = load_adata(cfg.ADATA_PROCESSED_PATH)

    # ── Store raw counts ──────────────────────────────────────
    # Save raw (filtered) counts before normalization
    adata.raw = adata

    # ── Normalization ─────────────────────────────────────────
    print(f"[norm] Normalizing to {cfg.NORMALIZATION['target_sum']:.0f} counts per cell...")
    sc.pp.normalize_total(adata, target_sum=cfg.NORMALIZATION["target_sum"])

    # Log1p transformation
    if cfg.NORMALIZATION["log1p"]:
        print("[norm] Applying log1p transformation...")
        sc.pp.log1p(adata)

    # ── Highly Variable Genes (HVGs) ──────────────────────────
    print(f"[norm] Selecting top {cfg.NORMALIZATION['n_top_genes']} highly variable genes...")

    flavor = cfg.NORMALIZATION["hvg_flavor"]
    if flavor == "seurat_v3":
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=cfg.NORMALIZATION["n_top_genes"],
            flavor="seurat_v3",
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            min_mean=cfg.NORMALIZATION["hvg_min_mean"],
            max_mean=cfg.NORMALIZATION["hvg_max_mean"],
            min_disp=cfg.NORMALIZATION["hvg_min_disp"],
            flavor=flavor,
        )

    n_hvg = adata.var["highly_variable"].sum()
    print(f"[norm] Found {n_hvg} highly variable genes.")

    # Plot HVGs
    plot_hvg(adata, save_path=os.path.join(fig_dir, "hvg_plot.png"))

    # ── Subset to HVGs ────────────────────────────────────────
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"[norm] Data after HVG subset: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── Regress out confounders (optional) ────────────────────
    # Regress out total counts and % mitochondrial reads to reduce noise
    print("[norm] Regressing out confounders (total_counts, pct_counts_mt)...")
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

    # ── Scale features ────────────────────────────────────────
    # Zero-mean, unit-variance; clip at max_value=10 to limit outlier influence
    print("[norm] Scaling to unit variance (max_value=10)...")
    sc.pp.scale(adata, max_value=10)

    # ── Save ──────────────────────────────────────────────────
    save_adata(adata, cfg.ADATA_NORMALIZED_PATH)
    print("\n=== Step 2 Complete ===\n")


if __name__ == "__main__":
    run_normalization()
