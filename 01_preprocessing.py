"""
scripts/01_preprocessing.py
─────────────────────────────
Step 1: Load raw data, compute QC metrics, filter low-quality
cells and genes, and save the preprocessed AnnData.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scanpy as sc
import config as cfg
from utils import (
    load_data, save_adata, ensure_dirs,
    compute_qc_metrics, filter_cells_and_genes, detect_doublets,
    plot_qc_metrics, plot_qc_scatter, plot_highest_expressed_genes,
)


def run_preprocessing():
    # ── Setup ────────────────────────────────────────────────
    sc.settings.verbosity = cfg.VERBOSITY
    sc.settings.set_figure_params(dpi=cfg.FIGURE["dpi"], facecolor="white")

    ensure_dirs(cfg.RESULTS_DIR, cfg.FIGURES_DIR)
    fig_dir = os.path.join(cfg.FIGURES_DIR, "01_qc")
    ensure_dirs(fig_dir)

    # ── Load Data ────────────────────────────────────────────
    print("\n=== Step 1: Preprocessing & Quality Control ===\n")
    cfg_dict = {k: getattr(cfg, k) for k in dir(cfg) if not k.startswith("_")}
    adata = load_data(cfg_dict)

    # Make gene names unique
    adata.var_names_make_unique()

    print(f"\nRaw data: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── QC Metrics ───────────────────────────────────────────
    adata = compute_qc_metrics(adata, mito_prefix=cfg.QC["mito_prefix"])

    # Plot QC before filtering
    plot_qc_metrics(adata, save_path=os.path.join(fig_dir, "qc_violin_before"))
    plot_qc_scatter(adata, save_path=os.path.join(fig_dir, "qc_scatter_before"))
    plot_highest_expressed_genes(
        adata, n_top=20,
        save_path=os.path.join(fig_dir, "highest_expressed_genes.png")
    )

    # ── Filter Cells & Genes ─────────────────────────────────
    adata = filter_cells_and_genes(adata, cfg.QC)

    # Plot QC after filtering
    plot_qc_metrics(adata, save_path=os.path.join(fig_dir, "qc_violin_after"))
    plot_qc_scatter(adata, save_path=os.path.join(fig_dir, "qc_scatter_after"))

    # ── (Optional) Doublet Detection ─────────────────────────
    adata = detect_doublets(adata)

    # Remove predicted doublets if column exists
    if "predicted_doublet" in adata.obs.columns:
        n_before = adata.n_obs
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        print(f"[qc] Removed {n_before - adata.n_obs} doublets.")

    # ── Save ─────────────────────────────────────────────────
    print(f"\nAfter QC: {adata.n_obs} cells × {adata.n_vars} genes")
    save_adata(adata, cfg.ADATA_PROCESSED_PATH)
    print("\n=== Step 1 Complete ===\n")


if __name__ == "__main__":
    run_preprocessing()
