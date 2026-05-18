"""
scripts/03_dim_reduction.py
────────────────────────────
Step 3: PCA → Neighborhood Graph → UMAP & t-SNE
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scanpy as sc
import matplotlib.pyplot as plt
import config as cfg
from utils import load_adata, save_adata, ensure_dirs, plot_pca_variance, plot_umap, plot_tsne


def run_dim_reduction():
    sc.settings.verbosity = cfg.VERBOSITY
    sc.settings.set_figure_params(dpi=cfg.FIGURE["dpi"], facecolor="white")

    fig_dir = os.path.join(cfg.FIGURES_DIR, "03_dim_reduction")
    ensure_dirs(fig_dir)

    # ── Load ──────────────────────────────────────────────────
    print("\n=== Step 3: Dimensionality Reduction ===\n")
    adata = load_adata(cfg.ADATA_NORMALIZED_PATH)

    dr = cfg.DIM_REDUCTION

    # ── PCA ───────────────────────────────────────────────────
    print(f"[dr] Running PCA with {dr['n_comps_pca']} components...")
    sc.tl.pca(
        adata,
        n_comps=dr["n_comps_pca"],
        svd_solver="arpack",
        random_state=dr["random_state"],
    )

    # Scree / elbow plot
    plot_pca_variance(
        adata,
        n_pcs=dr["n_comps_pca"],
        save_path=os.path.join(fig_dir, "pca_variance"),
    )

    # PCA scatter colored by QC metrics
    sc.pl.pca(
        adata,
        color=["total_counts", "pct_counts_mt"],
        show=False,
    )
    plt.savefig(os.path.join(fig_dir, "pca_qc_coloring.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # ── Neighborhood Graph ────────────────────────────────────
    print(f"[dr] Building k-NN graph (k={dr['n_neighbors']}, n_pcs={dr['n_pcs']})...")
    sc.pp.neighbors(
        adata,
        n_neighbors=dr["n_neighbors"],
        n_pcs=dr["n_pcs"],
        random_state=dr["random_state"],
    )

    # ── UMAP ──────────────────────────────────────────────────
    print("[dr] Computing UMAP...")
    sc.tl.umap(
        adata,
        min_dist=dr["umap_min_dist"],
        spread=dr["umap_spread"],
        random_state=dr["random_state"],
    )

    plot_umap(
        adata,
        color=["total_counts", "pct_counts_mt", "n_genes_by_counts"],
        save_path=os.path.join(fig_dir, "umap_qc.png"),
        ncols=3,
    )

    # ── t-SNE ─────────────────────────────────────────────────
    print("[dr] Computing t-SNE...")
    sc.tl.tsne(
        adata,
        n_pcs=dr["n_pcs"],
        perplexity=dr["tsne_perplexity"],
        random_state=dr["random_state"],
    )

    plot_tsne(
        adata,
        color=["total_counts", "pct_counts_mt"],
        save_path=os.path.join(fig_dir, "tsne_qc.png"),
    )

    # ── Save ──────────────────────────────────────────────────
    save_adata(adata, cfg.ADATA_REDUCED_PATH)
    print("\n=== Step 3 Complete ===\n")


if __name__ == "__main__":
    run_dim_reduction()
