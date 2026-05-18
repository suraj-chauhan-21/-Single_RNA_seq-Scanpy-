"""
utils/plot_utils.py — Reusable plotting functions for scRNA-seq analysis.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import scanpy as sc
import anndata as ad
from typing import Optional, List, Tuple


# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────

def save_fig(fig: plt.Figure, path: str, dpi: int = 150, fmt: str = "png") -> None:
    """Save a matplotlib figure, creating parent directory if needed."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    full_path = f"{path}.{fmt}" if not path.endswith(f".{fmt}") else path
    fig.savefig(full_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] Saved → {full_path}")


# ──────────────────────────────────────────────
# QC Plots
# ──────────────────────────────────────────────

def plot_qc_metrics(
    adata: ad.AnnData,
    save_path: Optional[str] = None,
    figsize: Tuple = (14, 4),
) -> None:
    """
    Plot violin plots of key QC metrics:
      - n_genes_by_counts
      - total_counts
      - pct_counts_mt
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    labels  = ["Genes per Cell", "Total UMI Counts", "% Mitochondrial"]

    for ax, metric, label in zip(axes, metrics, labels):
        data = adata.obs[metric]
        ax.violinplot(data, showmedians=True)
        ax.set_title(label, fontsize=12)
        ax.set_ylabel(label)
        ax.set_xticks([])
        sns.despine(ax=ax)

    fig.suptitle("Quality Control Metrics", fontsize=14, y=1.02)
    plt.tight_layout()

    if save_path:
        save_fig(fig, save_path)
    else:
        plt.show()


def plot_qc_scatter(
    adata: ad.AnnData,
    save_path: Optional[str] = None,
    figsize: Tuple = (12, 4),
) -> None:
    """
    Scatter plots to identify low-quality cells:
      - total_counts vs n_genes_by_counts
      - total_counts vs pct_counts_mt
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    axes[0].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        s=1, alpha=0.3, c="steelblue"
    )
    axes[0].set_xlabel("Total UMI Counts")
    axes[0].set_ylabel("Number of Genes")
    axes[0].set_title("UMI Counts vs Genes Detected")

    axes[1].scatter(
        adata.obs["total_counts"],
        adata.obs["pct_counts_mt"],
        s=1, alpha=0.3, c="tomato"
    )
    axes[1].set_xlabel("Total UMI Counts")
    axes[1].set_ylabel("% Mitochondrial")
    axes[1].set_title("UMI Counts vs % Mitochondrial")

    for ax in axes:
        sns.despine(ax=ax)

    plt.tight_layout()
    if save_path:
        save_fig(fig, save_path)
    else:
        plt.show()


def plot_highest_expressed_genes(
    adata: ad.AnnData,
    n_top: int = 20,
    save_path: Optional[str] = None,
) -> None:
    """Plot top N most highly expressed genes."""
    sc.pl.highest_expr_genes(adata, n_top=n_top, show=False)
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()


# ──────────────────────────────────────────────
# HVG Plot
# ──────────────────────────────────────────────

def plot_hvg(
    adata: ad.AnnData,
    save_path: Optional[str] = None,
) -> None:
    """Plot dispersion vs mean expression highlighting HVGs."""
    sc.pl.highly_variable_genes(adata, show=False)
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"[plot] Saved → {save_path}")


# ──────────────────────────────────────────────
# PCA Plots
# ──────────────────────────────────────────────

def plot_pca_variance(
    adata: ad.AnnData,
    n_pcs: int = 30,
    save_path: Optional[str] = None,
    figsize: Tuple = (7, 4),
) -> None:
    """Scree plot: explained variance per PC."""
    var_ratio = adata.uns["pca"]["variance_ratio"][:n_pcs]
    cumulative = np.cumsum(var_ratio)

    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()

    ax1.bar(range(1, n_pcs + 1), var_ratio, color="steelblue", alpha=0.7, label="Per-PC")
    ax2.plot(range(1, n_pcs + 1), cumulative, color="tomato", marker="o",
             markersize=3, linewidth=1.5, label="Cumulative")

    ax1.set_xlabel("Principal Component")
    ax1.set_ylabel("Explained Variance Ratio", color="steelblue")
    ax2.set_ylabel("Cumulative Variance", color="tomato")
    ax1.set_title("PCA — Explained Variance")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")

    plt.tight_layout()
    if save_path:
        save_fig(fig, save_path)
    else:
        plt.show()


# ──────────────────────────────────────────────
# Embedding Plots
# ──────────────────────────────────────────────

def plot_umap(
    adata: ad.AnnData,
    color: List[str],
    save_path: Optional[str] = None,
    title: str = "UMAP",
    **kwargs,
) -> None:
    """UMAP colored by one or more features/obs columns."""
    sc.pl.umap(adata, color=color, show=False, title=title, **kwargs)
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"[plot] Saved → {save_path}")


def plot_tsne(
    adata: ad.AnnData,
    color: List[str],
    save_path: Optional[str] = None,
    **kwargs,
) -> None:
    """t-SNE colored by one or more features/obs columns."""
    sc.pl.tsne(adata, color=color, show=False, **kwargs)
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"[plot] Saved → {save_path}")


# ──────────────────────────────────────────────
# Clustering / Marker Gene Plots
# ──────────────────────────────────────────────

def plot_marker_dotplot(
    adata: ad.AnnData,
    marker_genes: dict,
    groupby: str = "leiden",
    save_path: Optional[str] = None,
) -> None:
    """
    Dot plot of curated marker genes per cluster.

    Parameters
    ----------
    marker_genes : dict, e.g. {"T cells": ["CD3D", "CD3E"], "B cells": ["MS4A1"]}
    """
    sc.pl.dotplot(adata, marker_genes, groupby=groupby, show=False)
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()


def plot_rank_genes_heatmap(
    adata: ad.AnnData,
    n_genes: int = 10,
    groupby: str = "leiden",
    key: str = "rank_genes_groups",
    save_path: Optional[str] = None,
) -> None:
    """Heatmap of top marker genes per cluster."""
    sc.pl.rank_genes_groups_heatmap(
        adata, n_genes=n_genes, groupby=groupby, key=key,
        show_gene_labels=True, show=False
    )
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()


def plot_rank_genes_dotplot(
    adata: ad.AnnData,
    n_genes: int = 5,
    key: str = "rank_genes_groups",
    save_path: Optional[str] = None,
) -> None:
    """Dot plot of ranked marker genes."""
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, key=key, show=False)
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
