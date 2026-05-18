"""
utils — Utility modules for scRNA-seq Scanpy pipeline.
"""
from .io_utils import load_data, save_adata, load_adata, ensure_dirs
from .qc_utils import compute_qc_metrics, filter_cells_and_genes, detect_doublets
from .plot_utils import (
    plot_qc_metrics, plot_qc_scatter, plot_highest_expressed_genes,
    plot_hvg, plot_pca_variance, plot_umap, plot_tsne,
    plot_marker_dotplot, plot_rank_genes_heatmap, plot_rank_genes_dotplot,
)

__all__ = [
    "load_data", "save_adata", "load_adata", "ensure_dirs",
    "compute_qc_metrics", "filter_cells_and_genes", "detect_doublets",
    "plot_qc_metrics", "plot_qc_scatter", "plot_highest_expressed_genes",
    "plot_hvg", "plot_pca_variance", "plot_umap", "plot_tsne",
    "plot_marker_dotplot", "plot_rank_genes_heatmap", "plot_rank_genes_dotplot",
]
