# Results Directory

This folder is auto-populated when the pipeline runs.

## Expected Output Files

| File | Created By | Description |
|------|-----------|-------------|
| `adata_preprocessed.h5ad` | Step 1 | QC-filtered AnnData |
| `adata_normalized.h5ad` | Step 2 | Normalized + HVG-subset AnnData |
| `adata_reduced.h5ad` | Step 3 | + PCA, kNN, UMAP, t-SNE |
| `adata_clustered.h5ad` | Step 4 | + Leiden clusters |
| `adata_final.h5ad` | Step 5 | + Differential expression results |
| `cluster_annotations.csv` | Step 4 | Cluster → cell type mapping |
| `all_cluster_markers.csv` | Step 5 | All ranked marker genes |
| `top25_markers_per_cluster.csv` | Step 5 | Top 25 markers per cluster |

## Figures

Figures are saved to `results/figures/` organised by step:

```
figures/
├── 01_qc/
│   ├── qc_violin_before.png
│   ├── qc_violin_after.png
│   ├── qc_scatter_before.png
│   ├── qc_scatter_after.png
│   └── highest_expressed_genes.png
├── 02_normalization/
│   └── hvg_plot.png
├── 03_dim_reduction/
│   ├── pca_variance.png
│   ├── pca_qc_coloring.png
│   ├── umap_qc.png
│   └── tsne_qc.png
├── 04_clustering/
│   ├── umap_leiden.png
│   ├── tsne_leiden.png
│   ├── umap_known_markers.png
│   ├── dotplot_known_markers.png
│   └── umap_cell_type.png         (if annotated)
└── 05_diff_expression/
    ├── dotplot_top5_markers.png
    ├── heatmap_top5_markers.png
    ├── violin_cluster*_top_markers.png
    └── umap_top_markers_first3clusters.png
```
