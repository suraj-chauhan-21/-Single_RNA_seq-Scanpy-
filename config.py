"""
config.py — Central configuration for the scRNA-seq Scanpy pipeline.
Edit these parameters before running the pipeline.
"""

import os

# ─────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")

# Input: choose ONE of the formats below
# 10x Genomics (CellRanger output folder with matrix.mtx, barcodes.tsv, features.tsv)
DATA_10X_PATH = os.path.join(DATA_DIR, "10x_output")

# Alternatively: a pre-built H5AD file
H5AD_PATH = os.path.join(DATA_DIR, "raw_counts.h5ad")

# Saved AnnData objects between steps
ADATA_PROCESSED_PATH  = os.path.join(RESULTS_DIR, "adata_preprocessed.h5ad")
ADATA_NORMALIZED_PATH = os.path.join(RESULTS_DIR, "adata_normalized.h5ad")
ADATA_REDUCED_PATH    = os.path.join(RESULTS_DIR, "adata_reduced.h5ad")
ADATA_CLUSTERED_PATH  = os.path.join(RESULTS_DIR, "adata_clustered.h5ad")
ADATA_FINAL_PATH      = os.path.join(RESULTS_DIR, "adata_final.h5ad")

# ─────────────────────────────────────────────
# INPUT FORMAT
# ─────────────────────────────────────────────
INPUT_FORMAT = "10x"  # Options: "10x" | "h5ad" | "csv"

# ─────────────────────────────────────────────
# QUALITY CONTROL PARAMETERS
# ─────────────────────────────────────────────
QC = {
    "min_genes":          200,    # Min genes per cell
    "max_genes":          6000,   # Max genes per cell (filter doublets/multiplets)
    "min_cells":          3,      # Min cells a gene must appear in
    "max_pct_mito":       20,     # Max % mitochondrial reads per cell (%)
    "mito_prefix":        "MT-",  # Prefix for mitochondrial genes (human: "MT-", mouse: "mt-")
    "min_counts":         500,    # Minimum UMI counts per cell
    "max_counts":         50000,  # Maximum UMI counts per cell
}

# ─────────────────────────────────────────────
# NORMALIZATION PARAMETERS
# ─────────────────────────────────────────────
NORMALIZATION = {
    "target_sum":         1e4,    # Normalize each cell to this total count
    "log1p":              True,   # Apply log1p transformation after normalization
    "n_top_genes":        2000,   # Number of highly variable genes to select
    "hvg_flavor":         "seurat",  # "seurat" | "cell_ranger" | "seurat_v3"
    "hvg_min_mean":       0.0125,
    "hvg_max_mean":       3,
    "hvg_min_disp":       0.5,
}

# ─────────────────────────────────────────────
# DIMENSIONALITY REDUCTION PARAMETERS
# ─────────────────────────────────────────────
DIM_REDUCTION = {
    "n_comps_pca":        50,     # Number of PCA components
    "n_neighbors":        15,     # k-nearest neighbors
    "n_pcs":              40,     # PCs to use for neighborhood graph
    "umap_min_dist":      0.3,
    "umap_spread":        1.0,
    "tsne_perplexity":    30,
    "random_state":       42,
}

# ─────────────────────────────────────────────
# CLUSTERING PARAMETERS
# ─────────────────────────────────────────────
CLUSTERING = {
    "method":             "leiden",   # "leiden" | "louvain"
    "resolution":         0.5,        # Higher = more clusters
    "random_state":       42,
}

# ─────────────────────────────────────────────
# DIFFERENTIAL EXPRESSION PARAMETERS
# ─────────────────────────────────────────────
DIFF_EXPR = {
    "method":             "wilcoxon",  # "wilcoxon" | "t-test" | "logreg"
    "n_genes":            25,          # Top marker genes per cluster
    "groupby":            "leiden",    # Column in adata.obs to group by
    "use_raw":            True,
    "key_added":          "rank_genes_groups",
}

# ─────────────────────────────────────────────
# FIGURE SETTINGS
# ─────────────────────────────────────────────
FIGURE = {
    "dpi":                150,
    "format":             "png",   # "png" | "svg" | "pdf"
    "figsize":            (8, 6),
    "color_map":          "viridis",
}

# ─────────────────────────────────────────────
# MISC
# ─────────────────────────────────────────────
RANDOM_SEED   = 42
N_JOBS        = -1   # Parallel jobs (-1 = use all CPUs)
VERBOSITY     = 1    # Scanpy verbosity: 0=errors, 1=warnings, 2=info, 3=hints
