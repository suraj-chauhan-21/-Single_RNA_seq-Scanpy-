"""
utils/io_utils.py — Input/Output helper functions for scRNA-seq pipeline.
"""

import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np


def load_data(config: dict) -> ad.AnnData:
    """
    Load raw count data into an AnnData object.

    Supports:
      - 10x Genomics CellRanger output directory
      - H5AD (AnnData) files
      - CSV/TSV count matrices

    Parameters
    ----------
    config : dict
        Must contain keys: INPUT_FORMAT, DATA_10X_PATH, H5AD_PATH

    Returns
    -------
    adata : AnnData
    """
    fmt = config["INPUT_FORMAT"]

    if fmt == "10x":
        path = config["DATA_10X_PATH"]
        print(f"[io] Loading 10x Genomics data from: {path}")
        adata = sc.read_10x_mtx(
            path,
            var_names="gene_symbols",
            cache=True,
        )
        adata.var_names_make_unique()

    elif fmt == "h5ad":
        path = config["H5AD_PATH"]
        print(f"[io] Loading H5AD from: {path}")
        adata = sc.read_h5ad(path)

    elif fmt == "csv":
        path = config.get("CSV_PATH", "data/counts.csv")
        print(f"[io] Loading CSV from: {path}")
        df = pd.read_csv(path, index_col=0)
        adata = ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index),
                           var=pd.DataFrame(index=df.columns))

    else:
        raise ValueError(f"Unknown INPUT_FORMAT: {fmt}. Choose '10x', 'h5ad', or 'csv'.")

    print(f"[io] Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def save_adata(adata: ad.AnnData, path: str) -> None:
    """Save AnnData to H5AD file, creating parent directories if needed."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path)
    print(f"[io] Saved AnnData → {path}  ({adata.n_obs} cells × {adata.n_vars} genes)")


def load_adata(path: str) -> ad.AnnData:
    """Load AnnData from H5AD file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"AnnData file not found: {path}")
    adata = sc.read_h5ad(path)
    print(f"[io] Loaded AnnData ← {path}  ({adata.n_obs} cells × {adata.n_vars} genes)")
    return adata


def ensure_dirs(*paths: str) -> None:
    """Create directories if they don't exist."""
    for p in paths:
        os.makedirs(p, exist_ok=True)
