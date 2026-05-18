"""
utils/qc_utils.py — Quality control helper functions.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional


def compute_qc_metrics(adata: ad.AnnData, mito_prefix: str = "MT-") -> ad.AnnData:
    """
    Compute QC metrics and add them to adata.obs.

    Adds:
        - n_genes_by_counts  : number of genes with >0 counts per cell
        - total_counts       : total UMI counts per cell
        - pct_counts_mt      : % mitochondrial counts

    Parameters
    ----------
    adata       : AnnData
    mito_prefix : str, prefix for mitochondrial genes (human: "MT-", mouse: "mt-")

    Returns
    -------
    adata with QC metrics added to .obs and .var
    """
    # Flag mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)

    # Flag ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # Compute QC metrics via Scanpy
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    print("[qc] QC metrics computed.")
    print(f"     Median genes/cell  : {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"     Median counts/cell : {adata.obs['total_counts'].median():.0f}")
    print(f"     Median %mito       : {adata.obs['pct_counts_mt'].median():.2f}%")
    return adata


def filter_cells_and_genes(adata: ad.AnnData, qc_params: dict) -> ad.AnnData:
    """
    Filter cells and genes based on QC thresholds.

    Parameters
    ----------
    adata      : AnnData with QC metrics already computed
    qc_params  : dict with keys from config.QC

    Returns
    -------
    Filtered AnnData
    """
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=qc_params["min_cells"])

    # Filter cells by gene count
    sc.pp.filter_cells(adata, min_genes=qc_params["min_genes"])
    adata = adata[adata.obs["n_genes_by_counts"] < qc_params["max_genes"], :]

    # Filter cells by total counts
    adata = adata[adata.obs["total_counts"] >= qc_params["min_counts"], :]
    adata = adata[adata.obs["total_counts"] <= qc_params["max_counts"], :]

    # Filter cells by mitochondrial %
    adata = adata[adata.obs["pct_counts_mt"] < qc_params["max_pct_mito"], :]

    print(f"[qc] Cells  : {n_cells_before} → {adata.n_obs}  "
          f"(removed {n_cells_before - adata.n_obs})")
    print(f"[qc] Genes  : {n_genes_before} → {adata.n_vars}  "
          f"(removed {n_genes_before - adata.n_vars})")
    return adata


def detect_doublets(adata: ad.AnnData, expected_rate: float = 0.06) -> ad.AnnData:
    """
    Optionally run Scrublet for doublet detection (if installed).
    Adds 'doublet_score' and 'predicted_doublet' to adata.obs.

    Parameters
    ----------
    adata         : AnnData
    expected_rate : float, expected doublet rate (default 6%)
    """
    try:
        import scrublet as scr
        print("[qc] Running Scrublet doublet detection...")
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_rate)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets
        n_doublets = predicted_doublets.sum()
        print(f"[qc] Detected {n_doublets} predicted doublets "
              f"({n_doublets / adata.n_obs * 100:.1f}%)")
    except ImportError:
        print("[qc] Scrublet not installed — skipping doublet detection. "
              "Install with: pip install scrublet")
    return adata
