"""
scripts/04_clustering.py
─────────────────────────
Step 4: Leiden/Louvain clustering, cluster visualization,
and optional manual cell-type annotation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import config as cfg
from utils import (
    load_adata, save_adata, ensure_dirs,
    plot_umap, plot_tsne, plot_marker_dotplot,
)


# ─────────────────────────────────────────────────────────
# Optional: Edit this dict to annotate clusters after
# inspecting marker genes from step 5.
# Keys = cluster IDs (strings), Values = cell-type labels.
# ─────────────────────────────────────────────────────────
CLUSTER_ANNOTATIONS = {
    # "0": "CD4+ T Cells",
    # "1": "CD8+ T Cells",
    # "2": "B Cells",
    # "3": "NK Cells",
    # "4": "Monocytes",
    # "5": "Dendritic Cells",
}

# Well-known marker genes for PBMC cell types (edit for your tissue)
KNOWN_MARKERS = {
    "T Cells":        ["CD3D", "CD3E", "CD3G"],
    "CD4+ T Cells":   ["CD4", "IL7R", "CCR7"],
    "CD8+ T Cells":   ["CD8A", "CD8B"],
    "B Cells":        ["MS4A1", "CD79A", "CD19"],
    "NK Cells":       ["GNLY", "NKG7", "KLRB1"],
    "Monocytes":      ["CD14", "LYZ", "CST3"],
    "DC":             ["FCER1A", "CST3"],
    "Megakaryocytes": ["PPBP", "PF4"],
}


def run_clustering():
    sc.settings.verbosity = cfg.VERBOSITY
    sc.settings.set_figure_params(dpi=cfg.FIGURE["dpi"], facecolor="white")

    fig_dir = os.path.join(cfg.FIGURES_DIR, "04_clustering")
    ensure_dirs(fig_dir)

    # ── Load ──────────────────────────────────────────────────
    print("\n=== Step 4: Clustering ===\n")
    adata = load_adata(cfg.ADATA_REDUCED_PATH)

    cl = cfg.CLUSTERING
    method = cl["method"]

    # ── Clustering ────────────────────────────────────────────
    print(f"[cl] Running {method} clustering (resolution={cl['resolution']})...")
    if method == "leiden":
        sc.tl.leiden(
            adata,
            resolution=cl["resolution"],
            random_state=cl["random_state"],
            key_added="leiden",
        )
        cluster_key = "leiden"
    elif method == "louvain":
        sc.tl.louvain(
            adata,
            resolution=cl["resolution"],
            random_state=cl["random_state"],
            key_added="louvain",
        )
        cluster_key = "louvain"
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    n_clusters = adata.obs[cluster_key].nunique()
    print(f"[cl] Found {n_clusters} clusters.")

    # Print cluster sizes
    cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
    print("\n[cl] Cluster sizes:")
    for cid, cnt in cluster_counts.items():
        print(f"     Cluster {cid:>3s}: {cnt:>5d} cells")

    # ── UMAP with clusters ────────────────────────────────────
    plot_umap(
        adata,
        color=[cluster_key],
        save_path=os.path.join(fig_dir, f"umap_{cluster_key}.png"),
        legend_loc="on data",
        title=f"UMAP — {method.capitalize()} Clusters",
    )

    # ── t-SNE with clusters ───────────────────────────────────
    plot_tsne(
        adata,
        color=[cluster_key],
        save_path=os.path.join(fig_dir, f"tsne_{cluster_key}.png"),
        legend_loc="on data",
    )

    # ── Known Marker Genes on UMAP ────────────────────────────
    # Plot only markers present in the dataset
    all_markers = [g for genes in KNOWN_MARKERS.values() for g in genes]
    present = [g for g in all_markers if g in adata.var_names]
    if present:
        sc.pl.umap(
            adata, color=present[:min(12, len(present))],
            ncols=4, show=False, use_raw=False
        )
        plt.savefig(os.path.join(fig_dir, "umap_known_markers.png"),
                    dpi=150, bbox_inches="tight")
        plt.close()
        print(f"[cl] Plotted {len(present)} known marker genes on UMAP.")

    # ── Dot Plot of Known Markers ─────────────────────────────
    # Filter dict to genes present in dataset
    markers_present = {
        ct: [g for g in genes if g in adata.var_names]
        for ct, genes in KNOWN_MARKERS.items()
    }
    markers_present = {ct: g for ct, g in markers_present.items() if g}
    if markers_present:
        plot_marker_dotplot(
            adata,
            marker_genes=markers_present,
            groupby=cluster_key,
            save_path=os.path.join(fig_dir, "dotplot_known_markers.png"),
        )

    # ── Optional: Annotate Clusters ───────────────────────────
    if CLUSTER_ANNOTATIONS:
        print("\n[cl] Applying cluster annotations...")
        adata.obs["cell_type"] = (
            adata.obs[cluster_key]
            .map(CLUSTER_ANNOTATIONS)
            .fillna("Unknown")
            .astype("category")
        )
        plot_umap(
            adata,
            color=["cell_type"],
            save_path=os.path.join(fig_dir, "umap_cell_type.png"),
            title="UMAP — Cell Types",
            legend_loc="on data",
        )

        # Save annotation table
        annot_df = pd.DataFrame({
            "cluster": list(CLUSTER_ANNOTATIONS.keys()),
            "cell_type": list(CLUSTER_ANNOTATIONS.values()),
        })
        annot_path = os.path.join(cfg.RESULTS_DIR, "cluster_annotations.csv")
        annot_df.to_csv(annot_path, index=False)
        print(f"[cl] Annotations saved → {annot_path}")

    # ── Save ──────────────────────────────────────────────────
    save_adata(adata, cfg.ADATA_CLUSTERED_PATH)
    print("\n=== Step 4 Complete ===\n")


if __name__ == "__main__":
    run_clustering()
