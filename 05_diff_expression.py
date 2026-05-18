"""
scripts/05_diff_expression.py
──────────────────────────────
Step 5: Find marker genes per cluster using Wilcoxon rank-sum (or other tests),
generate summary tables and visualizations.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import config as cfg
from utils import (
    load_adata, save_adata, ensure_dirs,
    plot_rank_genes_heatmap, plot_rank_genes_dotplot,
)


def run_diff_expression():
    sc.settings.verbosity = cfg.VERBOSITY
    sc.settings.set_figure_params(dpi=cfg.FIGURE["dpi"], facecolor="white")

    fig_dir = os.path.join(cfg.FIGURES_DIR, "05_diff_expression")
    ensure_dirs(fig_dir)

    # ── Load ──────────────────────────────────────────────────
    print("\n=== Step 5: Differential Expression & Marker Genes ===\n")
    adata = load_adata(cfg.ADATA_CLUSTERED_PATH)

    de = cfg.DIFF_EXPR
    groupby   = de["groupby"]
    method    = de["method"]
    n_genes   = de["n_genes"]
    key_added = de["key_added"]

    # ── Rank Genes per Group ───────────────────────────────────
    print(f"[de] Ranking marker genes per '{groupby}' cluster using '{method}'...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        use_raw=de["use_raw"],
        key_added=key_added,
        n_genes=n_genes,
    )

    # ── Visualizations ────────────────────────────────────────
    # 1. Dot plot of top N markers per cluster
    plot_rank_genes_dotplot(
        adata,
        n_genes=5,
        key=key_added,
        save_path=os.path.join(fig_dir, "dotplot_top5_markers.png"),
    )

    # 2. Heatmap
    plot_rank_genes_heatmap(
        adata,
        n_genes=5,
        groupby=groupby,
        key=key_added,
        save_path=os.path.join(fig_dir, "heatmap_top5_markers.png"),
    )

    # 3. Violin plots of top 3 markers per cluster
    groups = adata.obs[groupby].cat.categories.tolist()
    for group in groups[:min(len(groups), 6)]:   # limit to first 6 clusters
        try:
            top_genes = sc.get.rank_genes_groups_df(
                adata, group=str(group), key=key_added
            ).head(3)["names"].tolist()
            genes_present = [g for g in top_genes if g in adata.var_names]
            if genes_present:
                sc.pl.violin(adata, genes_present, groupby=groupby, show=False)
                plt.suptitle(f"Top Markers — Cluster {group}", y=1.02)
                plt.savefig(
                    os.path.join(fig_dir, f"violin_cluster{group}_top_markers.png"),
                    dpi=150, bbox_inches="tight"
                )
                plt.close()
        except Exception as e:
            print(f"[de] Skipping violin for cluster {group}: {e}")

    # 4. UMAP colored by top marker of first 3 clusters
    top_marker_genes = []
    for group in groups[:3]:
        try:
            df = sc.get.rank_genes_groups_df(adata, group=str(group), key=key_added)
            top_gene = df.iloc[0]["names"]
            if top_gene in adata.var_names:
                top_marker_genes.append(top_gene)
        except Exception:
            pass
    if top_marker_genes:
        sc.pl.umap(adata, color=top_marker_genes, show=False, use_raw=False, ncols=3)
        plt.savefig(
            os.path.join(fig_dir, "umap_top_markers_first3clusters.png"),
            dpi=150, bbox_inches="tight"
        )
        plt.close()

    # ── Export Marker Gene Tables ─────────────────────────────
    print("[de] Exporting marker gene tables...")
    results_list = []

    for group in groups:
        try:
            df = sc.get.rank_genes_groups_df(
                adata, group=str(group), key=key_added
            )
            df.insert(0, "cluster", str(group))
            results_list.append(df)
        except Exception as e:
            print(f"[de] Could not export group {group}: {e}")

    if results_list:
        all_markers = pd.concat(results_list, ignore_index=True)

        # Save full table
        all_markers_path = os.path.join(cfg.RESULTS_DIR, "all_cluster_markers.csv")
        all_markers.to_csv(all_markers_path, index=False)
        print(f"[de] All markers → {all_markers_path}")

        # Save top-N per cluster
        top_n = (
            all_markers.groupby("cluster")
            .head(n_genes)
            .reset_index(drop=True)
        )
        top_n_path = os.path.join(cfg.RESULTS_DIR, f"top{n_genes}_markers_per_cluster.csv")
        top_n.to_csv(top_n_path, index=False)
        print(f"[de] Top-{n_genes} markers → {top_n_path}")

        # Print top 5 per cluster to console
        print(f"\n[de] Top 5 marker genes per cluster:")
        for group in groups:
            grp_df = all_markers[all_markers["cluster"] == str(group)].head(5)
            genes_str = ", ".join(grp_df["names"].tolist())
            print(f"     Cluster {group:>3s}: {genes_str}")

    # ── Save Final AnnData ────────────────────────────────────
    save_adata(adata, cfg.ADATA_FINAL_PATH)
    print("\n=== Step 5 Complete ===")
    print(f"\n✓ Final AnnData saved:  {cfg.ADATA_FINAL_PATH}")
    print(f"✓ Figures in:           {cfg.FIGURES_DIR}")
    print(f"✓ Results in:           {cfg.RESULTS_DIR}\n")


if __name__ == "__main__":
    run_diff_expression()
