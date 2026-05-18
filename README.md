# Single-Cell RNA-seq Analysis with Scanpy

A complete pipeline for single-cell RNA sequencing (scRNA-seq) data analysis using [Scanpy](https://scanpy.readthedocs.io/), an open-source Python toolkit for analyzing single-cell gene expression data.

---

## 📁 Project Structure

```
scrna_scanpy/
├── README.md
├── requirements.txt
├── environment.yml
├── config.py
├── data/
│   └── README_data.md
├── notebooks/
│   ├── 01_preprocessing_qc.ipynb
│   ├── 02_normalization_feature_selection.ipynb
│   ├── 03_dimensionality_reduction.ipynb
│   ├── 04_clustering.ipynb
│   └── 05_differential_expression.ipynb
├── scripts/
│   ├── 01_preprocessing.py
│   ├── 02_normalization.py
│   ├── 03_dim_reduction.py
│   ├── 04_clustering.py
│   ├── 05_diff_expression.py
│   └── run_pipeline.py
├── utils/
│   ├── __init__.py
│   ├── qc_utils.py
│   ├── plot_utils.py
│   └── io_utils.py
└── results/
    └── README_results.md
```

---

## 🔬 Analysis Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_preprocessing.py` | Quality control, filtering low-quality cells/genes |
| 2 | `02_normalization.py` | Normalization, log-transformation, HVG selection |
| 3 | `03_dim_reduction.py` | PCA, UMAP, t-SNE |
| 4 | `04_clustering.py` | Leiden/Louvain clustering, cluster annotation |
| 5 | `05_diff_expression.py` | Marker genes & differential expression |

---

## 🚀 Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/suraj-chauhan-21/-Single_RNA_seq-Scanpy-
cd -Single_RNA_seq-Scanpy-
```

### 2. Set up environment
```bash
conda env create -f environment.yml
conda activate scrna_scanpy
```

Or with pip:
```bash
pip install -r requirements.txt
```

### 3. Configure paths
Edit `config.py` to set your data paths and analysis parameters.

### 4. Run the full pipeline
```bash
python scripts/run_pipeline.py
```

Or run individual steps:
```bash
python scripts/01_preprocessing.py
python scripts/02_normalization.py
python scripts/03_dim_reduction.py
python scripts/04_clustering.py
python scripts/05_diff_expression.py
```

---

## 📊 Input Data

The pipeline accepts:
- **10x Genomics** output (matrix.mtx, barcodes.tsv, features.tsv)
- **H5AD** files (AnnData format)
- **CSV/TSV** count matrices

Place your data in the `data/` directory and update `config.py` accordingly.

---

## 📦 Dependencies

- Python ≥ 3.8
- scanpy ≥ 1.9
- anndata
- pandas, numpy, matplotlib, seaborn
- leidenalg (for clustering)

See `requirements.txt` or `environment.yml` for full dependency list.

---

## 📝 Citation

If you use Scanpy in your research, please cite:

> Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0
