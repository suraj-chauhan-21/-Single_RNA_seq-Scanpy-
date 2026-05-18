# Data Directory

Place your raw input data here. The pipeline supports three formats:

## Option A — 10x Genomics CellRanger Output (recommended)

Place the CellRanger output folder here:
```
data/10x_output/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```
Set `INPUT_FORMAT = "10x"` and `DATA_10X_PATH` in `config.py`.

## Option B — H5AD (AnnData file)

```
data/raw_counts.h5ad
```
Set `INPUT_FORMAT = "h5ad"` and `H5AD_PATH` in `config.py`.

## Option C — CSV/TSV Count Matrix

```
data/counts.csv     # rows = cells, columns = genes
```
Set `INPUT_FORMAT = "csv"` and `CSV_PATH` in `config.py`.

---

## Public Datasets for Testing

| Dataset | Cells | Description | Source |
|---------|-------|-------------|--------|
| PBMC 3k | 2,700 | Peripheral blood mononuclear cells | [10x Genomics](https://www.10xgenomics.com/resources/datasets) |
| PBMC 68k | 68,000 | Large PBMC dataset | [10x Genomics](https://www.10xgenomics.com/resources/datasets) |
| Mouse brain | 1.3M | Large mouse brain atlas | [10x Genomics](https://www.10xgenomics.com/resources/datasets) |

### Quick Download (PBMC 3k):
```bash
wget https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
mv filtered_gene_bc_matrices/hg19 data/10x_output
```
