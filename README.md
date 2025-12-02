# Single_RNA_Sequencing 
  using Scanpy
Imagine you have a basket of fruit (your tissue sample) containing strawberries, bananas, and grapes.

    Bulk RNA Sequencing is like putting all that fruit into a blender.

        The Result: A smoothie.

        The Problem: You can taste that it's sweet, but you can't tell if the sweetness came from 50 ripe strawberries or just one really sugary banana. You lose the identity of the individual ingredients.

    Single-Cell RNA Sequencing (scRNA-seq) is like making a fruit salad.

        The Result: A bowl of separate fruit pieces.

        The Benefit: You can pick up each piece and say, "This is a grape," "This is a strawberry," and "Hey, this specific banana slice is rotten!"

Why do we need it?

If you are studying a brain tumor (the fruit basket), Bulk Sequencing might tell you "There are cancer genes here."

But Single-Cell Sequencing tells you:

    Who is who: "These are immune cells, these are healthy neurons, and these specific guys are the tumor cells."

    The Rare Suspect: It finds the "rotten banana" hiding in the corner—a rare cell type that might be driving the disease but would be invisible in a smoothie.


    To explain this using **Scanpy** (Single-Cell Analysis in Python), we need to understand that Scanpy uses a specific data structure called **AnnData** (Annotated Data). Think of AnnData as a smart Excel sheet that holds your gene counts, cell information (metadata), and gene information all in one object.

Here is how we perform the "Fruit Salad" analysis (scRNA-seq) using Scanpy commands, step-by-step:

### **1. Loading and QC (The Inspection)**
**Goal:** Load the data and throw away the "rotten fruit" (low-quality cells).

* **The Command:** `sc.read_10x_mtx()`
    * This loads your data into the AnnData object (often called `adata`).
* **The QC Check:** `sc.pp.calculate_qc_metrics(adata)`
    * Scanpy counts how many genes each cell has and calculates the percentage of mitochondrial reads.
* **The Filter:**
    * `sc.pp.filter_cells(adata, min_genes=200)`: "Throw away cells with fewer than 200 genes (empty droplets)."
    * `adata = adata[adata.obs['pct_counts_mt'] < 20]` : "Throw away cells with >20% mitochondrial reads (dying cells)."

### **2. Preprocessing (The Wash & Cut)**
**Goal:** Normalize the data so big cells don't overshadow small cells.

* **Normalization:** `sc.pp.normalize_total(adata, target_sum=1e4)`
    * This scales every cell to have 10,000 counts. It’s like weighing every fruit slice to make sure they are equal size before comparing.
* **Log Transformation:** `sc.pp.log1p(adata)`
    * This applies a logarithm to the data ($\log(x+1)$). It turns massive numbers (like 10,000 counts of a housekeeping gene) into manageable numbers so they don't skew the analysis.

### **3. Feature Selection (Picking the Distinctive Flavors)**
**Goal:** Ignore "boring" genes that are the same in every cell (housekeeping genes) and focus on genes that make cells different.

* **The Command:** `sc.pp.highly_variable_genes(adata)`
    * Scanpy identifies genes that vary significantly between cells. Instead of analyzing all 20,000 genes, we focus on the top ~2,000 "Highly Variable Genes" (HVGs) because they contain the real biological signal.

### **4. Dimensionality Reduction (The Map Sketch)**
**Goal:** Simplify the complex data into a basic map structure.

* **PCA (Principal Component Analysis):** `sc.tl.pca(adata)`
    * This rotates the data to find the main axes of variation. It compresses the 2,000 HVGs into roughly 30-50 "Principal Components."
* **Neighbors:** `sc.pp.neighbors(adata)`
    * This calculates which cells are similar to each other based on the PCA results. It builds a "graph" connecting every cell to its nearest neighbors.

### **5. Integration (The Blender - Optional but often needed)**
**Goal:** Remove batch effects (technical differences).

* *Note:* Standard Scanpy uses tools like **Ingest** or **BBKNN** for this. In the video you watched, they used **scvi-tools** (a separate Deep Learning package) at this stage because it is more powerful for complex data.
* If using pure Scanpy: `sc.external.pp.bbknn(adata)` or `sc.external.pp.harmony_integrate(adata)`.

### **6. Clustering and Visualization (The Plating)**
**Goal:** Group the cells and draw the final picture.

* **UMAP:** `sc.tl.umap(adata)`
    * This takes the "Neighbors" graph and squashes it down to 2D coordinates so we can plot it.
* **Clustering:** `sc.tl.leiden(adata)`
    * This uses the Leiden algorithm to detect "communities" in the neighbor graph. It assigns every cell a cluster number (0, 1, 2, etc.).
* **Plotting:** `sc.pl.umap(adata, color=['leiden'])`
    * This generates the colorful scatter plot you see in papers.

### **7. Annotation (The Taste Test)**
**Goal:** Identify what "Cluster 0" actually is.

* **Finding Markers:** `sc.tl.rank_genes_groups(adata, 'leiden')`
    * Scanpy asks: "What genes are unique to Cluster 0 compared to all other clusters?"
* **Checking Genes:** `sc.pl.dotplot(adata, var_names, groupby='leiden')`
    * You visualize these marker genes. If Cluster 0 has high expression of *CD3D*, you rename it "T-Cells."



This diagram summarizes the flow: from the raw matrix $\to$ AnnData object $\to$ Preprocessing $\to$ PCA/UMAP $\to$ Final Clusters.
