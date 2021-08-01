# Scorch: quick integration of multi-sample data with Harmony and a preprocessing Scanpy recipe

**Main use case**: multiple 10x RNA-seq samples with `.h5` inputs. 

Annotation of cell types working from a clustering solution can be difficult if the clustering solution is not ideal. This problem is amplified when using data integrated and clustered across samples and/or batches. To make this easier, I wrapped some key preprocessing functions from [***Scanpy***](https://scanpy.readthedocs.io/en/stable/) and functions from a [**pytorch implementation**](https://github.com/lilab-bcb/harmony-pytorch) of [***Harmony***](https://github.com/immunogenomics/harmony) [Korsunsky, et al., *Nat. Methods* **2019**](https://www.nature.com/articles/s41592-019-0619-0).

The Scorch class uses the recipe from Scanpy and pytorch-enabled Harmony implementation to quickly preprocess and integrate a directory of `.h5` 10x output matrices. 

***Note on runtime:*** ~22k cells across four samples took about 2.3 mins to preprocess, integrate, and visualize. The full notebook and outputs is [here](https://github.com/mvinyard/vintools/blob/main/notebooks/Scorch_example_4xBrain_10x.ipynb). This is actually sort of a bad example; each brain overlays quite well with minimal sample-to-sample difference. I'll eventually find a better native example. 

## Basic Usage:

This workflow has five main functions:
**1**. Read data
**2**. Harmony and UMAP embedding
**3**. Visualize Harmony-integrated results
**4**. Adjust clustering and re-visualize
**5**. Use marker genes to annotate cell types

### (1) Read data (.h5 files)

```python
from vintools.single_cell import Scorch

scorch = Scorch()

data_dir = "./neurons_h5_files/*.h5"
scorch.read_plural_10x_h5(data_dir)
```
>Harmony version: 0.1.6
Scanpy version: 1.8.1


The basic Scanpy preprocessing and data QC commands can be funneled into three general catagories. These functions have many arguments. The defaults are set but can be changed during execution of each function. 

```python
scorch.scanpy_qc()
scorch.scanpy_counts_analysis()
scorch.scanpy_dimensional_reduction_clustering()
```
At this point, check the output UMAP plot; samples will likely cluster alone. This is the "before" picture.

### (2) Harmony and UMAP embedding

Harmonize and create a UMAP embedding of the Harmony-generated dimensional reduction. 
```python
scorch.harmonize()
scorch.umap()
```

### (3) Visualize Harmony-integrated results
```python
scorch.visualize(plot_by="sample", separate_colors=True)

viz_genes = ["Ndnf", "Cldn5", "Ctss", "Gad1"]

for gene in viz_genes:
    scorch.visualize(plot_by=gene)
```

### (4) Adjust clustering and re-visualize
The default clustering solution may not be ideal for downstream cell type annotation. The real purpose of why I assembled this recipe / workflow is that I found it difficult to rapidly re-cluster, re-embed, and re-annotate cell types. This section and the following section will illustrate how easy this workflow makes this process.  To facilitate quick adjustment and subsequent visualization of various clustering solutions using the harmony-derived UMAP embedding, we can use the following functions to tweak the resolution: 

```
scorch.adjust_clustering_resolution(desired_resolution=0.2)
```

### (5) Use marker genes to annotate cell types

```python
hippocampus_marker_genes = [
    "Gad1",
    "Gad2",
     ...
    "Ctss",
    "C1qa",
]

scorch.cluster_dotplot(marker_genes=hippocampus_marker_genes)
```
