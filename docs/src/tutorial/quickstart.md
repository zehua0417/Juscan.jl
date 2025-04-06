# Quick Start Guide for Juscan.jl

This guide walks you through a basic single-cell RNA-seq analysis pipeline using [`Juscan.jl`](https://github.com/zehua0417/Juscan.jl). The dataset used in this example comes from [doi:10.6084/m9.figshare.22716739.v1](https://doi.org/10.6084/m9.figshare.22716739.v1).

## 📦 Prerequisites

Make sure you have installed `Juscan.jl`, `Muon.jl`, and their dependencies:

```julia
using Pkg
Pkg.activate(".")
Pkg.add(url="https://github.com/zehua0417/Juscan.jl")
```

## 📁 Load Data

```julia
using Juscan
using Muon
using DataFrames

adata = Juscan.readh5ad("data/data.h5ad")
```

## 🔍 Quality Control

```julia
Juscan.Pp.filter_cells!(adata, min_genes=100)
Juscan.Pp.filter_genes!(adata, min_cells=3)
Juscan.Pp.filter_cells!(adata, max_genes=8000)
Juscan.Pp.filter_cells!(adata, max_counts=140000)

adata.var.mt = startswith.(adata.var_names, "MT-")
adata.var.ribo = startswith.(adata.var_names, "RPS") .| startswith.(adata.var_names, "RPL")
adata.var.hb = occursin.(r"^HB[^P]", adata.var_names)

Juscan.Pp.calculate_qc_metrics!(adata, qc_vars=["mt", "ribo", "hb"])
```

### 📊 Visualize QC Metrics

```julia
Juscan.Pl.violin(
  adata,
  ["pct_counts_mt", "n_genes_by_counts", "total_counts"];
  width=300,
  height=800,
  fill_alpha=0.7,
  savefig="/home/lihuax/Pictures/Juscan/qc_violin.png",
)

Juscan.Pl.scatter(
  adata,
  "total_counts",
  "n_genes_by_counts",
  color_key="pct_counts_mt",
  width=800,
  height=800,
  colormap_name="magma",
  savefig="/home/lihuax/Pictures/Juscan/qc_scatter.png",
)
```

## 🔬 Normalization

```julia
adata.layers["normalized"] = deepcopy(adata.X)
Juscan.Pp.normalize_total!(adata, target_sum=1000, layer="normalized")
Juscan.Tl.logp1_transform!(adata, layer="normalized", key_added="normalized_logp1")
adata.layers["normalized_logp1"] = Float64.(adata.layers["normalized_logp1"])
```

## 🧬 Highly Variable Genes

```julia
Juscan.Tl.highly_variable_genes!(adata, n_top_genes=2000, layer="normalized_logp1")
Juscan.Pl.hvg_scatter(adata, savefig="/home/lihuax/Pictures/Juscan/hvg_scatter.png")
```

## 📉 Dimensionality Reduction

```julia
Juscan.Tl.pca!(adata; key_added="pca", n_pcs=15)
Juscan.Tl.pca!(adata; layer="normalized_logp1", key_added="pca", n_pcs=15)
Juscan.Pl.plot_variance_ratio(adata, savefig="/home/lihuax/Pictures/Juscan/variance_ratio.png")

Juscan.Tl.subset_to_hvg!(adata; layer="normalized_logp1", n_top_genes=2000)
Juscan.Pp.filter_cells!(adata, min_counts=3000)
Juscan.Pp.filter_cells!(adata, max_counts=10000)
```

## 🔗 Clustering & UMAP

```julia
Juscan.Tl.clustering!(adata, method="km", use_pca=15, cluster_K=4, dist="Euclidean")
Juscan.Tl.umap!(
  adata;
  layer="normalized_logp1",
  key_added="umap",
  n_pcs=15,
  min_dist=0.5,
  n_neighbors=50,
)
fig = Juscan.Pl.plot_umap(adata, color_by="clusters_0.5")
```

---

This pipeline demonstrates the essential steps of scRNA-seq analysis using `Juscan.jl`. For more detailed usage, please refer to the [official documentation](https://zehua0417.github.io/Juscan.jl/).

Happy analyzing! 🧬✨
