# Plotting Module

The `plots` module provides a set of visualization tools designed to make single-cell data exploration intuitive and informative. These functions enable users to create publication-ready figures for quality control, feature exploration, dimensionality reduction, and clustering results.

## Violin and Scatter Plots

Violin and scatter plots are essential for examining cell-level metrics, such as gene counts, total expression, or mitochondrial content. The `violin` function offers compact summaries of distributions, while `scatter` allows for flexible 2D comparisons between features, optionally colored by additional metadata.

## Highly Variable Genes Visualization

The `hvg_scatter` function provides a dual-panel visualization of gene variability metrics, distinguishing highly variable genes (HVGs) from background genes. This plot is especially useful for evaluating gene selection prior to dimensionality reduction.

## Dimensionality Reduction Plots

The `plot_variance_ratio` function visualizes the explained variance of each principal component, helping users decide how many PCs to retain. `plot_umap` projects cells into a low-dimensional embedding using UMAP, colored by user-specified labels to reveal structure and cluster separation.

## Color Palettes

Color plays a vital role in visual clarity. Juscan.jl includes palette expansion and colormap utilities to generate well-balanced, customizable color schemes. These utilities ensure consistency across all plots.

By combining visual elegance with analytical depth, the `plots` module empowers users to communicate insights clearly and effectively.

---

```@index
Pages = ["plots.md"]
```

## QC and Feature Plots

```@docs
violin
scatter
```

## HVG Visualization

```@docs
hvg_scatter
```

## Dimensionality Reduction

```@docs
plot_variance_ratio
plot_umap
```

## Color Utilities

```@docs
expand_palette
get_continuous_colormap
```
