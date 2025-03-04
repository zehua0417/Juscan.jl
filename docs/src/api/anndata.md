# Anndata utils

The `AnnData` struct is imported from [`Muon.jl`](https://github.com/scverse/Muon.jl). The package provides read and write functions for `.h5ad` and `.h5mu` files, the typical H5-based format for storing Python `anndata` objects. The `AnnData` object stores datasets together with metadata, such as information on the variables (genes in scRNA-seq data) and observations (cells), as well as different kinds of annotations and transformations of the original count matrix, such as PCA or UMAP embeddings, or graphs of observations or variables.

For details on the Julia implementation in `Muon.jl`, see the [documentation](https://scverse.github.io/Muon.jl/dev/).

For more details on the original Python implementation of the `anndata` object, see the [documentation](https://anndata.readthedocs.io/en/latest/) and [preprint](https://doi.org/10.1101/2021.12.16.473007).

![anndata](../assets/anndata_schema.svg)

```@index
Pages = ["anndata.md"]
```

## celltypes

```@docs
get_celltypes
```

## subset adata

```@docs
subset_adata
subset_adata!
```

## insert

```@docs
insert_obs!
insert_var!
```
