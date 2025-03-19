#=========================================
fileanme: pca
author: Lihuax
date: 2025-02-23 17:45
description: This file implements Principal Component Analysis (PCA) and UMAP dimensionality reduction methods for AnnData objects.
=========================================#

import ..Pp: normalize_total, normalize_total!
using Statistics, UMAP, TSVD

"""
    pca!(adata::Muon.AnnData; layer::String="log_transformed", n_pcs::Int=1000, verbose::Bool=true)

Perform Principal Component Analysis (PCA) on the specified data layer of an `AnnData` object.

If the specified `layer` does not exist in `adata.layers`, the function checks for a `"log_transformed"` layer.
If it is not found, it checks for a `"normalized"` layer to perform log-transformation.
If neither is available, the function will normalize and log-transform the original data `adata.X`.
After processing, the function computes the PCA using the first `n_pcs` columns (or the minimum available dimension if `n_pcs` exceeds the data dimensions)
and stores the resulting principal components in `adata.obsm["PCA"]`.

# Arguments
- `adata`: An instance of `Muon.AnnData` containing the data matrix and associated metadata.
- `layer`: The key of the data layer to be used for PCA (default: `"log_transformed"`).
- `n_pcs`: The number of principal components to compute (default: 1000).
- `verbose`: If `true`, prints informational messages (default: `true`).

# Returns
The updated `adata` with the PCA result stored in `adata.obsm["PCA"]`.

# Examples
```julia
adata = ...  # initialize AnnData
pca!(adata; layer="log_transformed", n_pcs=50, verbose=false)
```
"""
function pca!(
  adata::Muon.AnnData;
  layer::String="log_transformed",
  n_pcs::Int=1000,
  verbose::Bool=true,
)
  if !haskey(adata.layers, layer)
    @warn "layer $(layer) not found in `adata.layers`"
    if haskey(adata.layers, "log_transformed")
      verbose && @info "calculating PCA on log-transformed layer..."
    elseif haskey(adata.layers, "normalized")
      verbose && @info "log-transforming normalized counts before applying PCA..."
      log_transform!(adata, layer="normalized", verbose=verbose)
      verbose && @info "calculating PCA on log-transformed normalized counts..."
    else
      verbose && @info "normalizing and log-transforming before applying PCA..."
      normalize_total(adata)
      log_transform!(adata, layer="normalized", verbose=verbose)
      verbose && @info "calculating PCA on log-transformed normalized counts..."
    end
    X = adata.layers["log_transformed"]
  else
    X = adata.layers[layer]
    verbose && @info "calculating PCA on layer $(layer)..."
  end

  if n_pcs > minimum(size(X))
    @warn "not enough cells or variables available for calculating $(n_pcs) principal components, using $(minimum(size(X))) PCS."
    n_pcs = minimum(size(X))
  end

  pcs = prcomps(X[:, 1:n_pcs])

  adata.obsm["PCA"] = pcs[:, 1:n_pcs]
  return adata
end

"""
    umap!(adata::Muon.AnnData; layer::String="log_transformed", use_pca_init::Bool=false, n_pcs::Int=100, verbose::Bool=true, kwargs...)

Compute a UMAP embedding for the data in an `AnnData` object.

If `use_pca_init` is `true`, PCA is performed using `pca!(...)` to initialize the embedding.
If the specified `layer` does not exist in `adata.layers`, the function will check for a `"normalized"` layer and perform
normalization and log-transformation as needed.
The resulting data (either from PCA initialization or directly from the specified layer) is then passed to the UMAP algorithm.
The UMAP result, including the low-dimensional embedding, nearest neighbors (`knns`), distances (`knn_dists`),
and the fuzzy neighbor graph, is stored in the `adata` object.

# Arguments
- `adata`: An instance of `Muon.AnnData` containing the data matrix and associated metadata.
- `layer`: The key of the data layer to be used for UMAP (default: `"log_transformed"`).
- `use_pca_init`: If `true`, perform PCA initialization prior to UMAP (default: `false`).
- `n_pcs`: The number of principal components to compute if using PCA initialization (default: 100).
- `verbose`: If `true`, prints informational messages (default: `true`).
- `kwargs`: Additional keyword arguments passed to the UMAP algorithm.

# Returns
The updated `adata` with the UMAP embedding stored in `adata.obsm["umap"]`, along with nearest neighbor and graph information.

# Examples
```julia
adata = ...  # initialize AnnData
umap!(adata; layer="log_transformed", use_pca_init=true, n_pcs=50, verbose=false, n_neighbors=15)
```
"""
function umap!(
  adata::Muon.AnnData;
  layer::String="log_transformed",
  use_pca_init::Bool=false,
  n_pcs::Int=100,
  verbose::Bool=true,
  kwargs...,
)
  if use_pca_init
    pca!(adata; layer=layer, n_pcs=n_pcs, verbose=verbose)
    X = adata.obsm["PCA"]
  elseif !haskey(adata.layers, layer)
    @warn "layer $(layer) not found in `adata.layers`, calculating log-transformed normalized counts..."
    if !haskey(adata.layers, "normalized")
      normalize_total!(adata)
    end
    if !haskey(adata.layers, "log_transformed")
      log_transform!(adata, layer="normalized", verbose=verbose)
    end
    X = adata.layers["log_transformed"]
  else
    X = adata.layers[layer]
  end
  # calculate UMAP

  umap_result = UMAP.UMAP_(X'; kwargs...)

  # store results
  adata.obsm["umap"] = umap_result.embedding'
  adata.obsm["knns"] = umap_result.knns'
  adata.obsm["knn_dists"] = umap_result.dists'

  adata.obsp["fuzzy_neighbor_graph"] = umap_result.graph

  return adata
end

function log_transform!(adata::Muon.AnnData; layer::String="normalized", verbose::Bool=false)
  if !haskey(adata.layers, layer)
    @warn "layer $(layer) not found in `adata.layers`, defaulting to log + 1 transformation on `adata.X`..."
    logp1_transform!(adata; verbose=verbose)
    adata.layers["log_transformed"] = adata.layers["logp1_transformed"]
    return adata
  else
    verbose && @info "performing log transformation on layer $(layer)..."
    X = adata.layers[layer]
    adata.layers["log_transformed"] = log.(X .+ eps(eltype(X)))
    return adata
  end
end

function logp1_transform!(
  adata::Muon.AnnData;
  layer::Union{String, Nothing}=nothing,
  verbose::Bool=false,
)
  if haskey(adata.layers, layer)
    verbose && @info "performing log + 1 transformation on layer $(layer)..."
    X = adata.layers[layer]
  else
    verbose && @info "performing log +1 transformation on X..."
    X = adata.X
  end

  adata.layers["logp1_transformed"] = log.(X .+ one(eltype(X)))
  return adata
end

function standardize(x)
  (x .- Statistics.mean(x, dims=1)) ./ Statistics.std(x, dims=1)
end

function prcomps(mat::AbstractMatrix, standardizeinput=true)
  if standardizeinput
    mat = standardize(mat)
  end
  u, s, v = svd(mat)
  return u * Diagonal(s)
end

function prcomps(mat::SparseArrays.SparseMatrixCSC, standardizeinput=true, n_pcs=50)
  if standardizeinput
    mat = standardize(mat)
  end
  U, S, Vt = TSVD.tsvd(mat, n_pcs)
  return U * Diagonal(S)
end
