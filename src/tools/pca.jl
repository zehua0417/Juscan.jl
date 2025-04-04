#=========================================
fileanme: pca
author: Lihuax
date: 2025-02-23 17:45
description: This file implements Principal Component Analysis (PCA) and UMAP dimensionality reduction methods for AnnData objects.
=========================================#

import ..Pp: normalize_total, normalize_total!, _set_obs_rep!
using Statistics, UMAP, TSVD

"""
    pca!(adata::Muon.AnnData; layer="log_transformed", n_pcs=1000, key_added="pca", verbose=true)

Performs Principal Component Analysis (PCA) on the specified layer of an `AnnData` object and stores the result in `adata.obsm`.

If the specified layer is missing, the function will attempt to log-transform a normalized layer, or normalize and log-transform the raw counts if needed.

# Arguments
- `adata::Muon.AnnData`: The annotated data object on which to perform PCA.

## Keyword Arguments
- `layer::String = "log_transformed"`: The data layer to use for PCA. Defaults to `"log_transformed"`.
- `n_pcs::Int = 1000`: The number of principal components to compute. Automatically clipped to the smallest matrix dimension if too large.
- `key_added::String = "pca"`: The key under which to store the PCA result in `adata.obsm`.
- `verbose::Bool = true`: Whether to print progress messages.

# Returns
The modified `AnnData` object with PCA results stored in `adata.obsm[key_added]`.

# Notes
- This function performs automatic preprocessing if the requested layer is not present.
- The PCA is computed via SVD on standardized data.

"""
function pca!(
  adata::Muon.AnnData;
  layer::String="log_transformed",
  n_pcs::Int=1000,
  key_added::String="pca",
  verbose::Bool=true,
)
  if !haskey(adata.layers, layer)
    @warn "layer $(layer) not found in `adata.layers`"
    if haskey(adata.layers, "log_transformed")
      verbose && @info "calculating PCA on log-transformed layer..."
    elseif haskey(adata.layers, "normalized")
      verbose && @info "log-transforming normalized counts before applying PCA..."
      log_transform!(adata, layer="normalized", key_added="log_transformed", verbose=verbose)
      verbose && @info "calculating PCA on log-transformed normalized counts..."
    else
      verbose && @info "normalizing and log-transforming before applying PCA..."
      norm_res = normalize_total(adata)
      _set_obs_rep!(adata, norm_res["X"], layer="normalized")
      log_transform!(adata, layer="normalized", key_added="log_transformed", verbose=verbose)
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

  adata.obsm[key_added] = pcs[:, 1:n_pcs]
  return adata
end

"""
    umap!(adata::Muon.AnnData; layer="log_transformed", use_pca=nothing, n_pcs=100, key_added="umap", verbose=true, kwargs...)

Computes a UMAP embedding from the data in the specified layer or PCA representation and stores the result in `adata.obsm`.

If PCA is requested via `use_pca`, it will be computed automatically if not already present.

# Arguments
- `adata::Muon.AnnData`: The annotated data object on which to compute UMAP.

## Keyword Arguments
- `layer::String = "log_transformed"`: The layer to use as input if `use_pca` is not specified.
- `use_pca::Union{String, Nothing} = nothing`: If specified, use this key in `adata.obsm` as PCA input. If missing, it will be computed.
- `n_pcs::Int = 100`: Number of principal components to use if PCA needs to be computed.
- `key_added::String = "umap"`: The key under which to store the UMAP embedding.
- `verbose::Bool = true`: Whether to print progress messages.
- `kwargs...`: Additional keyword arguments passed to `UMAP.UMAP_()`.

# Returns
The modified `AnnData` object with UMAP results stored in:
- `adata.obsm[key_added]`: The UMAP embedding.
- `adata.obsm["knns"]`: K-nearest neighbors matrix.
- `adata.obsm["knn_dists"]`: KNN distance matrix.
- `adata.obsp["fuzzy_neighbor_graph"]`: Fuzzy graph representation.

# Notes
- Automatically performs normalization and log transformation if necessary.
- Uses the `UMAP.jl` package under the hood.

"""
function umap!(
  adata::Muon.AnnData;
  layer::String="log_transformed",
  use_pca::Union{AbstractString, Nothing}=nothing,
  n_pcs::Int=100,
  key_added::String="umap",
  verbose::Bool=true,
  kwargs...,
)
  if !isnothing(use_pca)
    pca!(adata; layer=layer, n_pcs=n_pcs, verbose=verbose)
    X = adata.obsm[use_pca]
  elseif !haskey(adata.layers, layer)
    @warn "layer $(layer) not found in `adata.layers`, calculating log-transformed normalized counts..."
    if !haskey(adata.layers, "normalized")
      normalize_total!(adata)
    end
    if !haskey(adata.layers, "log_transformed")
      log_transform!(adata, layer="normalized", key_added="log_transformed", verbose=verbose)
    end
    X = adata.layers["log_transformed"]
  else
    X = adata.layers[layer]
  end
  # calculate UMAP

  umap_result = UMAP.UMAP_(X'; kwargs...)

  # store results
  adata.obsm[key_added] = umap_result.embedding'
  adata.obsm["knns"] = umap_result.knns'
  adata.obsm["knn_dists"] = umap_result.dists'

  adata.obsp["fuzzy_neighbor_graph"] = umap_result.graph

  return adata
end

"""
    log_transform!(adata::Muon.AnnData; layer="normalized", key_added="log_transformed", verbose=false)

Applies a log transformation to the specified data layer of an `AnnData` object and stores the result in `adata.layers`.

If the specified layer is missing, the function defaults to applying a log(1 + x) transformation to `adata.X`.

# Arguments
- `adata::Muon.AnnData`: The data object to transform.

## Keyword Arguments
- `layer::String = "normalized"`: The layer to transform. Must exist in `adata.layers`.
- `key_added::String = "log_transformed"`: The key to store the result under in `adata.layers`.
- `verbose::Bool = false`: Whether to print messages during the process.

# Returns
The modified `AnnData` object with the log-transformed data added to `adata.layers[key_added]`.

# Notes
- The transformation is log(x + ϵ), where ϵ is a small constant to avoid log(0).
- For default fallback behavior, see `logp1_transform!()`.

"""
function log_transform!(
  adata::Muon.AnnData;
  layer::String="normalized",
  key_added::AbstractString="log_transformed",
  verbose::Bool=false,
)
  if !haskey(adata.layers, layer)
    @warn "layer $(layer) not found in `adata.layers`, defaulting to log + 1 transformation on `adata.X`..."
    logp1_transform!(adata; key_added="logp1_transformed", verbose=verbose)
    adata.layers[key_added] = adata.layers["logp1_transformed"]
    return adata
  else
    verbose && @info "performing log transformation on layer $(layer)..."
    X = adata.layers[layer]
    adata.layers[key_added] = log.(X .+ eps(eltype(X)))
    return adata
  end
end

"""
    logp1_transform!(adata::Muon.AnnData; layer=nothing, key_added="log1_transformed", verbose=false)

Applies a log(1 + x) transformation to the specified layer or the main data matrix `adata.X` in an `AnnData` object. The result is stored in `adata.layers[key_added]`.

# Arguments
- `adata::Muon.AnnData`: The annotated data object to transform.

## Keyword Arguments
- `layer::Union{String, Nothing} = nothing`: The name of the data layer to transform. If `nothing` or the layer is missing, uses `adata.X`.
- `key_added::AbstractString = "log1_transformed"`: The name under which to store the transformed result in `adata.layers`.
- `verbose::Bool = false`: Whether to print transformation messages.

# Returns
The modified `AnnData` object with the log(1 + x) transformed data stored in `adata.layers[key_added]`.

# Notes
- This transformation is commonly used to stabilize variance and reduce the effect of outliers.
- Compared to `log_transform!`, this version adds 1 to the data before taking the logarithm, making it more robust to zero entries.

"""
function logp1_transform!(
  adata::Muon.AnnData;
  layer::Union{String, Nothing}=nothing,
  key_added::AbstractString="log1_transformed",
  verbose::Bool=false,
)
  if haskey(adata.layers, layer)
    verbose && @info "performing log + 1 transformation on layer $(layer)..."
    X = adata.layers[layer]
  else
    verbose && @info "performing log +1 transformation on X..."
    X = adata.X
  end

  adata.layers[key_added] = log.(X .+ one(eltype(X)))
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
