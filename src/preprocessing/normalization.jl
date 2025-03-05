#=========================================
fileanme: normalization
author: Lihuax
date: 2025-02-15 15:04
description:
=========================================#

using Statistics
import ..Juscan: _get_obs_rep, _set_obs_rep!

function _normalize_data(
  X::AbstractMatrix{<:Number},
  counts::AbstractVector{<:Real};
  after::Union{Real, Nothing}=nothing,
  copy::Bool=false,
)
  X = copy ? deepcopy(X) : X
  if eltype(X) <: Integer
    X = Float32.(X)
  end
  if isnothing(after)
    counts_pos = counts[counts .> 0]
    after = median(counts_pos)
  end
  counts = counts ./ after
  # call helper function: divide each row (cell) by the corresponding element in counts
  return axis_mul_or_truediv!(X, counts, 0, /, allow_divide_by_zero=false)
end

"""
    normalize_total(
        adata::AnnData;
        target_sum::Union{Real, Nothing}=nothing,
        exclude_highly_expressed::Bool=false,
        max_fraction::Float64=0.05,
        key_added::Union{String, Nothing}=nothing,
        layer::Union{String, Nothing}=nothing,
        layers::Union{String, Vector{String}, Nothing}=nothing,
        layer_norm::Union{String, Nothing}=nothing,
        copy::Bool=false
    ) -> Union{AnnData, Dict{String, Any}}

Normalize total counts per cell to a target sum.

# Arguments
- `adata::AnnData`: The single-cell dataset to be normalized.
- `target_sum::Union{Real, Nothing}=nothing`: Target total count per cell. If `nothing`, normalization is performed to the median of counts per cell.
- `exclude_highly_expressed::Bool=false`: If `true`, highly expressed genes are excluded from the normalization factor calculation.
- `max_fraction::Float64=0.05`: A gene is considered highly expressed if it accounts for more than `max_fraction` of the total counts in any cell.
- `key_added::Union{String, Nothing}=nothing`: If provided, the computed normalization factors are stored in `adata.obs[key_added]`.
- `layer::Union{String, Nothing}=nothing`: The layer to normalize. If `nothing`, the main count matrix (`adata.X`) is used.
- `layers::Union{String, Vector{String}, Nothing}=nothing`: Specifies multiple layers to normalize. If `"all"`, all layers in `adata.layers` are normalized.
- `layer_norm::Union{String, Nothing}=nothing`: Defines the normalization reference for additional layers. Can be `"after"`, `"X"`, or `nothing`.
- `copy::Bool=false`: If `true`, returns a copy of `adata` with normalized counts; otherwise, modifies `adata` in place.

# Returns
- If `copy=true`, returns a new `AnnData` object with normalized counts.
- If `copy=false`, modifies `adata` in place and returns a dictionary containing:
  - `"X"`: The normalized count matrix.
  - `"norm_factor"`: The computed normalization factors.
  - Additional layers if `layers` is specified.

# Notes
- If `exclude_highly_expressed=true`, genes that exceed `max_fraction` in any cell are ignored when computing normalization factors.
- Cells with zero counts will generate a warning.
- The function supports layer-wise normalization if `layers` is specified.

# Example
```julia
adata = AnnData(X)
normalize_total(adata; target_sum=10000)
```
"""
function normalize_total(
  adata::AnnData;
  target_sum::Union{Real, Nothing}=nothing,
  exclude_highly_expressed::Bool=false,
  max_fraction::Float64=0.05,
  key_added::Union{String, Nothing}=nothing,
  layer::Union{String, Nothing}=nothing,
  layers::Union{String, Vector{String}, Nothing}=nothing,
  layer_norm::Union{String, Nothing}=nothing,
  copy::Bool=false,
)
  if copy
    adata = deepcopy(adata)
  end
  if max_fraction < 0 || max_fraction > 1
    error("Choose max_fraction between 0 and 1.")
  end

  # deal with deprecated arguments layers and layer_norm
  if layers !== nothing
    if isa(layers, String) && layers != "all"
      error("`layers` needs to be an array of strings or \"all\", not $(layers)")
    end
  end

  # get the data matrix to be normalized
  x = _get_obs_rep(adata; layer=layer)
  # calc the total count of each cell (sum along rows)
  counts_per_cell = vec(sum(x, dims=2))
  gene_subset::Union{Vector{Bool}, Nothing} = nothing
  msg = "normalizing counts per cell"

  if exclude_highly_expressed
    # if the expression of a gene in any cell exceeds max_fraction of the total count of that cell
    gene_subset = [sum(x[:, j] .> (counts_per_cell .* max_fraction)) == 0 for j in 1:size(x, 2)]
    highly_expressed_genes = adata.var_names[.!gene_subset]
    msg *= ". The following highly-expressed genes are not considered during normalization factor computation:\n$(highly_expressed_genes)"
    # only use genes that are not excluded to recalculate counts_per_cell
    x = x[:, gene_subset]
    counts_per_cell = vec(sum(x, dims=2))
  end

  @info msg

  cell_subset = counts_per_cell .> 0
  if !all(cell_subset)
    @warn "Some cells have zero counts"
  end

  if copy
    if key_added !== nothing
      adata.obs[!, key_added] = counts_per_cell
    end
    _set_obs_rep!(
      adata,
      _normalize_data(x, counts_per_cell; after=target_sum, copy=false);
      layer=layer,
    )
    return adata
  else
    dat = Dict(
      "X" => _normalize_data(x, counts_per_cell; after=target_sum, copy=true),
      "norm_factor" => counts_per_cell,
    )
    after_local::Union{Real, Nothing} = nothing
    if layer_norm == "after"
      after_local = target_sum
    elseif layer_norm == "X"
      after_local = median(counts_per_cell[cell_subset])
    elseif isnothing(layer_norm)
      after_local = nothing
    else
      error("layer_norm should be \"after\", \"X\" or nothing")
    end

    if layers !== nothing
      layers = layers == "all" ? collect(keys(adata.layers)) : layers
      for layer_to_norm::String in layers
        _ = normalize_total(adata; layer=layer_to_norm, target_sum=after_local, inplace=inplace)
        dat[layer_to_norm] = _get_obs_rep(adata; layer=layer_to_norm)
      end
    end

    @info "finished normalization"
    if key_added !== nothing
      @debug "and added $(key_added), counts per cell before normalization (adata.obs)"
    end
    return dat
  end

  @info "finished normalization"
  if key_added !== nothing
    @debug "and added $(key_added), counts per cell before normalization (adata.obs)"
  end

  return adata
end

"""
    normalize_total!(
        adata::AnnData;
        target_sum::Union{Real, Nothing}=nothing,
        exclude_highly_expressed::Bool=false,
        max_fraction::Float64=0.05,
        key_added::Union{String, Nothing}=nothing,
        layer::Union{String, Nothing}=nothing,
        layers::Union{String, Vector{String}, Nothing}=nothing
    ) -> Nothing

Normalize total counts per cell **in-place**.

# Arguments
- `adata::AnnData`: The single-cell dataset to be normalized.
- `target_sum::Union{Real, Nothing}=nothing`: Target total count per cell. If `nothing`, normalization is performed to the median of counts per cell.
- `exclude_highly_expressed::Bool=false`: If `true`, highly expressed genes are excluded from the normalization factor calculation.
- `max_fraction::Float64=0.05`: A gene is considered highly expressed if it accounts for more than `max_fraction` of the total counts in any cell.
- `key_added::Union{String, Nothing}=nothing`: If provided, the computed normalization factors are stored in `adata.obs[key_added]`.
- `layer::Union{String, Nothing}=nothing`: The layer to normalize. If `nothing`, the main count matrix (`adata.X`) is used.
- `layers::Union{String, Vector{String}, Nothing}=nothing`: Specifies multiple layers to normalize. If `"all"`, all layers in `adata.layers` are normalized.

# Returns
- This function **modifies** `adata` in place and does not return a new object.

# Notes
- If `exclude_highly_expressed=true`, genes that exceed `max_fraction` in any cell are ignored when computing normalization factors.
- Cells with zero counts will generate a warning.
- The function supports layer-wise normalization if `layers` is specified.

# Example
```julia
adata = AnnData(X)
normalize_total!(adata; target_sum=10000)
```
"""
function normalize_total!(
  adata::AnnData;
  target_sum::Union{Real, Nothing}=nothing,
  exclude_highly_expressed::Bool=false,
  max_fraction::Float64=0.05,
  key_added::Union{String, Nothing}=nothing,
  layer::Union{String, Nothing}=nothing,
  layers::Union{String, Vector{String}, Nothing}=nothing,
)
  if max_fraction < 0 || max_fraction > 1
    error("Choose max_fraction between 0 and 1.")
  end

  # deal with deprecated arguments layers and layer_norm
  if layers !== nothing
    if isa(layers, String) && layers != "all"
      error("`layers` needs to be an array of strings or \"all\", not $(layers)")
    end
  end

  # get the data matrix to be normalized
  x = _get_obs_rep(adata; layer=layer)
  # calc the total count of each cell (sum along rows)
  counts_per_cell = vec(sum(x, dims=2))
  gene_subset::Union{Vector{Bool}, Nothing} = nothing
  msg = "normalizing counts per cell"

  if exclude_highly_expressed
    # if the expression of a gene in any cell exceeds max_fraction of the total count of that cell
    gene_subset = [sum(x[:, j] .> (counts_per_cell .* max_fraction)) == 0 for j in 1:size(x, 2)]
    highly_expressed_genes = adata.var_names[.!gene_subset]
    msg *= ". The following highly-expressed genes are not considered during normalization factor computation:\n$(highly_expressed_genes)"
    # 仅使用未被排除的基因重新计算 counts_per_cell
    # only use genes that are not excluded to recalculate counts_per_cell
    x = x[:, gene_subset]
    counts_per_cell = vec(sum(x, dims=2))
  end

  @info msg

  cell_subset = counts_per_cell .> 0
  if !all(cell_subset)
    @warn "Some cells have zero counts"
  end

  if key_added !== nothing
    adata.obs[!, key_added] = counts_per_cell
  end
  _set_obs_rep!(
    adata,
    _normalize_data(x, counts_per_cell; after=target_sum, copy=false);
    layer=layer,
  )
end
