#-------------------------------------------------------------------------------------
# AnnData struct
#-------------------------------------------------------------------------------------
#=
"""
    mutable struct AnnData

Minimal Julia implementation of the Python `AnnData` object (see [package documentation](https://anndata.readthedocs.io/en/latest/)
and [Github repository](https://github.com/scverse/anndata)): An `AnnData` object stores a `countmatrix` together with annotations 
of observations `obs` (`obsm`, `obsp`), variables `var` (`varm`, `varp`), and unstructured annotations `uns`.

**Keyword arguments** 
---------------------
 - `countmatrix::Matrix`: countmatrix in (cell x gene) shape
 - `layers::Union{Dict,Nothing}=nothing`: dictionary of other layers (e.g., normalized counts) in the same shape as the countmatrix
 - `obs::Union{DataFrame,Nothing}=nothing`: dataframe of information about cells, e.g., celltypes
 - `obsm::Union{Dict, Nothing}=nothing`: dictionary of observation-level matrices, e.g., a UMAP embedding. The first dimension of the matrix has to correspond to the number of cells.
 - `obsp::Union{Dict, Nothing}=nothing`: dictionary of (observation x observation) matrices, e.g., representing cell graphs. 
 - `var::Union{DataFrame, Nothing}=nothing`: dataframe of information about genes/features, e.g., gene names or highly variable genes
 - `varm::Union{DataFrame, Nothing}=nothing`: dictionary of variable-level matrices. The first dimension of the matrix has to correspond to the number of genes.
 - `obsp::Union{Dict, Nothing}=nothing`: dictionary of (variable x variable) matrices, e.g., representing gene graphs. 
 - `celltypes=nothing`: vector of cell type names, shorthand for `adata.obs["cell_type"]`
 - `uns::Union{Dict, Nothing}=nothing`: dictionary of unstructured annotation. 

"""
Base.@kwdef mutable struct AnnData
    countmatrix::Matrix # shape: cells by genes 
    layers::Union{Dict,Nothing}=nothing
    obs::Union{DataFrame,Nothing}=nothing
    obsm::Union{Dict,Nothing}=nothing
    obsp::Union{Dict,Nothing}=nothing
    var::Union{DataFrame, Nothing}=nothing
    varm::Union{Dict, Nothing}=nothing
    varp::Union{Dict, Nothing}=nothing
    uns::Union{Dict, Nothing}=nothing
    celltypes=nothing
    #is_trained::Bool=false
end
=#

#Base.size(a::AnnData) = size(a.X)
#Base.size(a::AnnData, ind::Int) = size(a.X, ind)
ncells(a::AnnData) = size(a, 1)
ngenes(a::AnnData) = size(a, 2)

"""
    get_celltypes(a::AnnData)

Tries to infer the cell types of cells in an `AnnData` object. 

Returns a vector of cell type names if the cell types are stored in 
`adata.obs["cell_type"]`, `adata.obs["celltype"]`, `adata.obs["celltypes"]`, or `adata.obs["cell_types"]`. 
Otherwise, returns `nothing`.
"""
function get_celltypes(a::AnnData)
  celltypes = nothing
  if !isnothing(a.obs)
    if hasproperty(a.obs, :cell_type)
      celltypes = a.obs.cell_type
    elseif hasproperty(a.obs, :celltype)
      celltypes = a.obs.celltype
    elseif hasproperty(a.obs, :celltypes)
      celltypes = a.obs.celltypes
    elseif hasproperty(a.obs, :cell_types)
      celltypes = a.obs.cell_types
    end
  end
  return celltypes
end

function Base.show(io::IO, a::AnnData)
  #ncells, ngenes = size(adata.X)
  println(
    io,
    "$(typeof(a)) object with a countmatrix with $(ncells(a)) cells and $(ngenes(a)) genes",
  )
  length(a.layers) > 0 && println(io, "   layers dict with the following keys: $(keys(a.layers))")
  nrow(a.obs) > 0 && println(io, "   information about cells: $(first(a.obs,3))")
  nrow(a.var) > 0 && println(io, "   information about genes: $(first(a.var,3))")
  #!hasproperty(a.obs, "celltype") && println(io, "   unique celltypes: $(unique(a.obs.cell_type))")
  #!isnothing(a.celltypes) && println(io, "   unique celltypes: $(unique(a.celltypes))")
  #a.is_trained ? println(io, "    training status: trained") : println(io, "   training status: not trained")
  nothing
end

#=========
# getter #
=========#
#     _get_obs_rep(adata::Muon.AnnData; layer=nothing, obsm=nothing, obsp=nothing)
# 
# Retrieve the observation representation from an AnnData object.
# 
# # Arguments
# - `adata::Muon.AnnData`: The AnnData object containing single-cell data.
# - `layer::Union{Nothing, AbstractString}=nothing`: The specific layer from which to extract data.
# - `obsm::Union{Nothing, AbstractString}=nothing`: The key in `obsm` containing the observation-level representation.
# - `obsp::Union{Nothing, AbstractString}=nothing`: The key in `obsp` containing the pairwise observation-level representation.
# 
# # Returns
# - The corresponding observation representation, either from `adata.X`, `adata.layers`, `adata.obsm`, or `adata.obsp`.
# 
# # Errors
# - Throws an assertion error if more than one of `layer`, `obsm`, or `obsp` is specified.
function _get_obs_rep(
  adata::Muon.AnnData;
  layer::Union{Nothing, AbstractString}=nothing,
  obsm::Union{Nothing, AbstractString}=nothing,
  obsp::Union{Nothing, AbstractString}=nothing,
)
  is_layer = layer !== nothing
  is_obsm = obsm !== nothing
  is_obsp = obsp !== nothing
  choices_made = sum((is_layer, is_obsm, is_obsp))
  @assert choices_made in (0, 1) "only one of layer, use_raw, obsm, obsp can be specified"
  if choices_made == 0
    return adata.X
  elseif is_layer
    return adata.layers[layer]
  elseif is_obsm
    return adata.obsm[obsm]
  elseif is_obsp
    return adata.obsp[obsp]
  end
end

#=========
# setter #
=========#
#     _set_obs_rep!(adata::Muon.AnnData, val::Any; use_raw=false, layer=nothing, obsm=nothing, obsp=nothing)
# 
# Set the observation representation in an AnnData object.
# 
# # Arguments
# - `adata::Muon.AnnData`: The AnnData object to modify.
# - `val::Any`: The new observation representation to set.
# - `use_raw::Bool=false`: If `true`, modifies `adata.raw.X` instead of other fields.
# - `layer::Union{String, Nothing}=nothing`: The specific layer in `adata.layers` to modify.
# - `obsm::Union{String, Nothing}=nothing`: The key in `obsm` to set the observation-level representation.
# - `obsp::Union{String, Nothing}=nothing`: The key in `obsp` to set the pairwise observation-level representation.
# 
# # Returns
# - Nothing, modifies `adata` in place.
# 
# # Errors
# - Throws an error if more than one of `use_raw`, `layer`, `obsm`, or `obsp` is specified.
# - Throws an unexpected error message if an invalid case occurs.
function _set_obs_rep!(
  adata::Muon.AnnData,
  val::Any;
  use_raw::Bool=false,
  layer::Union{String, Nothing}=nothing,
  obsm::Union{String, Nothing}=nothing,
  obsp::Union{String, Nothing}=nothing,
)
  is_layer = !isnothing(layer)
  is_raw_flag = use_raw
  is_obsm = !isnothing(obsm)
  is_obsp = !isnothing(obsp)

  choices_made = (is_layer ? 1 : 0) + (is_raw_flag ? 1 : 0) + (is_obsm ? 1 : 0) + (is_obsp ? 1 : 0)
  if choices_made > 1
    error("Multiple conflicting arguments provided.")
  end

  if choices_made == 0
    adata.X = val
  elseif is_layer
    adata.layers[layer] = val
  elseif is_raw_flag
    adata.raw.X = val
  elseif is_obsm
    adata.obsm[obsm] = val
  elseif is_obsp
    adata.obsp[obsp] = val
  else
    error("Invalid arguments provided.")
  end
  return true
end

# subsetting and copying 

#adata = read_h5ad("cortex_julia_anndata.h5ad")
#adata[1:10, 1:20]
#adata[2,3]

#=
import Base.getindex

function getindex(adata::AnnData, inds...)
    adata_sub = subset_adata(adata, (inds[1], inds[2]))
end

=#
"""
    subset_adata(adata::AnnData, subset_inds::Tuple, dims::Symbol=:both)

Subset an `AnnData` object by indices passed as a tuple of vectors of integers, UnitRanges, or vectors of Booleans.
If `dims` is set to :both, the first element of `subset_inds` is used to subset cells and the second element is used to subset genes.

# Arguments
- `adata`: `AnnData` object to subset
- `subset_inds`: tuple of vectors of integers, UnitRanges, or vectors of Booleans
- `dims`: dimension to subset, either `:cells`, `:genes`, or `:both`

# Returns
- a copy of the `AnnData` object with the subsetted data
"""
function subset_adata(adata::AnnData, subset_inds::Tuple, dims::Symbol=:both)
  adata_new = deepcopy(adata)
  subset_adata!(adata_new, subset_inds, dims)
  return adata_new
end

"""
    subset_adata(adata::AnnData, subset_inds::Union{Int, Vector{Int}, UnitRange, Vector{Bool}}, dims::Symbol)

Subset an `AnnData` object by indices passed as a vector of integers or booleans or as a UnitRange. 
The `dims` argument can be set to either `:cells` or `:genes` to specify which dimension to subset. 

# Arguments
- `adata`: `AnnData` object to subset
- `subset_inds`: vector of integers or booleans or UnitRange
- `dims`: dimension to subset, either `:cells` or `:genes`

# Returns
- a copy of the `AnnData` object with the subsetted data
"""
function subset_adata(
  adata::AnnData,
  subset_inds::Union{Int, Vector{Int}, UnitRange, Vector{Bool}},
  dims::Symbol,
)
  adata_new = deepcopy(adata)
  subset_adata!(adata_new, subset_inds, dims)
  return adata_new
end

"""
    subset_adata!(adata::AnnData, subset_inds, dims::Symbol)

In-place version of `subset_adata`, see `?subset_adata` for more details. 

For `subset_inds`, either a tuple of vectors or ranges can be passed with `dims` set to :both, for subsetting 
both cells and genes, or a single vector or range can be passed with `dims` set to either :cells or :genes.

# Arguments
- `adata`: `AnnData` object to subset
- `subset_inds`: tuple of vectors of integers, UnitRanges, or vectors of Booleans or vector of integers or booleans or UnitRange
- `dims`: dimension to subset, either `:cells`, `:genes`, or `:both`

# Returns
- the `AnnData` object with the subsetted data
"""
function subset_adata!(adata::AnnData, subset_inds, dims::Symbol)
  return subset_adata!(adata, subset_inds, Val(dims))
end

function subset_adata!(adata::AnnData, subset_inds::Tuple, ::Val{:both})
  subset_adata!(adata, subset_inds[1], :cells)
  subset_adata!(adata, subset_inds[2], :genes)
  return adata
end

subset_adata!(adata::AnnData, subset_inds::Tuple, ::Val{:cells}) =
  subset_adata!(adata, subset_inds[1], :cells)

subset_adata!(adata::AnnData, subset_inds::Tuple, ::Val{:genes}) =
  subset_adata!(adata, subset_inds[2], :genes)

#subset_adata!(adata::AnnData, subset_inds::Union{Int, Vector{Int}, UnitRange, BitVector}, ::Val{:cells}) = adata[subset_inds, :]

#subset_adata!(adata::AnnData, subset_inds::Union{Int, Vector{Int}, UnitRange, BitVector}, ::Val{:genes}) = adata[:, subset_inds]

function subset_adata!(
  adata::AnnData,
  subset_inds::Union{Int, Vector{Int}, Vector{Bool}, UnitRange},
  ::Val{:cells},
)
  #adata.ncells = length(subset_inds)

  adata.X = adata.X[subset_inds, :]

  adata.obs_names = adata.obs_names[subset_inds]

  if nrow(adata.obs) > 0 #&& nrow(adata.var) > 0
    adata.obs = adata.obs[subset_inds, :]
  end

  if length(adata.layers) > 0
    for key in keys(adata.layers)
      adata.layers[key] = setindex!(adata.layers, adata.layers[key][subset_inds, :], key)
    end
  end

  if length(adata.obsm) > 0
    for key in keys(adata.obsm)
      adata.obsm[key] = adata.obsm[key][subset_inds, :]
    end
  end

  if length(adata.obsp) > 0
    for key in keys(adata.obsp)
      adata.obsp[key] = adata.obsp[key][subset_inds, subset_inds]
    end
  end

  return adata
end

function subset_adata!(
  adata::AnnData,
  subset_inds::Union{Int, Vector{Int}, Vector{Bool}, UnitRange},
  ::Val{:genes},
)
  adata.X = adata.X[:, subset_inds]

  if length(adata.var_names) > 0
    adata.var_names = adata.var_names[subset_inds]
  end

  if nrow(adata.var) > 0 #&& nrow(adata.var) > 0
    adata.var = adata.var[subset_inds, :]
  end

  if length(adata.layers) > 0
    for key in keys(adata.layers)
      adata.layers[key] = setindex!(adata.layers, adata.layers[key][:, subset_inds], key)
    end
  end

  if length(adata.varm) > 0
    for key in keys(adata.varm)
      adata.varm[key] = adata.varm[key][:, subset_inds]
    end
  end

  if length(adata.varp) > 0
    for key in keys(adata.varp)
      adata.varp[key] = adata.varp[key][subset_inds, subset_inds]
    end
  end

  return adata
end

"""
    insert_obs!(adata::AnnData, obs::DataFrame)

Insert new columns into the `obs` DataFrame of an `AnnData` object.

# Arguments
- `adata`: The `AnnData` object to modify.
- `obs`: A DataFrame containing the new columns to be inserted. The number of rows in `obs` must match the number of cells in `adata`.

# Returns
- `nothing`: Modifies the `adata.obs` DataFrame in place.

# Notes
- This function assumes that the number of rows in `obs` matches the number of cells in `adata`. If they do not match, an error will be thrown.
- This function modifies the `adata.obs` DataFrame directly. If you want to keep the original `adata` unchanged, make a copy before calling this function.
"""
function insert_obs!(adata::AnnData, obs::DataFrame)
  for col in names(obs)
    adata.obs[!, col] = obs[!, col]
  end
end

"""
    insert_var!(adata::AnnData, var::DataFrame)

Insert new columns into the `var` DataFrame of an `AnnData` object.

# Arguments
- `adata`: The `AnnData` object to modify.
- `var`: A DataFrame containing the new columns to be inserted. The number of rows in `var` must match the number of genes in `adata`.

# Returns
- `nothing`: Modifies the `adata.var` DataFrame in place.

# Notes
- This function assumes that the number of rows in `var` matches the number of genes in `adata`. If they do not match, an error will be thrown.
- This function modifies the `adata.var` DataFrame directly. If you want to keep the original `adata` unchanged, make a copy before calling this function.
"""
function insert_var!(adata::AnnData, var::DataFrame)
  for col in names(var)
    adata.var[!, col] = var[!, col]
  end
end
