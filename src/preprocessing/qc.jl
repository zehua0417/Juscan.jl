#=========================================
filename: qc
author: Lihuax
date: 2025-03-02 15:37
description: This file contains functions for calculating quality control (QC) metrics for single-cell or multi-omics data. It includes utilities for computing per-cell and per-gene metrics, as well as support functions for matrix manipulation and QC calculations.
=========================================#

"""
    calculate_qc_metrics!(adata; expr_type="counts", var_type="genes", qc_vars=String[], 
                          percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, 
                          use_log1p=true, parallel=nothing) -> Tuple{DataFrame, DataFrame}

Calculate quality control metrics and return.

# Arguments
- `adata`: `Muon.AnnData`  
  Annotated data matrix containing single-cell or multi-omics data.
- `expr_type`: `String` (default: `"counts"`)  
  The type of expression data to use, such as raw counts or normalized values.
- `var_type`: `String` (default: `"genes"`)  
  Specifies the variable type, e.g., `"genes"` for gene expression.
- `qc_vars`: `Union{Vector{String}, String}` (default: `String[]`)  
  A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.
- `percent_top`: `Union{Vector{Int}, Nothing}` (default: `[50, 100, 200, 500]`)  
  A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.
- `layer`: `Union{String, Nothing}` (default: `nothing`)  
  Specifies the layer in `adata` to use. If `nothing`, the default matrix is used.
- `use_raw`: `Bool` (default: `false`)  
  Whether to use the raw counts matrix (`adata.raw`) instead of the processed data.
- `use_log1p`: `Bool` (default: `true`)  
  Whether to apply `log1p` transformation to the data before computing quality control metrics.
- `parallel`: `Union{Bool, Nothing}` (default: `nothing`)  
  **(Deprecated)** This argument is ignored but retained for backward compatibility.

# Returns
- A tuple of two `DataFrames.DataFrame`:  
  1. The first `DataFrame` contains per-cell quality control metrics.  
  2. The second `DataFrame` contains per-gene quality control metrics.
"""
function calculate_qc_metrics(
  adata::Muon.AnnData;
  expr_type::String="counts",
  var_type::String="genes",
  qc_vars::Union{Vector{String}, String}=String[],
  percent_top::Union{Vector{Int}, Nothing}=[50, 100, 200, 500],
  layer::Union{String, Nothing}=nothing,
  use_raw::Bool=false,
  use_log1p::Bool=true,
  parallel::Union{Bool, Nothing}=nothing,
)::Tuple{DataFrames.DataFrame, DataFrames.DataFrame}
  if !isnothing(parallel)
    @warn("Argument `parallel` is deprecated, and currently has no effect.")
  end
  # get X matrix
  X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
  if X isa
     LinearAlgebra.Adjoint{T, SparseArrays.SparseMatrixCSC{T, I}} where {T <: Real, I <: Integer}
    X = SparseArrays.sparse(X)
  end
  if isa(qc_vars, String)
    qc_vars = [qc_vars]
  end

  obs_metrics = describe_obs(
    adata, # Muon.AnnData
    expr_type=expr_type, # String
    var_type=var_type,
    qc_vars=qc_vars,
    percent_top=percent_top,
    X=X,
    use_log1p=use_log1p,
  )

  var_metrics = describe_var(
    adata, # Muon.AnnData
    expr_type=expr_type, # String
    var_type=var_type,
    X=X,
    use_log1p=use_log1p,
  )
  return (obs_metrics, var_metrics)
end

"""
    calculate_qc_metrics!(adata; expr_type="counts", var_type="genes", qc_vars=String[], 
                          percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, 
                          use_log1p=true, parallel=nothing) -> Nothing

Calculate and store quality control (QC) metrics directly in the `obs` and `var` attributes of `adata`.

# Arguments
- `adata`: `Muon.AnnData`  
  Annotated data matrix containing single-cell or multi-omics data.  
  The QC metrics will be stored in `adata.obs` (per-cell) and `adata.var` (per-gene).
- `expr_type`: `String` (default: `"counts"`)  
  The type of expression data to use, such as raw counts or normalized values.
- `var_type`: `String` (default: `"genes"`)  
  Specifies the variable type, e.g., `"genes"` for gene expression.
- `qc_vars`: `Union{Vector{String}, String}` (default: `String[]`)  
  A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.
- `percent_top`: `Union{Vector{Int}, Nothing}` (default: `[50, 100, 200, 500]`)  
  A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.
- `layer`: `Union{String, Nothing}` (default: `nothing`)  
  Specifies the layer in `adata` to use. If `nothing`, the default matrix is used.
- `use_raw`: `Bool` (default: `false`)  
  Whether to use the raw counts matrix (`adata.raw`) instead of the processed data.
- `use_log1p`: `Bool` (default: `true`)  
  Whether to apply `log1p` transformation to the data before computing quality control metrics.
- `parallel`: `Union{Bool, Nothing}` (default: `nothing`)  
  **(Deprecated)** This argument is ignored but retained for backward compatibility.

# Returns
- `Nothing`:  
  This function modifies `adata` in place by adding QC metrics to `adata.obs` and `adata.var`.

# Notes
- Unlike `calculate_qc_metrics`, which returns QC metrics as `DataFrame` objects,  
  this function directly updates `adata` without returning any values.
"""
function calculate_qc_metrics!(
  adata::Muon.AnnData;
  expr_type::String="counts",
  var_type::String="genes",
  qc_vars::Union{Vector{String}, String}=String[],
  percent_top::Union{Vector{Int}, Nothing}=[50, 100, 200, 500],
  layer::Union{String, Nothing}=nothing,
  use_raw::Bool=false,
  use_log1p::Bool=true,
  parallel::Union{Bool, Nothing}=nothing,
)::Nothing
  if !isnothing(parallel)
    @warn("Argument `parallel` is deprecated, and currently has no effect.")
  end
  # get X matrix
  X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
  if X isa
     LinearAlgebra.Adjoint{T, SparseArrays.SparseMatrixCSC{T, I}} where {T <: Real, I <: Integer}
    X = SparseArrays.sparse(X)
  end
  if isa(qc_vars, String)
    qc_vars = [qc_vars]
  end

  describe_obs!(
    adata, # Muon.AnnData
    expr_type=expr_type, # String
    var_type=var_type,
    qc_vars=qc_vars,
    percent_top=percent_top,
    X=X,
    use_log1p=use_log1p,
  )

  describe_var!(
    adata, # Muon.AnnData
    expr_type=expr_type, # String
    var_type=var_type,
    X=X,
    use_log1p=use_log1p,
  )
end

function _choose_mtx_rep(
  adata::Muon.AnnData;
  use_raw::Bool=false,
  layer::Union{String, Nothing}=nothing,
)::Union{
  LinearAlgebra.Adjoint{
    T,
    SparseArrays.SparseArrays.SparseMatrixCSC{T, I},
  } where {T <: Real, I <: Integer},
  SparseArrays.SparseMatrixCSC{T, I} where {T <: Real, I <: Integer},
}
  is_layer = !isnothing(layer)

  if use_raw && is_layer
    throw(
      ArgumentError(
        "Cannot use expression from both layer and raw. You provided: use_raw=\$(use_raw) and layer=\$(layer)",
      ),
    )
  end

  if is_layer
    return adata.layers[layer]
  elseif use_raw
    return adata.raw.X
  else
    return adata.X
  end
end

"""
    describe_obs(adata; expr_type="counts", var_type="genes", qc_vars=String[], 
                 percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, 
                 use_log1p=true, X=nothing, parallel=nothing) -> DataFrames.DataFrame

Compute per-cell quality control (QC) metrics from the expression matrix and return them as a `DataFrame`.

# Arguments
- `adata`: `Muon.AnnData`  
  Annotated data matrix containing single-cell or multi-omics data.
- `expr_type`: `String` (default: `"counts"`)  
  The type of expression data to use, such as raw counts or normalized values.
- `var_type`: `String` (default: `"genes"`)  
  Specifies the variable type, e.g., `"genes"` for gene expression.
- `qc_vars`: `Union{Vector{String}, String}` (default: `String[]`)  
  A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.
- `percent_top`: `Union{Vector{Int}, Nothing}` (default: `[50, 100, 200, 500]`)  
  A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.
- `layer`: `Union{String, Nothing}` (default: `nothing`)  
  Specifies the layer in `adata` to use. If `nothing`, the default matrix is used.
- `use_raw`: `Bool` (default: `false`)  
  Whether to use the raw counts matrix (`adata.raw`) instead of the processed data.
- `use_log1p`: `Bool` (default: `true`)  
  Whether to apply `log1p` transformation to the data before computing quality control metrics.
- `X`: `Union{Nothing, SparseMatrixCSC}` (default: `nothing`)  
  The expression matrix to use. If `nothing`, the matrix is determined by `_choose_mtx_rep`.
- `parallel`: `Union{Bool, Nothing}` (default: `nothing`)  
  **(Deprecated)** This argument is ignored but retained for backward compatibility.

# Returns
- `DataFrames.DataFrame`:  
  A `DataFrame` containing per-cell QC metrics, including:
  - `"n_\$(var_type)_by_\$(expr_type)"`: Number of detected features per cell.
  - `"log1p_n_\$(var_type)_by_\$(expr_type)"`: Log-transformed number of features per cell.
  - `"total_\$(expr_type)"`: Total expression count per cell.
  - `"log1p_total_\$(expr_type)"`: Log-transformed total expression count per cell.
  - `"pct_\$(expr_type)_in_top_N_\$(var_type)"`: Percentage of expression from the top `N` features.
  - `"total_\$(expr_type)_\$(qc_var)"`: Total expression for each QC variable per cell.
  - `"log1p_total_\$(expr_type)_\$(qc_var)"`: Log-transformed total expression for each QC variable per cell.
  - `"pct_\$(expr_type)_\$(qc_var)"`: Percentage of total expression contributed by each QC variable.

# Notes
- This function computes per-cell QC metrics based on expression data.  
  For per-gene metrics, see `describe_var`.
"""
function describe_obs(
  adata::Muon.AnnData;
  expr_type::String="counts",
  var_type::String="genes",
  qc_vars::Union{Vector{String}, String}=String[],
  percent_top::Union{Vector{Int}, Nothing}=[50, 100, 200, 500],
  layer::Union{String, Nothing}=nothing,
  use_raw::Bool=false,
  use_log1p::Bool=true,
  X::Union{
    Nothing,
    SparseArrays.SparseArrays.SparseMatrixCSC{T, I} where {T <: Real, I <: Integer},
  }=nothing,
  parallel::Union{Bool, Nothing}=nothing,
)::Union{DataFrames.DataFrame, Nothing}
  if !isnothing(parallel)
    @warn("Argument `parallel` is deprecated, and currently has no effect.")
  end
  if isnothing(X)
    X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
    if X isa
       LinearAlgebra.Adjoint{T, SparseArrays.SparseMatrixCSC{T, I}} where {T <: Real, I <: Integer}
      X = SparseArrays.sparse(X)
    end
  end
  obs_metrics = DataFrames.DataFrame(index=adata.obs_names)

  # (logip) non zero counts for each cell
  obs_metrics[!, "n_$(var_type)_by_$(expr_type)"] = vec(axis_nnz(X, axis=2))
  if use_log1p
    obs_metrics[!, "log1p_n_$(var_type)_by_$(expr_type)"] =
      log1p.(obs_metrics[!, "n_$(var_type)_by_$(expr_type)"])
  end

  # (log1p) total counts for each cell
  obs_metrics[!, "total_$(expr_type)"] = vec(axis_sum(X, axis=2))
  if use_log1p
    obs_metrics[!, "log1p_total_$(expr_type)"] = log1p.(obs_metrics[!, "total_$(expr_type)"])
  end

  # the percentage of total expression that comes from the top n most highly expressed genes in each cell
  if !isnothing(percent_top)
    percent_top = sort(percent_top)
    # Cumulative proportions of the top [percent_top] genes in each sample
    proportions = top_segment_proportions(X, percent_top)
    for (i, n) in enumerate(percent_top)
      obs_metrics[!, "pct_$(expr_type)_in_top_$(n)_$(var_type)"] = proportions[:, i] .* 100
    end
  end

  for qc_var in qc_vars
    # (log1p) total counts for each cell for each qc_var
    obs_metrics[!, "total_$(expr_type)_$(qc_var)"] =
      axis_sum(X[adata.var[!, qc_var] .== true, :], 1)
    if use_log1p
      obs_metrics[!, "log1p_total_$(expr_type)_$(qc_var)"] =
        log1p.(obs_metrics[!, "total_$(expr_type)_$(qc_var)"])
    end
    # the percentage of total expression that comes from a specific type of gene in each cell
    obs_metrics[!, "pct_$(expr_type)_$(qc_var)"] =
      obs_metrics[!, "total_$(expr_type)_$(qc_var)"] ./ obs_metrics[!, "total_$(expr_type)"] .* 100
  end
  return obs_metrics
end

"""
    describe_obs!(adata; expr_type="counts", var_type="genes", qc_vars=String[], 
                  percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, 
                  use_log1p=true, X=nothing, parallel=nothing) -> Nothing

Compute per-cell quality control (QC) metrics and store them in `adata.obs`.

# Arguments
- `adata`: `Muon.AnnData`  
  Annotated data matrix containing single-cell or multi-omics data.
- `expr_type`: `String` (default: `"counts"`)  
  The type of expression data to use, such as raw counts or normalized values.
- `var_type`: `String` (default: `"genes"`)  
  Specifies the variable type, e.g., `"genes"` for gene expression.
- `qc_vars`: `Union{Vector{String}, String}` (default: `String[]`)  
  A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.
- `percent_top`: `Union{Vector{Int}, Nothing}` (default: `[50, 100, 200, 500]`)  
  A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.
- `layer`: `Union{String, Nothing}` (default: `nothing`)  
  Specifies the layer in `adata` to use. If `nothing`, the default matrix is used.
- `use_raw`: `Bool` (default: `false`)  
  Whether to use the raw counts matrix (`adata.raw`) instead of the processed data.
- `use_log1p`: `Bool` (default: `true`)  
  Whether to apply `log1p` transformation to the data before computing quality control metrics.
- `X`: `Union{Nothing, SparseMatrixCSC}` (default: `nothing`)  
  The expression matrix to use. If `nothing`, the matrix is determined by `_choose_mtx_rep`.
- `parallel`: `Union{Bool, Nothing}` (default: `nothing`)  
  **(Deprecated)** This argument is ignored but retained for backward compatibility.

# Returns
- `Nothing`:  
  This function modifies `adata` in place by adding per-cell QC metrics to `adata.obs`.

# Notes
- This function is similar to [`describe_obs`](@ref) but modifies `adata` directly.  
  It computes QC metrics using `describe_obs` and stores the results in `adata.obs`.  
  - For per-gene QC metrics, see [`describe_var!`](@ref).
"""
function describe_obs!(
  adata::Muon.AnnData;
  expr_type::String="counts",
  var_type::String="genes",
  qc_vars::Union{Vector{String}, String}=String[],
  percent_top::Union{Vector{Int}, Nothing}=[50, 100, 200, 500],
  layer::Union{String, Nothing}=nothing,
  use_raw::Bool=false,
  use_log1p::Bool=true,
  X::Union{
    Nothing,
    SparseArrays.SparseArrays.SparseMatrixCSC{T, I} where {T <: Real, I <: Integer},
  }=nothing,
  parallel::Union{Bool, Nothing}=nothing,
)::Nothing
  obs_metrics = describe_obs(
    adata,
    expr_type=expr_type,
    var_type=var_type,
    qc_vars=qc_vars,
    percent_top=percent_top,
    layer=layer,
    use_raw=use_raw,
    use_log1p=use_log1p,
    X=X,
    parallel=parallel,
  )
  insert_obs!(adata, obs_metrics)
end

"""
    describe_var(adata; expr_type="counts", var_type="genes", layer=nothing, 
                 use_raw=false, use_log1p=true, X=nothing) -> DataFrame

Compute per-gene quality control (QC) metrics and return as a DataFrame.

# Arguments
- `adata`: `Muon.AnnData`  
  Annotated data matrix containing single-cell or multi-omics data.
- `expr_type`: `String` (default: `"counts"`)  
  The type of expression data to use, such as raw counts or normalized values.
- `var_type`: `String` (default: `"genes"`)  
  Specifies the variable type, e.g., `"genes"` for gene expression.
- `layer`: `Union{String, Nothing}` (default: `nothing`)  
  Specifies the layer in `adata` to use. If `nothing`, the default matrix is used.
- `use_raw`: `Bool` (default: `false`)  
  Whether to use the raw counts matrix (`adata.raw`) instead of the processed data.
- `use_log1p`: `Bool` (default: `true`)  
  Whether to apply `log1p` transformation to the data before computing quality control metrics.
- `X`: `Union{Nothing, SparseMatrixCSC}` (default: `nothing`)  
  The expression matrix to use. If `nothing`, the matrix is determined by `_choose_mtx_rep`.

# Returns
- `DataFrames.DataFrame`:  
  A DataFrame where each row represents a gene (or variable), containing the following QC metrics:
  - `"n_cells_by_<expr_type>"`: Number of cells in which each gene is detected.
  - `"mean_<expr_type>"`: Mean expression level of each gene.
  - `"log1p_mean_<expr_type>"`: Log-transformed mean expression level.
  - `"pct_dropout_by_<expr_type>"`: Percentage of cells in which each gene is not detected.
  - `"total_<expr_type>"`: Total expression level of each gene.
  - `"log1p_total_<expr_type>"`: Log-transformed total expression level.

# Notes
- This function calculates quality control metrics per gene (or other variables).  
  - To store these metrics directly in `adata.var`, use [`describe_var!`](@ref).  
  - To compute per-cell QC metrics, see [`describe_obs`](@ref).
"""
function describe_var(
  adata::Muon.AnnData;
  expr_type::String="counts",
  var_type::String="genes",
  layer::Union{String, Nothing}=nothing,
  use_raw::Bool=false,
  use_log1p::Bool=true,
  X::Union{
    Nothing,
    SparseArrays.SparseArrays.SparseMatrixCSC{T, I} where {T <: Real, I <: Integer},
  }=nothing,
)::DataFrames.DataFrame
  if isnothing(X)
    X = _choose_mtx_rep(adata, use_raw=use_raw, layer=layer)
    if X isa
       LinearAlgebra.Adjoint{T, SparseArrays.SparseMatrixCSC{T, I}} where {T <: Real, I <: Integer}
      X = SparseArrays.sparse(X)
    end
  end

  var_metrics = DataFrames.DataFrame(index=adata.var_names)

  # number of cells in which each gene is detected
  var_metrics[!, "n_cells_by_$(expr_type)"] = vec(axis_nnz(X, axis=1))

  # (log1p) mean expression of each gene
  var_metrics[!, "mean_$(expr_type)"] = vec(_get_mean_var(X, axis=2)[1])
  if use_log1p
    var_metrics[!, "log1p_mean_$(expr_type)"] = log1p.(var_metrics[!, "mean_$(expr_type)"])
  end

  # the percentage of cells in which each gene is not detected
  var_metrics[!, "pct_dropout_by_$(expr_type)"] =
    (1 .- var_metrics[!, "n_cells_by_$(expr_type)"] ./ size(X, 1)) .* 100

  # (log1p) total expression of each gene
  var_metrics[!, "total_$(expr_type)"] = vec(axis_sum(X, axis=1))
  if use_log1p
    var_metrics[!, "log1p_total_$(expr_type)"] = log1p.(var_metrics[!, "total_$(expr_type)"])
  end

  return var_metrics
end

"""
    describe_var!(adata; expr_type="counts", var_type="genes", layer=nothing, 
                  use_raw=false, use_log1p=true, X=nothing) -> Nothing

Compute per-gene quality control (QC) metrics and store them in `adata.var`.

# Arguments
- `adata`: `Muon.AnnData`  
  Annotated data matrix containing single-cell or multi-omics data.
- `expr_type`: `String` (default: `"counts"`)  
  The type of expression data to use, such as raw counts or normalized values.
- `var_type`: `String` (default: `"genes"`)  
  Specifies the variable type, e.g., `"genes"` for gene expression.
- `layer`: `Union{String, Nothing}` (default: `nothing`)  
  Specifies the layer in `adata` to use. If `nothing`, the default matrix is used.
- `use_raw`: `Bool` (default: `false`)  
  Whether to use the raw counts matrix (`adata.raw`) instead of the processed data.
- `use_log1p`: `Bool` (default: `true`)  
  Whether to apply `log1p` transformation to the data before computing quality control metrics.
- `X`: `Union{Nothing, SparseMatrixCSC}` (default: `nothing`)  
  The expression matrix to use. If `nothing`, the matrix is determined by `_choose_mtx_rep`.

# Returns
- `Nothing`:  
  This function modifies `adata` in place by adding per-gene QC metrics to `adata.var`.

# Notes
- This function is similar to [`describe_var`](@ref) but modifies `adata` directly.  
  It computes QC metrics using `describe_var` and stores the results in `adata.var`.  
  - For per-cell QC metrics, see [`describe_obs!`](@ref).
"""
function describe_var!(
  adata::Muon.AnnData;
  expr_type::String="counts",
  var_type::String="genes",
  layer::Union{String, Nothing}=nothing,
  use_raw::Bool=false,
  use_log1p::Bool=true,
  X::Union{
    Nothing,
    SparseArrays.SparseArrays.SparseMatrixCSC{T, I} where {T <: Real, I <: Integer},
  }=nothing,
)::Nothing
  var_metrics = describe_var(
    adata,
    expr_type=expr_type,
    var_type=var_type,
    layer=layer,
    use_raw=use_raw,
    use_log1p=use_log1p,
    X=X,
  )
  insert_var!(adata, var_metrics)
end

#==
Calculate the proportion of top segment values for each column in a sparse matrix.

This function takes a sparse matrix and a list of segment sizes (ns), then computes
the cumulative proportion of the top values in each column for each segment size.
The result is normalized by the total column sum.

# Arguments:
- mtx::SparseArrays.SparseMatrixCSC{T, I}: The input sparse matrix.
- ns::Union{Vector{<:Integer}, Nothing}: A list of segment sizes (e.g., [50, 100, 200]).
  If nothing is provided, the default list [50, 100, 200, 500] is used.

# Returns:
- Matrix{T}: A matrix where each column corresponds to the cumulative proportion
  of the top segment values for each segment size.

# Notes:
- The function assumes that the input matrix is sparse and that `ns` is sorted in ascending order.
- The result is normalized by the column sums, so each value represents a proportion.
==#
function top_segment_proportions(
  mtx::SparseArrays.SparseMatrixCSC{T, I},
  ns::Union{Vector{<:Integer}, Nothing}=[50, 100, 200, 500],
)::Matrix{T} where {T <: Real, I <: Integer}
  # Ensure ns is sorted in ascending order
  ns = sort(ns)
  num_rows, num_cols = size(mtx)

  # Calculate the sum of each row
  row_sums = sum(mtx, dims=2)[:]

  # Initialize the result matrix
  values = zeros(T, num_rows, length(ns))

  # For each row, extract the top `ns[-1]` values
  for j in 1:num_rows
    row = mtx[j, :] |> collect  # Convert sparse column to dense vector
    sorted_row = sort(row, rev=true)  # Sort in descending order
    acc = zero(T)  # Accumulator for cumulative sum
    prev = 0
    for (i, n) in enumerate(ns)
      if num_cols < n
        throw(ArgumentError("The number of rows in the matrix is less than the segment size"))
      end
      acc += sum(sorted_row[(prev + 1):n])  # Sum the top `n` values
      values[j, i] = acc
      prev = n
    end
  end

  # Normalize the result by column sums
  return values ./ row_sums
end
