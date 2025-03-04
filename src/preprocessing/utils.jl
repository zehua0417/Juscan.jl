#=========================================
filename: utils
author: Lihuax
date: 2025-03-02 15:49
description: This file contains utility functions for data manipulation and statistical calculations, primarily designed to support single-cell or multi-omics data analysis. It includes functions for computing means and variances along specified axes of sparse matrices.
=========================================#

function _get_mean_var(
  X::SparseArrays.SparseMatrixCSC{T, I};
  axis::Int=1,
)::Tuple{LinearAlgebra.Matrix, LinearAlgebra.Matrix} where {T <: Real, I <: Integer}
  if axis == 1
    # calc mean & mean_sq for each row
    n = size(X, 2)
    mean = sum(X, dims=2) ./ n
    mean_sq = sum(X .* X, dims=2) ./ n
  elseif axis == 2
    # calc mean & mean_sq for each column
    n = size(X, 1)
    mean = sum(X, dims=1) ./ n
    mean_sq = sum(X .* X, dims=1) ./ n
  else
    throw(ArgumentError("axis must be 1 (row-wise) or 2 (column-wise)"))
  end
  var = mean_sq .- mean .^ 2

  # Unbiased estimate adjustment
  if n > 1
    var .= var .* (n / (n - 1))
  end

  return mean, var
end
