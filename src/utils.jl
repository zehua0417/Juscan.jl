#=========================================
filename: utils
author: Lihuax
date: 2025-03-02 17:34
description: This file contains utility functions for performing common operations on sparse matrices, particularly designed to support single-cell or multi-omics data analysis. It includes functions for counting non-zero elements, computing sums, and calculating means along specified axes.
=========================================#

# Function to count non-zero elements along a specified axis
function axis_nnz(X::SparseArrays.SparseMatrixCSC{T, I}; axis::Int) where {T <: Real, I <: Integer}
  if axis == 1 || axis == 2
    return sum(X .!= 0, dims=axis)
  else
    throw(ArgumentError("axis must be 1 or 2"))
  end
end

function axis_sum(X::SparseArrays.SparseMatrixCSC{T, I}; axis::Int) where {T <: Real, I <: Integer}
  if axis == 1 || axis == 2
    return sum(X, dims=axis)
  else
    throw(ArgumentError("axis must be 1 or 2"))
  end
end

function axis_mean(X::SparseArrays.SparseMatrixCSC{T, I}; axis::Int) where {T <: Real, I <: Integer}
  if axis == 1 || axis == 2
    return sum(X, dims=axis) ./ axis_nnz(X, axis=axis)
  else
    throw(ArgumentError("axis must be 1 or 2"))
  end
end
