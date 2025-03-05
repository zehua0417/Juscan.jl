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

function check_op(op::Function)
  if op !== (*) && op !== (/)
    error("op must be either multiplication (*) or division (/)")
  end
end

# helper function to broadcast a 1D array to a shape that can be broadcasted
function broadcast_axis(arr::AbstractVector, axis::Int)
  if axis == 0
    return reshape(arr, :, 1)
  elseif axis == 1
    return reshape(arr, 1, :)
  else
    error("axis must be 0 or 1")
  end
end

# calc X * scaling_array or X / scaling_array along the specified axis
function axis_mul_or_truediv!(
  X::AbstractArray{<:Number},
  scaling_array::AbstractArray{<:Number},
  axis::Int,
  op::Function;
  allow_divide_by_zero::Bool=true,
  out::Union{AbstractArray, Nothing}=nothing,
)
  check_op(op)
  scaling_array_b = broadcast_axis(scaling_array, axis)
  if op === (*)
    result = X .* scaling_array_b
    if out !== nothing
      out .= result
      return out
    else
      return result
    end
  elseif op === (/)
    if !allow_divide_by_zero
      # for each element, add 1 if it is 0 to avoid division by zero
      scaling_array_b = copy(scaling_array_b) .+ Int.(scaling_array_b .== 0)
    end
    result = X ./ scaling_array_b
    if out !== nothing
      out .= result
      return out
    else
      return result
    end
  else
    error("op must be either multiplication (*) or division (/)")
  end
end
