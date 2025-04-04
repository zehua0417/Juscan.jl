#=========================================
filename: snn
author: Lihuax
date: 2025-03-21 11:28
description:
=========================================#

function SNN(nn_rank::AbstractMatrix, prune::AbstractFloat)
  rn, cn = size(nn_rank)
  snn = zeros(Float64, rn, rn)

  @inbounds @fastmath @simd for idx in CartesianIndices(nn_rank)
    snn[idx[1], Int(nn_rank[idx])] = 1
  end

  snn = snn * transpose(snn)
  index = snn .!= 0
  snn[index] .= TransPrune.(snn[index], cn, prune)

  return snn |> SparseArrays.sparse
end

function TransPrune(x, cn::Integer, prune::AbstractFloat)::AbstractFloat
  y = x / (2 * cn - x)
  return y < prune ? 0 : y
end
