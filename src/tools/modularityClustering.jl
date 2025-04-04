#=========================================
filename: ModularityClustering
author: Lihuax
date: 2025-03-21 11:29
description:
=========================================#

function ModClustering(
  nn::SparseArrays.AbstractSparseMatrix;
  modularity_fun::Integer=1,
  resolution::Real=0.8,
  algorithm::Integer=1,
  rsn::Integer=10,
  itn::Integer=10,
  seed::Integer=0,
)
  nn = LinearAlgebra.LowerTriangular(nn) |> SparseArrays.sparse
  @inbounds @fastmath @simd for i in 1:(nn.n)
    nn[i, i] = 0
  end
  SparseArrays.dropzeros!(nn)
  node2, node1, edge_weights = SparseArrays.findnz(nn)

  nodes_number = maximum([nn.m, nn.n])
  network = MatrixToNetwork(node1, node2, Float32.(edge_weights), modularity_fun, nodes_number)

  resolution2 =
    modularity_fun == 1 ?
    (resolution / (2 * TotalEdgeWeights(network) + network.total_weights_self)) : resolution

  # Community detection
  max_modularity = -Inf
  random = JavaRandom(seed)
  Seed!(random, seed)
  clustering = nothing
  modularity = nothing
  @inbounds for i in 1:rsn
    vosct = VOSct(network, resolution2)
    j = 1
    update = true
    while j <= itn && update
      if algorithm == 1
        update = Louvain!(vosct, random)
      else
        @warn "Now only support louvain modularity clustering!"
        update = Louvain!(vosct, random)
      end
      j += 1
      modularity = Quality(vosct)
    end

    if modularity > max_modularity
      clustering = vosct.clustering
      max_modularity = modularity
    end
  end
  @info "Max modularity: $max_modularity"

  return clustering
end
