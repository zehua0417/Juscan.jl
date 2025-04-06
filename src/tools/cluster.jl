#=========================================
filename: cluster
author: Lihuax
date: 2025-03-19 14:50
description:
=========================================#
using Clustering
using GLM
using StatsModels
using Graphs, NearestNeighbors, SimpleWeightedGraphs
using Distributed
using StatsBase
using Distances

import ..Juscan: _get_obs_rep

# We provide graph-based and kmedoids methods for clustering.
# We also support a automatic selection of k for kmedoids.

"""
    clustering!(data::Muon.AnnData; kwargs...)

Perform clustering on a `Muon.AnnData` object and store the results in place.

This function supports modularity-based graph clustering (`"mc"`) and K-means clustering (`"km"`). By default, modularity clustering is applied using a shared nearest neighbor (SNN) graph constructed from PCA-reduced data.

# Arguments
- `data::Muon.AnnData`: The annotated data matrix to cluster. The object will be modified in place.

## Keyword Arguments
- `method::AbstractString = "mc"`: Clustering method. Options are `"mc"` (modularity clustering) or `"km"` (K-means).
- `reduction::Union{AbstractString, Symbol} = :auto`: Dimensionality reduction method to use. Accepts `"pca"` or `"harmony"`. When `:auto`, it defaults to `"pca"` (support for `"harmony"` is planned).
- `use_pca::Union{AbstractString, Integer} = "pca_cut"`: Number of principal components to use, or the key in `obsm` specifying a PCA representation.
- `tree_K::Integer = 20`: Number of neighbors to use when building the SNN graph. Relevant only for `"mc"` clustering.
- `resolution::Union{Symbol, Real, AbstractRange} = :auto`: Resolution(s) for modularity optimization. `:auto` uses `0.2:0.1:2.0`.
- `cluster_K::Union{Nothing, Integer} = nothing`: Number of clusters for K-means. If `nothing`, it will be auto-determined.
- `cluster_K_max::Union{Nothing, Integer} = 30`: Maximum number of clusters to try for automatic K-means clustering.
- `dist::AbstractString = "Euclidean"`: Distance metric for K-means, e.g., `"Euclidean"`.
- `network::AbstractString = "SNN"`: Type of graph network to construct. Currently supports `"SNN"`.
- `random_starts_number::Integer = 10`: Number of random initializations for clustering (modularity clustering only).
- `iter_number::Integer = 10`: Maximum number of iterations for the clustering optimization.
- `prune::AbstractFloat = 1/15`: Pruning factor for the graph. Must be between 0 and 1.
- `seed::Integer = -1`: Random seed. Use a negative number to skip setting the seed.

# Returns
Nothing. The clustering results are stored in the `obs` field of the input `AnnData` object, typically under a key like `"clusters"`.

# Example
```julia
using Juscan

clustering!(adata; method="mc", use_pca="X_pca", resolution=1.0)
```

# Notes
- `"harmony"` support is planned but not yet implemented.
- The `"mc"` method builds a neighbor graph and performs community detection; `"km"` performs K-means clustering in reduced space.
- The `method` parameter only supports `"mc"` and `"km"`; invalid inputs will return a warning.

"""
function clustering!(
  data::Muon.AnnData;
  method::AbstractString="mc",
  reduction::Union{AbstractString, Symbol}=:auto,
  use_pca::Union{AbstractString, Integer}="pca_cut",
  tree_K::Integer=20,
  resolution::Union{Symbol, Real, AbstractRange}=:auto,
  cluster_K::Union{Nothing, Integer}=nothing,
  cluster_K_max::Union{Nothing, Integer}=30,
  dist::AbstractString="Euclidean",
  network::AbstractString="SNN",
  random_starts_number::Integer=10,
  iter_number::Integer=10,
  prune::AbstractFloat=1 / 15,
  seed::Integer=-1,
)

  # Support PCA or Harmony
  if reduction == :auto
    # TODO: Harmony
    # if "harmony" in keys(obj.meta)
    #   reduction = "harmony"
    # else
    reduction = "pca"
    # end
  end

  # Check prune
  if prune > 1 || prune < 0
    throw(DomainError(prune, "Must be 0 ~ 1."))
  end
  # Check method
  # Leiden package has some problems. It should be re-written in future!
  # For HDF5.jl reason, now only two characters sign is allowed.
  if !(method in ["mc", "km"])
    return "Nothing to do! 'method' is only support \"mc\"|\"km\"!"
  end

  # Check PCs
  if typeof(use_pca) <: AbstractString
    # use_pca = obj.meta[use_pca]
    # use_pca = _get_obs_rep(data, layer=use_pca)
    use_pca = data.uns[use_pca]
  end

  # Check resolution
  if resolution == :auto
    resolution = 0.2:0.1:2.0
  elseif typeof(resolution) <: Real
    resolution = resolution
  end

  # Do clustering
  if method == "mc"
    modularityClustering!(
      data;
      reduction=reduction,
      num_pca=use_pca,
      network=network,
      K=tree_K,
      resolution=resolution,
      prune=prune,
      random_starts_number=random_starts_number,
      iter_number=iter_number,
      seed=seed,
    )
  else
    kClustering!(
      data;
      reduction=reduction,
      num_pca=use_pca,
      dist_method=dist,
      K=cluster_K,
      max_K=cluster_K_max,
      seed=seed,
    )
  end
end

function kClustering!(
  data::Muon.AnnData;
  reduction::AbstractString="pca",
  num_pca::Integer=5,
  dist_method::AbstractString="CorrDist",
  K::Union{Nothing, Integer}=nothing,
  max_K::Union{Nothing, Integer}=50,
  seed::Integer=-1,
)
  if dist_method in [
    "Euclidean",
    "SqEuclidean",
    "Cityblock",
    "Chebyshev",
    "Jaccard",
    "CosineDist",
    "CorrDist",
    "RMSDeviation",
  ]
    if !isnothing(num_pca)
      dist =
        pairwise(getfield(Distances, Symbol(dist_method))(), data.obsm[reduction][:, 1:num_pca]')
    else
      dist = pairwise(getfield(Distances, Symbol(dist_method))(), data.obsm[reduction]')
    end
  end

  if isnothing(K) && isnothing(max_K)
    return "There must be a \"cluster_K\" or \"cluster_K_max\"!"
  end

  if !isnothing(K)
    Random.seed!(seed < 0 ? 1984 : seed)
    cell_assignments = kmedoids(dist, K).assignments
  else
    # Automatically choose resolution, based the silhouettes coefficients 
    # and elbow point by total costs
    #
    T = eltype(dist)
    clustering_vector = Vector{KmedoidsResult{T}}()
    resize!(clustering_vector, max_K - 1)
    Random.seed!(seed < 0 ? 1984 : seed)
    # for (i,k) in enumerate(2:max_K)
    #     clustering_vector[i] = kmedoids(dist,k)
    # end
    pmap(x -> clustering_vector[x[1]] = kmedoids(dist, x[2]), enumerate(2:max_K))

    # Check silhouettes coefficients
    sil_check = fill(Vector(undef, 4), length(clustering_vector))
    @inbounds for (i, c) in enumerate(clustering_vector)
      sh = silhouettes(c, dist)
      sh_avg = mean(sh)
      sh_std = std(sh)
      check_max = true
      check_min = true
      @inbounds for i in unique(c.assignments)
        idx = c.assignments .== i
        # Check all clusters have the max coefficient more than whole 
        # average
        check_max = check_max ? maximum(sh[idx]) <= sh_avg ? false : true : false
        # Check all clusters have no coefficient less than zero
        check_min = check_min ? minimum(sh[idx]) >= 0 ? false : true : false
        if !check_min && !check_max
          break
        end
      end
      sil_check[i] = Vector(undef, 4)
      sil_check[i][1] = sh_avg
      sil_check[i][2] = sh_std
      sil_check[i][3] = check_max
      sil_check[i][4] = check_min
    end
    # Order and find the best K
    idx = findall(x -> all(x[3:4]), sil_check)
    if isempty(idx)
      @warn "No suitable K found based on silhouette checks, falling back to default K = 5"
      silhouettes_K = 5
    else
      res = [idx .+ 1 [sil_check[i][1] for i in idx] [sil_check[i][2] for i in idx]]
      sh_coef_order = res[sortperm(res[:, 2]), 1]
      sh_std_order = res[sortperm(res[:, 3]; rev=true), 1]
      # Find the top order of average of silhouettes coefficients and bottom 
      # order of standard deviation of silhouettes coefficients
      silhouettes_K =
        [
          (i + findall(x -> x == sh_coef_order[i], sh_std_order)[1], sh_coef_order[i]) for
          i in eachindex(sh_coef_order)
        ] |> x -> x[findmax([i[1] for i in x])[2]][2] |> Int
    end

    # Find elbow of total costs
    tc_fit = glm(
      @formula(y ~ x + x^2 + x^3),
      DataFrames.DataFrame(x=2.0:max_K, y=[c.totalcost for c in clustering_vector]),
      Gamma(),
      InverseLink(),
    )
    tc_fit2 = lm(
      @formula(y ~ x),
      DataFrames.DataFrame(
        x=[2, max_K],
        y=[clustering_vector[1].totalcost, clustering_vector[end].totalcost],
      ),
    )
    totalcost_K = findmax(
      GLM.predict(tc_fit2, DataFrames.DataFrame(x=2:max_K)) |>
      x -> x - [c.totalcost for c in clustering_vector],
    )[2]

    K = round(Int, (silhouettes_K + totalcost_K) / 2)
    @info "Automatically choose K: $K"
    cell_assignments = clustering_vector[K - 1].assignments
  end

  SortAssignments!(cell_assignments)
  # Output
  data.obs[!, "clusters_" * dist_method * "_" * string(K)] = cell_assignments
  data.obs.clusters_latest = cell_assignments
end

function modularityClustering!(
  data::Muon.AnnData;
  modularity_fun::Integer=1,
  algorithm::Integer=1,
  reduction::AbstractString="pca",
  num_pca::Integer=5,
  network::AbstractString="SNN",
  K::Integer=20,
  resolution::Union{AbstractRange, Real}=0.2:0.1:1.5,
  prune::AbstractFloat=1 / 15,
  random_starts_number::Integer=10,
  iter_number::Integer=10,
  seed::Integer=-1,
)
  if !(modularity_fun in [1, 2])
    @error("modularity_fun must be 1 or 2")
  end
  if !(algorithm in [1, 2, 3, 4])
    @error("algorithm must be 1, 2, 3 or 4!")
  end
  if random_starts_number < 1
    @error("random_starts_number must be more than or equal to 1!")
  end
  if iter_number < 1
    @error("iter_number must be more than or equal to 1!")
  end
  if any(resolution .> 1) && modularity_fun == 2
    @error("resolution must be less than or equal to 1, when 'modularity_fun' is 2")
  end

  seed = seed < 0 ? 1984 : seed

  if !isnothing(num_pca)
    d = data.obsm[reduction][:, 1:num_pca]
  else
    d = data.obsm[reduction]
  end

  tree = KDTree(d')
  res = knn(tree, d', K, true)[1]

  if network == "KNN"
    # KNN clustering is faster than SNN clustering!

    # # From Seurat, faster but not adjacent/symmetric matrix!
    # # Clustering result may be not good!
    # i = vcat(res...)
    # knn_matrix = sparse((eachindex(i) . - 1) .÷ K .+ 1,i,1)

    # Adjacent matrix, but slower!
    # Clustering result may be better!
    g = SimpleWeightedGraph(size(d, 1))
    # for v in vertices(g)
    #     for k in 1:K
    #         # add_edge!(g,v,res[v][k],1)
    #         g.weights[v,res[v][k]] = 1
    #     end
    # end
    map(v -> map(k -> add_edge!(g, v, res[v][k], 1), 1:K), vertices(g))
    m = g.weights
  elseif network == "SNN"
    m = SNN(hcat(res...)', prune)
  else
    @error("network must be 'KNN' or 'SNN'!")
  end

  # Check empty network
  empty_network = true
  @inbounds for i in eachindex(m)
    if m[i] != 0
      empty_network = false
      break
    end
  end
  if empty_network
    @error("Matrix has no network inside!")
  end

  if typeof(resolution) <: Real
    cell_assignments = (
      resolution,
      ModClustering(
        m;
        resolution=resolution,
        rsn=random_starts_number,
        itn=iter_number,
        seed=seed,
        modularity_fun=modularity_fun,
        algorithm=algorithm,
      ).clusters,
    )
  else
    clustering_vector = Tuple{Float32, Vector{Int}}[]
    slice_time = 0
    @inbounds @fastmath for i in eachindex(resolution)
      push!(
        clustering_vector,
        (
          resolution[i],
          ModClustering(
            m;
            resolution=resolution[i],
            rsn=random_starts_number,
            itn=iter_number,
            seed=seed,
            modularity_fun=modularity_fun,
            algorithm=algorithm,
          ).clusters,
        ),
      )
      if i == 1
        continue
      end
      slice_number = Ref(0)
      # Automatically detect the resolution, so can NOT use @tturbo
      # @inbounds @fastmath @simd for j in 
      #         unique(clustering_vector[i - 1][2])
      #     idx = clustering_vector[i - 1][2] .== j
      #     counter = StatsBase.counts(clustering_vector[i][2][idx])
      #     percentage = [ i / sum(counter) for i in counter ]
      #     # # Frankly, I think the commited condition is reasonable!
      #     # if isnothing(findfirst(x -> x > 0.99,percentage))
      #     if !isnothing(findfirst(x -> x > 0.1,percentage))
      #         slice_number[] += 1
      #     end
      # end
      map(
        j -> res_detect!(i, j, clustering_vector, slice_number),
        unique(clustering_vector[i - 1][2]),
      )

      if slice_number[] >= 2 &&
         length(unique(clustering_vector[i][2])) == length(unique(clustering_vector[i - 1][2]))
        slice_time += 1
      end
      if slice_time == 2
        @info "Recommended resolution is $(resolution[i - 1])"
        break
      end
    end
    pos = length(clustering_vector) - 1
    cell_assignments = (clustering_vector[pos][1], clustering_vector[pos][2])
  end

  SortAssignments!(cell_assignments[2])
  # Output
  data.obs[!, "clusters_" * string(cell_assignments[1])] = cell_assignments[2]
  data.obs.clusters_latest = cell_assignments[2]
end

@inline function res_detect!(
  i::Integer,
  j::Integer,
  clustering_vector::Vector,
  slice_number::Base.RefValue,
)
  idx = clustering_vector[i - 1][2] .== j
  counter = StatsBase.counts(clustering_vector[i][2][idx])
  percentage = [i / sum(counter) for i in counter]
  # # Frankly, I think the commited condition is reasonable!
  # if isnothing(findfirst(x -> x > 0.99,percentage))
  if !isnothing(findfirst(x -> x > 0.1, percentage))
    slice_number[] += 1
  end
end

function SortAssignments!(ca::AbstractVector)
  x = Int[]
  m = countmap(ca)
  @inbounds @fastmath @simd for _ in 2:length(m)
    mm = findmax(m)[2]
    push!(x, mm)
    pop!(m, mm)
  end
  idx = [ca .== v for v in x]
  @inbounds for (i, id) in enumerate(idx)
    ca[id] .= i
  end
end
