#=========================================
filename: Juscan
author: Lihuax
date: 2025-03-04 16:12
description: Main module for Juscan
=========================================#

module Juscan

using DataFrames
using HDF5
using Random
using Muon
using Muon: AnnData
using SparseArrays, LinearAlgebra
using Distributed

include("anndata.jl")
export AnnData,
  get_celltypes, subset_adata, subset_adata!, insert_obs!, insert_var!, _get_obs_rep, _set_obs_rep!
export readh5ad, writeh5ad

include("preprocessing/pp.jl")
using .Pp
export calculate_qc_metrics,
  calculate_qc_metrics!,
  describe_obs,
  describe_obs!,
  describe_var,
  describe_var!,
  filter_cells,
  filter_genes,
  filter_genes!,
  filter_cells!,
  normalize_total!,
  normalize_total

include("tools/tl.jl")
using .Tl
export highly_variable_genes!,
  highly_variable_genes, subset_to_hvg!, pca!, umap!, log_transform!, logp1_transform!, clustering!

end # module Juscan
