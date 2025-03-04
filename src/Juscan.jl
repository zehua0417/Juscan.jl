#=========================================
filename: Juscan
author: Lihuax
date: 2025-03-04 16:12
description: Main module for Juscan
=========================================#

module Juscan

using DataFrames
using HDF5
using Muon
using Muon: AnnData
using SparseArrays, LinearAlgebra

include("anndata.jl")
export AnnData,
  get_celltypes, subset_adata, subset_adata!, insert_obs!, insert_var!, read_h5ad, write_h5ad

include("preprocessing/pp.jl")
using .Pp
export calculate_qc_metrics,
  calculate_qc_metrics!, describe_obs, describe_obs!, describe_var, describe_var!

end # module Juscan
