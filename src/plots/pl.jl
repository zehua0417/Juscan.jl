module Pl

using ..DataFrames
using ..Random, ..Distributed
using ..SparseArrays, ..LinearAlgebra
using ..Muon
using Plots, GLMakie, CairoMakie
using Colors, ColorSchemes, ColorVectorSpace

include("plots.jl")
include("colors.jl")

export violin, scatter, hvg_scatter, plot_variance_ratio, plot_umap
export expand_palette, get_continuous_colormap

end  # module Pl
