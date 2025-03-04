using Juscan
using Juscan.Pp
using Test
using DataFrames
using SparseArrays

@testset "Basic AnnData operations" begin
  @info "Basic AnnData operations"
  include("anndata.jl")
end

@testset "preprocesing Module" begin
  @info "preprocesing Module quality control metrics"
  include("preprocessing/qc.jl")
end
