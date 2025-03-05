using Juscan
using Juscan.Pp
using Test
using DataFrames
using SparseArrays

@testset "Basic AnnData operations" begin
  @info "Basic AnnData operations"
  include("anndata.jl")
end

@testset "preprocesing Module quality control" begin
  @info "preprocesing Module quality control metrics"
  include("preprocessing/qc.jl")
end

@testset "preprocesing Module filter cells" begin
  @info "preprocesing Module filter cells"
  include("preprocessing/filter.jl")
end

@testset "preprocesing Module normalization" begin
  @info "preprocesing Module normalization"
  include("preprocessing/normalization.jl")
end
