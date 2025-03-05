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
  @testset "quality control" begin
    @info "test quality control metrics"
    include("preprocessing/qc.jl")
  end

  @testset "filter cells" begin
    @info "test filter cells"
    include("preprocessing/filter.jl")
  end

  @testset "normalization" begin
    @info "test normalization"
    include("preprocessing/normalization.jl")
  end
end

@testset "preprocesing Module" begin
  @testset "highly variable genes" begin
    @info "test highly variable genes"
    include("tools/hvg.jl")
  end
end
