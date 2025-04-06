using Juscan
using Juscan.Pp
using Test
using DataFrames
using SparseArrays
using Muon

function runtests()
  all_results = Test.DefaultTestSet[]  # 存放每组测试结果

  function run_named_test(name::AbstractString, filepath::AbstractString)
    @info "Testing $name..."
    return @testset "$name" begin
      include(filepath)
    end
  end

  push!(all_results, @testset "Basic AnnData operations" begin
    run_named_test("AnnData", "anndata.jl")
  end)

  push!(all_results, @testset "preprocessing Module" begin
    run_named_test("quality control", "preprocessing/qc.jl")
    run_named_test("filter cells", "preprocessing/filter.jl")
    run_named_test("normalization", "preprocessing/normalization.jl")
  end)

  push!(all_results, @testset "tools Module" begin
    run_named_test("highly variable genes", "tools/hvg.jl")
    run_named_test("pca", "tools/pca.jl")
    run_named_test("cluster", "tools/cluster.jl")
  end)

  push!(all_results, @testset "plots Module" begin
    run_named_test("colors", "plots/colors.jl")
  end)

  for result in all_results
    println()
    Test.print_test_results(result)
  end
end

runtests()
