# @testset "test pca!" begin
#   # test if pca! works when log_transformed layer exists
#   adata = AnnData(
#     X=[2 0; 0 1; 10 2],
#     obs=DataFrame(cell_type=["a", "b", "c"]),
#     var=DataFrame(gene_name=["gene1", "gene2"]),
#     obsm=Dict("dim_red" => [1 2; 3 4; 5 6]),
#     varm=Dict("gene_vars" => [1 2 3; 3 4 5]),
#     varp=Dict("gene_corrs" => [1 2; 3 4]),
#     obsp=Dict("cell_corrs" => [1 2 3; 3 4 6; 5 6 7]),
#   )
#   pca!(adata; layer="counts", n_pcs=2, verbose=false)
#   @test haskey(adata.obsm, "PCA")
#   @test size(adata.obsm["PCA"]) == (size(adata.X, 1), 2)
# end
#
# X_large = rand(50, 20)
# adata1 = AnnData(
#   X=X_large,
#   obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
#   var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
#   layers=Dict("log_transformed" => log.(X_large)),
# )
# umap!(adata1; use_pca_init=true, n_pcs=10, n_neighbors=15, verbose=false)
# @test haskey(adata1.obsm, "umap")
# @test haskey(adata1.obsm, "knns")
# @test haskey(adata1.obsm, "knn_dists")
# @test haskey(adata1.obsp, "fuzzy_neighbor_graph")
# @test size(adata1.obsm["umap"], 1) == 2  # 2 维 embedding
#
# adata2 = AnnData(
#   X=X_large,
#   obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
#   var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
#   layers=Dict("log_transformed" => log.(X_large)),
# )
# umap!(adata2; use_pca_init=false, layer="log_transformed", n_neighbors=15, verbose=false)
# @test haskey(adata2.obsm, "umap")
# @test size(adata2.obsm["umap"], 1) == 2
#
# adata3 = AnnData(
#   X=X_large,
#   obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
#   var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
#   layers=Dict{String, Matrix{Float64}}(),  # 空的 layers
# )
# umap!(adata3; use_pca_init=false, layer="non_existing_layer", n_neighbors=15, verbose=false)
# @test haskey(adata3.layers, "normalized")
# @test haskey(adata3.layers, "log_transformed")
# @test haskey(adata3.obsm, "umap")
# @test haskey(adata3.obsp, "fuzzy_neighbor_graph")

@testset "test pca!" begin
  # test if pca! works when log_transformed layer exists
  adata = AnnData(
    X=[2 0; 0 1; 10 2],
    obs=DataFrame(cell_type=["a", "b", "c"]),
    var=DataFrame(gene_name=["gene1", "gene2"]),
    obsm=Dict("dim_red" => [1 2; 3 4; 5 6]),
    varm=Dict("gene_vars" => [1 2 3; 3 4 5]),
    varp=Dict("gene_corrs" => [1 2; 3 4]),
    obsp=Dict("cell_corrs" => [1 2 3; 3 4 6; 5 6 7]),
  )
  pca!(adata; layer="log_transformed", n_pcs=2, verbose=false)
  @test haskey(adata.obsm, "PCA")
  @test size(adata.obsm["PCA"]) == (size(adata.X, 1), 2)
end

@testset "test umap!" begin
  X_large = rand(50, 20)
  adata1 = AnnData(
    X=X_large,
    obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
    var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
    layers=Dict("log_transformed" => log.(X_large)),
  )
  umap!(adata1; use_pca_init=true, n_pcs=10, n_neighbors=15, verbose=false)
  @test haskey(adata1.obsm, "umap")
  @test haskey(adata1.obsm, "knns")
  @test haskey(adata1.obsm, "knn_dists")
  @test haskey(adata1.obsp, "fuzzy_neighbor_graph")
  @test size(adata1.obsm["umap"], 1) == 50  # 2 维 embedding

  adata2 = AnnData(
    X=X_large,
    obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
    var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
    layers=Dict("log_transformed" => log.(X_large)),
  )
  umap!(adata2; use_pca_init=false, layer="log_transformed", n_neighbors=15, verbose=false)
  @test haskey(adata2.obsm, "umap")
  @test size(adata2.obsm["umap"], 1) == 50

  adata3 = AnnData(
    X=X_large,
    obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
    var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
    layers=Dict{String, Matrix{Float64}}(),  # 空的 layers
  )
  umap!(adata3; use_pca_init=false, layer="non_existing_layer", n_neighbors=15, verbose=false)
  @test haskey(adata3.layers, "normalized")
  @test haskey(adata3.layers, "log_transformed")
  @test haskey(adata3.obsm, "umap")
  @test haskey(adata3.obsp, "fuzzy_neighbor_graph")
end
