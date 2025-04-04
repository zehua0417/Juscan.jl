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
  pca!(adata; layer="log_transformed", n_pcs=2, key_added="PCA", verbose=false)
  @test haskey(adata.obsm, "PCA")
  @test size(adata.obsm["PCA"]) == (size(adata.X, 1), 2)
end

@testset "test umap!" begin
  X_large = rand(50, 20)
  adata = AnnData(
    X=X_large,
    obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
    var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
    layers=Dict("log_transformed" => log.(X_large)),
  )
  umap!(
    adata;
    use_pca=nothing,
    layer="log_transformed",
    key_added="umap",
    n_neighbors=15,
    verbose=false,
  )
  @test haskey(adata.obsm, "umap")
  @test size(adata.obsm["umap"], 1) == 50
end
