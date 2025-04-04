@testset "Juscan Clustering Tests" begin
  X_large = rand(50, 20)
  adata = AnnData(
    X=X_large,
    obs=DataFrame(cell_type=["cell$(i)" for i in 1:size(X_large, 1)]),
    var=DataFrame(gene_name=["gene$(j)" for j in 1:size(X_large, 2)]),
    layers=Dict("log_transformed" => log.(X_large)),
  )

  adata.layers["counts"] = adata.X
  pca!(adata; layer="counts", n_pcs=7, key_added="pca", verbose=false)

  # 执行聚类
  clustering!(
    adata;
    # method="km",
    reduction="pca",
    use_pca=7,
    tree_K=5,
    resolution=0.5,
    seed=42,
    prune=0.1,
    random_starts_number=5,
    iter_number=5,
  )

  # 检查是否输出了聚类结果
  @test "clusters_0.5" in names(adata.obs)
  @test !any(ismissing, adata.obs[!, "clusters_0.5"])

  # 检查最新聚类结果是否同步到了 `clusters_latest`
  @test adata.obs.clusters_latest == adata.obs[!, "clusters_0.5"]
end
