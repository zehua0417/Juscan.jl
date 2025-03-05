X = [1 2 3; 4 5 6; 7 8 9]  # 3x3 矩阵
counts = sum(X, dims=2)[:]
adata = Muon.AnnData(X=X)

@testset "normalize_total! function" begin
  target_sum = 1000.0
  normalize_total!(adata, target_sum=target_sum)
  norm_counts = sum(adata.X, dims=2)[:]
  @test all(norm_counts .≈ target_sum)
end

@testset "normalize_total function" begin
  adata_copy = Muon.AnnData(X=X)
  target_sum = 500.0
  norm_data = normalize_total(adata_copy, target_sum=target_sum)
  norm_X = norm_data["X"]
  norm_counts = sum(norm_X, dims=2)[:]
  @test all(norm_counts .≈ target_sum)
end
