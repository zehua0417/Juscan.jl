@testset "Highly Variable Genes Tests" begin
  X = rand(1:100, 100, 200)  # 100 个细胞, 200 个基因的计数矩阵
  obs = DataFrame(batch=rand(["A", "B"], 100))  # 批次信息
  var = DataFrame(gene=string.(1:200))  # 基因信息
  adata = AnnData(X=X, obs=obs, var=var)

  # run hvg calculation
  hvg_info = highly_variable_genes(adata, n_top_genes=50, batch_key="batch")

  @test size(hvg_info, 1) == 200  # 确保所有基因都有计算结果
  @test :highly_variable in propertynames(hvg_info)  # 确保 HVG 信息存在
  @test count(hvg_info.highly_variable) == 50  # 确保选出的高变基因数量正确

  highly_variable_genes!(adata, n_top_genes=50, batch_key="batch")
  @test :highly_variable in propertynames(adata.var)  # 确保 inplace 计算更新 var
  @test count(adata.var.highly_variable) == 50  # 确保选出的高变基因数量正确

  subset_to_hvg!(adata, n_top_genes=50, batch_key="batch")
  @test size(adata.X, 2) == 50  # 确保矩阵被正确子集化
end
