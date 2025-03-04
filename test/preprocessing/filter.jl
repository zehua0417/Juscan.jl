X = adjoint(sparse([1 0 3; 0 5 0; 2 0 4]))  # 3 genes x 3 cells

cell_subset, num_counts = filter_cells(X; min_counts=4)
@test cell_subset == [false, true, true]  # 只保留 cell 2 和 3
@test num_counts == [3, 5, 7]

cell_subset, num_genes = filter_cells(X; min_genes=2)
@test cell_subset == [true, false, true]  # 只保留 cell 1 和 3

cell_subset, _ = filter_cells(X; max_counts=5)
@test cell_subset == [true, true, false]  # cell 3 被移除

adata = Muon.AnnData(X=X)
adata_filtered = filter_cells(adata; min_counts=4, copy=true)
@test size(adata_filtered.X) == (2, 3)  # 过滤后应只剩 2 个细胞

gene_subset, num_counts = filter_genes(X; min_counts=3)
@test gene_subset == [true, true, true]  # 只保留 gene 1 和 3

gene_subset, num_cells = filter_genes(X; min_cells=2)
@test gene_subset == [true, false, true]  # 只保留 gene 1 和 3

gene_subset, _ = filter_genes(X; max_counts=5)
@test gene_subset == [true, true, false]  # gene 3 被移除

adata = Muon.AnnData(X=X)
adata_filtered = filter_genes(adata; min_cells=2, copy=true)
@test size(adata_filtered.X) == (3, 2)  # 过滤后应只剩 2 个基因
