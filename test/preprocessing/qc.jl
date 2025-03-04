adata = AnnData(
  X=sparse([2 0; 0 1; 10 2]),
  obs=DataFrame(cell_type=["a", "b", "c"]),
  var=DataFrame(gene_name=["gene1", "gene2"]),
  obsm=Dict("dim_red" => [1 2; 3 4; 5 6]),
  varm=Dict("gene_vars" => [1 2 3; 3 4 5]),
  varp=Dict("gene_corrs" => [1 2; 3 4]),
  obsp=Dict("cell_corrs" => [1 2 3; 3 4 6; 5 6 7]),
)
orig_adata = deepcopy(adata)

# test describe_obs
obs_metrics = Juscan.describe_obs(adata, percent_top=[2])
@test size(obs_metrics, 1) == 3
@test "n_genes_by_counts" in names(obs_metrics)
@test "total_counts" in names(obs_metrics)

# test describe_var
var_metrics = Juscan.describe_var(adata)
@test size(var_metrics, 1) == 2
@test "n_cells_by_counts" in names(var_metrics)
@test "total_counts" in names(var_metrics)

# Test calculate_qc_metrics
(obs_metrics, var_metrics) = calculate_qc_metrics(adata, percent_top=[2])

# Test size of obs_metrics
@test size(obs_metrics, 1) == 3
@test "n_genes_by_counts" in names(obs_metrics)
@test "total_counts" in names(obs_metrics)

# Test size of var_metrics
@test size(var_metrics, 1) == 2
@test "n_cells_by_counts" in names(var_metrics)
@test "total_counts" in names(var_metrics)
