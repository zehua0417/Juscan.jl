var documenterSearchIndex = {"docs":
[{"location":"api/preprocessing/#Preprocessing-Modules","page":"Preprocessing","title":"Preprocessing Modules","text":"","category":"section"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"Preprocessing is a crucial step in single-cell and multi-omics data analysis. This module provides functions for quality control, filtering, and normalization to ensure that datasets are clean and ready for downstream analyses.","category":"page"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"Quality control (QC) metrics help identify low-quality cells and genes, such as those with extremely high or low counts, high dropout rates, or disproportionate expression levels. By computing these metrics, users can apply appropriate thresholds to retain only reliable data.","category":"page"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"Filtering enables users to remove unwanted cells or genes based on QC criteria, ensuring that low-quality features do not affect downstream analyses.","category":"page"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"Normalization methods adjust for differences in sequencing depth and technical variation, allowing meaningful comparisons across cells and conditions.","category":"page"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"These preprocessing functions ensure that data is well-structured and suitable for clustering, dimensionality reduction, and other analytical tasks.","category":"page"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"Pages = [\"preprocessing.md\"]","category":"page"},{"location":"api/preprocessing/#quality-control-metrics","page":"Preprocessing","title":"quality control metrics","text":"","category":"section"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"calculate_qc_metrics\ncalculate_qc_metrics!","category":"page"},{"location":"api/preprocessing/#Juscan.Pp.calculate_qc_metrics","page":"Preprocessing","title":"Juscan.Pp.calculate_qc_metrics","text":"calculate_qc_metrics!(adata; expr_type=\"counts\", var_type=\"genes\", qc_vars=String[], \n                      percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, \n                      use_log1p=true, parallel=nothing) -> Tuple{DataFrame, DataFrame}\n\nCalculate quality control metrics and return.\n\nArguments\n\nadata: Muon.AnnData   Annotated data matrix containing single-cell or multi-omics data.\nexpr_type: String (default: \"counts\")   The type of expression data to use, such as raw counts or normalized values.\nvar_type: String (default: \"genes\")   Specifies the variable type, e.g., \"genes\" for gene expression.\nqc_vars: Union{Vector{String}, String} (default: String[])   A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.\npercent_top: Union{Vector{Int}, Nothing} (default: [50, 100, 200, 500])   A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.\nlayer: Union{String, Nothing} (default: nothing)   Specifies the layer in adata to use. If nothing, the default matrix is used.\nuse_raw: Bool (default: false)   Whether to use the raw counts matrix (adata.raw) instead of the processed data.\nuse_log1p: Bool (default: true)   Whether to apply log1p transformation to the data before computing quality control metrics.\nparallel: Union{Bool, Nothing} (default: nothing)   (Deprecated) This argument is ignored but retained for backward compatibility.\n\nReturns\n\nA tuple of two DataFrames.DataFrame:  \nThe first DataFrame contains per-cell quality control metrics.  \nThe second DataFrame contains per-gene quality control metrics.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#Juscan.Pp.calculate_qc_metrics!","page":"Preprocessing","title":"Juscan.Pp.calculate_qc_metrics!","text":"calculate_qc_metrics!(adata; expr_type=\"counts\", var_type=\"genes\", qc_vars=String[], \n                      percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, \n                      use_log1p=true, parallel=nothing) -> Nothing\n\nCalculate and store quality control (QC) metrics directly in the obs and var attributes of adata.\n\nArguments\n\nadata: Muon.AnnData   Annotated data matrix containing single-cell or multi-omics data.   The QC metrics will be stored in adata.obs (per-cell) and adata.var (per-gene).\nexpr_type: String (default: \"counts\")   The type of expression data to use, such as raw counts or normalized values.\nvar_type: String (default: \"genes\")   Specifies the variable type, e.g., \"genes\" for gene expression.\nqc_vars: Union{Vector{String}, String} (default: String[])   A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.\npercent_top: Union{Vector{Int}, Nothing} (default: [50, 100, 200, 500])   A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.\nlayer: Union{String, Nothing} (default: nothing)   Specifies the layer in adata to use. If nothing, the default matrix is used.\nuse_raw: Bool (default: false)   Whether to use the raw counts matrix (adata.raw) instead of the processed data.\nuse_log1p: Bool (default: true)   Whether to apply log1p transformation to the data before computing quality control metrics.\nparallel: Union{Bool, Nothing} (default: nothing)   (Deprecated) This argument is ignored but retained for backward compatibility.\n\nReturns\n\nNothing:   This function modifies adata in place by adding QC metrics to adata.obs and adata.var.\n\nNotes\n\nUnlike calculate_qc_metrics, which returns QC metrics as DataFrame objects,   this function directly updates adata without returning any values.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#sub-functions","page":"Preprocessing","title":"sub functions","text":"","category":"section"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"describe_obs\ndescribe_obs!\ndescribe_var\ndescribe_var!","category":"page"},{"location":"api/preprocessing/#Juscan.Pp.describe_obs","page":"Preprocessing","title":"Juscan.Pp.describe_obs","text":"describe_obs(adata; expr_type=\"counts\", var_type=\"genes\", qc_vars=String[], \n             percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, \n             use_log1p=true, X=nothing, parallel=nothing) -> DataFrames.DataFrame\n\nCompute per-cell quality control (QC) metrics from the expression matrix and return them as a DataFrame.\n\nArguments\n\nadata: Muon.AnnData   Annotated data matrix containing single-cell or multi-omics data.\nexpr_type: String (default: \"counts\")   The type of expression data to use, such as raw counts or normalized values.\nvar_type: String (default: \"genes\")   Specifies the variable type, e.g., \"genes\" for gene expression.\nqc_vars: Union{Vector{String}, String} (default: String[])   A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.\npercent_top: Union{Vector{Int}, Nothing} (default: [50, 100, 200, 500])   A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.\nlayer: Union{String, Nothing} (default: nothing)   Specifies the layer in adata to use. If nothing, the default matrix is used.\nuse_raw: Bool (default: false)   Whether to use the raw counts matrix (adata.raw) instead of the processed data.\nuse_log1p: Bool (default: true)   Whether to apply log1p transformation to the data before computing quality control metrics.\nX: Union{Nothing, SparseMatrixCSC} (default: nothing)   The expression matrix to use. If nothing, the matrix is determined by _choose_mtx_rep.\nparallel: Union{Bool, Nothing} (default: nothing)   (Deprecated) This argument is ignored but retained for backward compatibility.\n\nReturns\n\nDataFrames.DataFrame:   A DataFrame containing per-cell QC metrics, including:\n\"n_$(var_type)_by_$(expr_type)\": Number of detected features per cell.\n\"log1p_n_$(var_type)_by_$(expr_type)\": Log-transformed number of features per cell.\n\"total_$(expr_type)\": Total expression count per cell.\n\"log1p_total_$(expr_type)\": Log-transformed total expression count per cell.\n\"pct_$(expr_type)_in_top_N_$(var_type)\": Percentage of expression from the top N features.\n\"total_$(expr_type)_$(qc_var)\": Total expression for each QC variable per cell.\n\"log1p_total_$(expr_type)_$(qc_var)\": Log-transformed total expression for each QC variable per cell.\n\"pct_$(expr_type)_$(qc_var)\": Percentage of total expression contributed by each QC variable.\n\nNotes\n\nThis function computes per-cell QC metrics based on expression data.   For per-gene metrics, see describe_var.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#Juscan.Pp.describe_obs!","page":"Preprocessing","title":"Juscan.Pp.describe_obs!","text":"describe_obs!(adata; expr_type=\"counts\", var_type=\"genes\", qc_vars=String[], \n              percent_top=[50, 100, 200, 500], layer=nothing, use_raw=false, \n              use_log1p=true, X=nothing, parallel=nothing) -> Nothing\n\nCompute per-cell quality control (QC) metrics and store them in adata.obs.\n\nArguments\n\nadata: Muon.AnnData   Annotated data matrix containing single-cell or multi-omics data.\nexpr_type: String (default: \"counts\")   The type of expression data to use, such as raw counts or normalized values.\nvar_type: String (default: \"genes\")   Specifies the variable type, e.g., \"genes\" for gene expression.\nqc_vars: Union{Vector{String}, String} (default: String[])   A list of quality control variables to consider. If a single string is provided, it will be converted to a vector.\npercent_top: Union{Vector{Int}, Nothing} (default: [50, 100, 200, 500])   A list of top feature (e.g., gene) counts to compute the percentage of total counts for each cell.\nlayer: Union{String, Nothing} (default: nothing)   Specifies the layer in adata to use. If nothing, the default matrix is used.\nuse_raw: Bool (default: false)   Whether to use the raw counts matrix (adata.raw) instead of the processed data.\nuse_log1p: Bool (default: true)   Whether to apply log1p transformation to the data before computing quality control metrics.\nX: Union{Nothing, SparseMatrixCSC} (default: nothing)   The expression matrix to use. If nothing, the matrix is determined by _choose_mtx_rep.\nparallel: Union{Bool, Nothing} (default: nothing)   (Deprecated) This argument is ignored but retained for backward compatibility.\n\nReturns\n\nNothing:   This function modifies adata in place by adding per-cell QC metrics to adata.obs.\n\nNotes\n\nThis function is similar to describe_obs but modifies adata directly.   It computes QC metrics using describe_obs and stores the results in adata.obs.  \nFor per-gene QC metrics, see describe_var!.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#Juscan.Pp.describe_var","page":"Preprocessing","title":"Juscan.Pp.describe_var","text":"describe_var(adata; expr_type=\"counts\", var_type=\"genes\", layer=nothing, \n             use_raw=false, use_log1p=true, X=nothing) -> DataFrame\n\nCompute per-gene quality control (QC) metrics and return as a DataFrame.\n\nArguments\n\nadata: Muon.AnnData   Annotated data matrix containing single-cell or multi-omics data.\nexpr_type: String (default: \"counts\")   The type of expression data to use, such as raw counts or normalized values.\nvar_type: String (default: \"genes\")   Specifies the variable type, e.g., \"genes\" for gene expression.\nlayer: Union{String, Nothing} (default: nothing)   Specifies the layer in adata to use. If nothing, the default matrix is used.\nuse_raw: Bool (default: false)   Whether to use the raw counts matrix (adata.raw) instead of the processed data.\nuse_log1p: Bool (default: true)   Whether to apply log1p transformation to the data before computing quality control metrics.\nX: Union{Nothing, SparseMatrixCSC} (default: nothing)   The expression matrix to use. If nothing, the matrix is determined by _choose_mtx_rep.\n\nReturns\n\nDataFrames.DataFrame:   A DataFrame where each row represents a gene (or variable), containing the following QC metrics:\n\"n_cells_by_<expr_type>\": Number of cells in which each gene is detected.\n\"mean_<expr_type>\": Mean expression level of each gene.\n\"log1p_mean_<expr_type>\": Log-transformed mean expression level.\n\"pct_dropout_by_<expr_type>\": Percentage of cells in which each gene is not detected.\n\"total_<expr_type>\": Total expression level of each gene.\n\"log1p_total_<expr_type>\": Log-transformed total expression level.\n\nNotes\n\nThis function calculates quality control metrics per gene (or other variables).  \nTo store these metrics directly in adata.var, use describe_var!.  \nTo compute per-cell QC metrics, see describe_obs.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#Juscan.Pp.describe_var!","page":"Preprocessing","title":"Juscan.Pp.describe_var!","text":"describe_var!(adata; expr_type=\"counts\", var_type=\"genes\", layer=nothing, \n              use_raw=false, use_log1p=true, X=nothing) -> Nothing\n\nCompute per-gene quality control (QC) metrics and store them in adata.var.\n\nArguments\n\nadata: Muon.AnnData   Annotated data matrix containing single-cell or multi-omics data.\nexpr_type: String (default: \"counts\")   The type of expression data to use, such as raw counts or normalized values.\nvar_type: String (default: \"genes\")   Specifies the variable type, e.g., \"genes\" for gene expression.\nlayer: Union{String, Nothing} (default: nothing)   Specifies the layer in adata to use. If nothing, the default matrix is used.\nuse_raw: Bool (default: false)   Whether to use the raw counts matrix (adata.raw) instead of the processed data.\nuse_log1p: Bool (default: true)   Whether to apply log1p transformation to the data before computing quality control metrics.\nX: Union{Nothing, SparseMatrixCSC} (default: nothing)   The expression matrix to use. If nothing, the matrix is determined by _choose_mtx_rep.\n\nReturns\n\nNothing:   This function modifies adata in place by adding per-gene QC metrics to adata.var.\n\nNotes\n\nThis function is similar to describe_var but modifies adata directly.   It computes QC metrics using describe_var and stores the results in adata.var.  \nFor per-cell QC metrics, see describe_obs!.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#filter-cells-and-genes","page":"Preprocessing","title":"filter cells and genes","text":"","category":"section"},{"location":"api/preprocessing/#filter-cells","page":"Preprocessing","title":"filter cells","text":"","category":"section"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"filter_cells\nfilter_cells!","category":"page"},{"location":"api/preprocessing/#Juscan.Pp.filter_cells","page":"Preprocessing","title":"Juscan.Pp.filter_cells","text":"filter_cells(data::Muon.AnnData; min_counts=nothing, min_genes=nothing, max_counts=nothing, max_genes=nothing, copy=false) -> Union{Tuple{BitVector, Vector}, Muon.AnnData}\n\nFilters cells based on the given threshold criteria and returns either a subset mask and count vector or a new filtered Muon.AnnData object.\n\nArguments\n\ndata::Muon.AnnData: The input single-cell data.\nmin_counts::Union{Int, Nothing}: Minimum total counts per cell.\nmin_genes::Union{Int, Nothing}: Minimum number of genes expressed per cell.\nmax_counts::Union{Int, Nothing}: Maximum total counts per cell.\nmax_genes::Union{Int, Nothing}: Maximum number of genes expressed per cell.\ncopy::Bool: If true, returns a filtered copy; otherwise, returns a mask and count vector.\n\nReturns\n\nIf copy == false, returns a tuple (cell_subset::BitVector, number_per_cell::Vector).\nIf copy == true, returns a filtered Muon.AnnData object.\n\n\n\n\n\nfilter_cells(data::AbstractMatrix; min_counts=nothing, min_genes=nothing, max_counts=nothing, max_genes=nothing) -> Tuple{BitVector, Vector}\n\nFilters cells from a count matrix based on the given threshold criteria.\n\nArguments\n\ndata::AbstractMatrix: Gene expression count matrix with cells as rows.\nmin_counts::Union{Int, Nothing}: Minimum total counts per cell.\nmin_genes::Union{Int, Nothing}: Minimum number of genes expressed per cell.\nmax_counts::Union{Int, Nothing}: Maximum total counts per cell.\nmax_genes::Union{Int, Nothing}: Maximum number of genes expressed per cell.\n\nReturns\n\ncell_subset::BitVector: A mask indicating cells that pass the filter.\nnumber_per_cell::Vector: A vector of counts or expressed gene numbers per cell.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#Juscan.Pp.filter_cells!","page":"Preprocessing","title":"Juscan.Pp.filter_cells!","text":"filter_cells!(data::Muon.AnnData; min_counts=nothing, min_genes=nothing, max_counts=nothing, max_genes=nothing) -> Nothing\n\nFilters cells in-place in a Muon.AnnData object.\n\nArguments\n\ndata::Muon.AnnData: The input single-cell data.\nmin_counts::Union{Int, Nothing}: Minimum total counts per cell.\nmin_genes::Union{Int, Nothing}: Minimum number of genes expressed per cell.\nmax_counts::Union{Int, Nothing}: Maximum total counts per cell.\nmax_genes::Union{Int, Nothing}: Maximum number of genes expressed per cell.\n\nReturns\n\nNothing. The filtering is applied in-place.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#filter-genes","page":"Preprocessing","title":"filter genes","text":"","category":"section"},{"location":"api/preprocessing/","page":"Preprocessing","title":"Preprocessing","text":"filter_genes\nfilter_genes!","category":"page"},{"location":"api/preprocessing/#Juscan.Pp.filter_genes","page":"Preprocessing","title":"Juscan.Pp.filter_genes","text":"filter_genes(data::Muon.AnnData; min_counts=nothing, min_cells=nothing, max_counts=nothing, max_cells=nothing, copy=false) -> Union{Tuple{BitVector, Vector}, Muon.AnnData}\n\nFilters genes based on the given threshold criteria and returns either a subset mask and count vector or a new filtered Muon.AnnData object.\n\nArguments\n\ndata::Muon.AnnData: The input single-cell data.\nmin_counts::Union{Int, Nothing}: Minimum total counts per gene.\nmin_cells::Union{Int, Nothing}: Minimum number of cells expressing the gene.\nmax_counts::Union{Int, Nothing}: Maximum total counts per gene.\nmax_cells::Union{Int, Nothing}: Maximum number of cells expressing the gene.\ncopy::Bool: If true, returns a filtered copy; otherwise, returns a mask and count vector.\n\nReturns\n\nIf copy == false, returns a tuple (gene_subset::BitVector, number_per_gene::Vector).\nIf copy == true, returns a filtered Muon.AnnData object.\n\n\n\n\n\nfilter_genes(data::AbstractMatrix; min_counts=nothing, min_cells=nothing, max_counts=nothing, max_cells=nothing) -> Tuple{BitVector, Vector}\n\nFilters genes from a count matrix based on the given threshold criteria.\n\nArguments\n\ndata::AbstractMatrix: Gene expression count matrix with genes as columns.\nmin_counts::Union{Int, Nothing}: Minimum total counts per gene.\nmin_cells::Union{Int, Nothing}: Minimum number of cells expressing the gene.\nmax_counts::Union{Int, Nothing}: Maximum total counts per gene.\nmax_cells::Union{Int, Nothing}: Maximum number of cells expressing the gene.\n\nReturns\n\ngene_subset::BitVector: A mask indicating genes that pass the filter.\nnumber_per_gene::Vector: A vector of counts or expressed cell numbers per gene.\n\n\n\n\n\n","category":"function"},{"location":"api/preprocessing/#Juscan.Pp.filter_genes!","page":"Preprocessing","title":"Juscan.Pp.filter_genes!","text":"filter_genes!(data::Muon.AnnData; min_counts=nothing, min_cells=nothing, max_counts=nothing, max_cells=nothing) -> Nothing\n\nFilters genes in-place in a Muon.AnnData object.\n\nArguments\n\ndata::Muon.AnnData: The input single-cell data.\nmin_counts::Union{Int, Nothing}: Minimum total counts per gene.\nmin_cells::Union{Int, Nothing}: Minimum number of cells expressing the gene.\nmax_counts::Union{Int, Nothing}: Maximum total counts per gene.\nmax_cells::Union{Int, Nothing}: Maximum number of cells expressing the gene.\n\nReturns\n\nNothing. The filtering is applied in-place.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/quickstart/#todo","page":"Quick Start","title":"todo","text":"","category":"section"},{"location":"#Juscan","page":"Home","title":"Juscan","text":"","category":"section"},{"location":"api/anndata/#Anndata-utils","page":"Anndata Utils","title":"Anndata utils","text":"","category":"section"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"The AnnData struct is imported from Muon.jl. The package provides read and write functions for .h5ad and .h5mu files, the typical H5-based format for storing Python anndata objects. The AnnData object stores datasets together with metadata, such as information on the variables (genes in scRNA-seq data) and observations (cells), as well as different kinds of annotations and transformations of the original count matrix, such as PCA or UMAP embeddings, or graphs of observations or variables.","category":"page"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"For details on the Julia implementation in Muon.jl, see the documentation.","category":"page"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"For more details on the original Python implementation of the anndata object, see the documentation and preprint.","category":"page"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"(Image: anndata)","category":"page"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"Pages = [\"anndata.md\"]","category":"page"},{"location":"api/anndata/#celltypes","page":"Anndata Utils","title":"celltypes","text":"","category":"section"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"get_celltypes","category":"page"},{"location":"api/anndata/#Juscan.get_celltypes","page":"Anndata Utils","title":"Juscan.get_celltypes","text":"get_celltypes(a::AnnData)\n\nTries to infer the cell types of cells in an AnnData object. \n\nReturns a vector of cell type names if the cell types are stored in  adata.obs[\"cell_type\"], adata.obs[\"celltype\"], adata.obs[\"celltypes\"], or adata.obs[\"cell_types\"].  Otherwise, returns nothing.\n\n\n\n\n\n","category":"function"},{"location":"api/anndata/#subset-adata","page":"Anndata Utils","title":"subset adata","text":"","category":"section"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"subset_adata\nsubset_adata!","category":"page"},{"location":"api/anndata/#Juscan.subset_adata","page":"Anndata Utils","title":"Juscan.subset_adata","text":"subset_adata(adata::AnnData, subset_inds::Tuple, dims::Symbol=:both)\n\nSubset an AnnData object by indices passed as a tuple of vectors of integers, UnitRanges, or vectors of Booleans. If dims is set to :both, the first element of subset_inds is used to subset cells and the second element is used to subset genes.\n\nArguments\n\nadata: AnnData object to subset\nsubset_inds: tuple of vectors of integers, UnitRanges, or vectors of Booleans\ndims: dimension to subset, either :cells, :genes, or :both\n\nReturns\n\na copy of the AnnData object with the subsetted data\n\n\n\n\n\nsubset_adata(adata::AnnData, subset_inds::Union{Int, Vector{Int}, UnitRange, Vector{Bool}}, dims::Symbol)\n\nSubset an AnnData object by indices passed as a vector of integers or booleans or as a UnitRange.  The dims argument can be set to either :cells or :genes to specify which dimension to subset. \n\nArguments\n\nadata: AnnData object to subset\nsubset_inds: vector of integers or booleans or UnitRange\ndims: dimension to subset, either :cells or :genes\n\nReturns\n\na copy of the AnnData object with the subsetted data\n\n\n\n\n\n","category":"function"},{"location":"api/anndata/#Juscan.subset_adata!","page":"Anndata Utils","title":"Juscan.subset_adata!","text":"subset_adata!(adata::AnnData, subset_inds, dims::Symbol)\n\nIn-place version of subset_adata, see ?subset_adata for more details. \n\nFor subset_inds, either a tuple of vectors or ranges can be passed with dims set to :both, for subsetting  both cells and genes, or a single vector or range can be passed with dims set to either :cells or :genes.\n\nArguments\n\nadata: AnnData object to subset\nsubset_inds: tuple of vectors of integers, UnitRanges, or vectors of Booleans or vector of integers or booleans or UnitRange\ndims: dimension to subset, either :cells, :genes, or :both\n\nReturns\n\nthe AnnData object with the subsetted data\n\n\n\n\n\n","category":"function"},{"location":"api/anndata/#insert","page":"Anndata Utils","title":"insert","text":"","category":"section"},{"location":"api/anndata/","page":"Anndata Utils","title":"Anndata Utils","text":"insert_obs!\ninsert_var!","category":"page"},{"location":"api/anndata/#Juscan.insert_obs!","page":"Anndata Utils","title":"Juscan.insert_obs!","text":"insert_obs!(adata::AnnData, obs::DataFrame)\n\nInsert new columns into the obs DataFrame of an AnnData object.\n\nArguments\n\nadata: The AnnData object to modify.\nobs: A DataFrame containing the new columns to be inserted. The number of rows in obs must match the number of cells in adata.\n\nReturns\n\nnothing: Modifies the adata.obs DataFrame in place.\n\nNotes\n\nThis function assumes that the number of rows in obs matches the number of cells in adata. If they do not match, an error will be thrown.\nThis function modifies the adata.obs DataFrame directly. If you want to keep the original adata unchanged, make a copy before calling this function.\n\n\n\n\n\n","category":"function"},{"location":"api/anndata/#Juscan.insert_var!","page":"Anndata Utils","title":"Juscan.insert_var!","text":"insert_var!(adata::AnnData, var::DataFrame)\n\nInsert new columns into the var DataFrame of an AnnData object.\n\nArguments\n\nadata: The AnnData object to modify.\nvar: A DataFrame containing the new columns to be inserted. The number of rows in var must match the number of genes in adata.\n\nReturns\n\nnothing: Modifies the adata.var DataFrame in place.\n\nNotes\n\nThis function assumes that the number of rows in var matches the number of genes in adata. If they do not match, an error will be thrown.\nThis function modifies the adata.var DataFrame directly. If you want to keep the original adata unchanged, make a copy before calling this function.\n\n\n\n\n\n","category":"function"}]
}
