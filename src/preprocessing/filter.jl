# ==========================================================
# Filename: filter.jl
# Author: Lihuax
# Date: 2025-01-13 14:11
# Description: This file contains functions for filtering cells and genes in high-dimensional biological data structures, including Muon.AnnData and sparse matrices. These functions implement customizable filtering based on various thresholds, enabling efficient preprocessing of single-cell RNA sequencing data.
# ==========================================================

import ..Juscan: subset_adata!

#==============
# filter cell #
==============#
"""
    filter_cells(
      data::Muon.AnnData;
      min_counts=nothing, min_genes=nothing,
      max_counts=nothing, max_genes=nothing,
      copy=false) -> Union{Tuple{BitVector, Vector}, Muon.AnnData}

Filters cells based on the given threshold criteria and returns either a subset mask and count vector or a new filtered `Muon.AnnData` object.

# Arguments
- `data::Muon.AnnData`: The input single-cell data.
- `min_counts::Union{Int, Nothing}`: Minimum total counts per cell.
- `min_genes::Union{Int, Nothing}`: Minimum number of genes expressed per cell.
- `max_counts::Union{Int, Nothing}`: Maximum total counts per cell.
- `max_genes::Union{Int, Nothing}`: Maximum number of genes expressed per cell.
- `copy::Bool`: If `true`, returns a filtered copy; otherwise, returns a mask and count vector.

# Returns
- If `copy == false`, returns a tuple `(cell_subset::BitVector, number_per_cell::Vector)`.
- If `copy == true`, returns a filtered `Muon.AnnData` object.
"""
function filter_cells(
  data::Muon.AnnData;
  min_counts::Union{Int, Nothing}=nothing,
  min_genes::Union{Int, Nothing}=nothing,
  max_counts::Union{Int, Nothing}=nothing,
  max_genes::Union{Int, Nothing}=nothing,
  copy::Bool=false,
)::Union{Tuple{BitVector, Vector}, Muon.AnnData}
  adata_copy = copy ? deepcopy(data) : data
  X = _choose_mtx_rep(adata_copy)

  cell_subset, number_per_cell = filter_cells(
    X,
    min_counts=min_counts,
    min_genes=min_genes,
    max_counts=max_counts,
    max_genes=max_genes,
  )
  if !copy
    return cell_subset, number_per_cell
  else
    if isnothing(min_genes) && isnothing(max_genes)
      adata_copy.obs[!, :n_counts] = number_per_cell
    else
      adata_copy.obs[!, :n_genes] = number_per_cell
    end
    # transform BitVector to something we can subset on
    if isa(cell_subset, BitVector)
      cell_subset = collect(cell_subset)
    end
    # filter adata in place
    subset_adata!(adata_copy, cell_subset, :cells)
    return adata_copy
  end
end

"""
    filter_cells(
      data::AbstractMatrix;
      min_counts=nothing, min_genes=nothing,
      max_counts=nothing, max_genes=nothing
    ) -> Tuple{BitVector, Vector}

Filters cells from a count matrix based on the given threshold criteria.

# Arguments
- `data::AbstractMatrix`: Gene expression count matrix with cells as rows.
- `min_counts::Union{Int, Nothing}`: Minimum total counts per cell.
- `min_genes::Union{Int, Nothing}`: Minimum number of genes expressed per cell.
- `max_counts::Union{Int, Nothing}`: Maximum total counts per cell.
- `max_genes::Union{Int, Nothing}`: Maximum number of genes expressed per cell.

# Returns
- `cell_subset::BitVector`: A mask indicating cells that pass the filter.
- `number_per_cell::Vector`: A vector of counts or expressed gene numbers per cell.
"""
function filter_cells(
  data::AbstractMatrix;
  min_counts::Union{Int, Nothing}=nothing,
  min_genes::Union{Int, Nothing}=nothing,
  max_counts::Union{Int, Nothing}=nothing,
  max_genes::Union{Int, Nothing}=nothing,
)::Tuple{BitVector, Vector}

  # Check that only one filtering option is provided
  options = [min_genes, min_counts, max_genes, max_counts]
  n_given_options = sum(!isnothing(option) for option in options)
  if n_given_options != 1
    throw(
      ArgumentError(
        "Only provide one of the optional parameters `min_counts`, " *
        "`min_genes`, `max_counts`, `max_genes` per call.",
      ),
    )
  end

  min_number = isnothing(min_genes) ? min_counts : min_genes
  max_number = isnothing(max_genes) ? max_counts : max_genes

  # Calculate the number of counts or genes per cell
  number_per_cell = if isnothing(min_genes) && isnothing(max_genes)
    vec(sum(data, dims=2))
  else
    vec(sum(data .> 0, dims=2))
  end

  # Apply filtering conditions
  if !isnothing(min_number)
    cell_subset = number_per_cell .>= min_number
  end
  if !isnothing(max_number)
    cell_subset = number_per_cell .<= max_number
  end

  # Output filtering information
  s = sum(.!cell_subset)
  if s > 0
    msg = "filtered out $s cells that have "
    if !isnothing(min_genes) || !isnothing(min_counts)
      msg *= "less than "
      msg *= min_counts === nothing ? "$min_genes genes expressed" : "$min_counts counts"
    end
    if !isnothing(max_genes) || !isnothing(max_counts)
      msg *= "more than "
      msg *= max_counts === nothing ? "$max_genes genes expressed" : "$max_counts counts"
    end
  else
    msg = "no cells passed the filtering threshold"
  end
  @info msg

  return cell_subset, number_per_cell
end

"""
    filter_cells!(
      data::Muon.AnnData;
      min_counts=nothing, min_genes=nothing,
      max_counts=nothing, max_genes=nothing
    ) -> Nothing

Filters cells in-place in a `Muon.AnnData` object.

# Arguments
- `data::Muon.AnnData`: The input single-cell data.
- `min_counts::Union{Int, Nothing}`: Minimum total counts per cell.
- `min_genes::Union{Int, Nothing}`: Minimum number of genes expressed per cell.
- `max_counts::Union{Int, Nothing}`: Maximum total counts per cell.
- `max_genes::Union{Int, Nothing}`: Maximum number of genes expressed per cell.

# Returns
- `Nothing`. The filtering is applied in-place.
"""
function filter_cells!(
  data::Muon.AnnData;
  min_counts::Union{Int, Nothing}=nothing,
  min_genes::Union{Int, Nothing}=nothing,
  max_counts::Union{Int, Nothing}=nothing,
  max_genes::Union{Int, Nothing}=nothing,
)::Nothing
  cell_subset, number = filter_cells(
    data;
    min_counts=min_counts,
    min_genes=min_genes,
    max_counts=max_counts,
    max_genes=max_genes,
  )
  # Update Muon.AnnData object
  if isnothing(min_genes) && isnothing(max_genes)
    data.obs[!, :n_counts] = number
  else
    data.obs[!, :n_genes] = number
  end

  # transform BitVector to something we can subset on
  if isa(cell_subset, BitVector)
    cell_subset = collect(cell_subset)
  end

  # filter adata in place
  data = subset_adata!(data, cell_subset, :cells)
  return nothing
end

#==============
# filter gene #
==============#
"""
    filter_genes(
      data::Muon.AnnData;
      min_counts=nothing, min_cells=nothing,
      max_counts=nothing, max_cells=nothing,
      copy=false
    ) -> Union{Tuple{BitVector, Vector}, Muon.AnnData}

Filters genes based on the given threshold criteria and returns either a subset mask and count vector or a new filtered `Muon.AnnData` object.

# Arguments
- `data::Muon.AnnData`: The input single-cell data.
- `min_counts::Union{Int, Nothing}`: Minimum total counts per gene.
- `min_cells::Union{Int, Nothing}`: Minimum number of cells expressing the gene.
- `max_counts::Union{Int, Nothing}`: Maximum total counts per gene.
- `max_cells::Union{Int, Nothing}`: Maximum number of cells expressing the gene.
- `copy::Bool`: If `true`, returns a filtered copy; otherwise, returns a mask and count vector.

# Returns
- If `copy == false`, returns a tuple `(gene_subset::BitVector, number_per_gene::Vector)`.
- If `copy == true`, returns a filtered `Muon.AnnData` object.
"""
function filter_genes(
  data::Muon.AnnData;
  min_counts::Union{Int, Nothing}=nothing,
  min_cells::Union{Int, Nothing}=nothing,
  max_counts::Union{Int, Nothing}=nothing,
  max_cells::Union{Int, Nothing}=nothing,
  copy::Bool=false,
)::Union{Tuple{BitVector, Vector}, Muon.AnnData}
  adata_copy = copy ? deepcopy(data) : data
  X = _choose_mtx_rep(adata_copy)

  gene_subset, number_per_gene = filter_genes(
    X;
    min_counts=min_counts,
    min_cells=min_cells,
    max_counts=max_counts,
    max_cells=max_cells,
  )
  if !copy
    return gene_subset, number_per_gene
  else
    if isnothing(min_cells) && isnothing(max_cells)
      adata_copy.var[!, :n_counts] = number_per_gene
    else
      adata_copy.var[!, :n_cells] = number_per_gene
    end
    # transform BitVector to something we can subset on
    if isa(gene_subset, BitVector)
      gene_subset = collect(gene_subset)
    end
    # filter adata in place 
    subset_adata!(adata_copy, gene_subset, :genes)
    return adata_copy
  end
end

"""
    filter_genes(
      data::AbstractMatrix;
      min_counts=nothing, min_cells=nothing,
      max_counts=nothing, max_cells=nothing
    ) -> Tuple{BitVector, Vector}

Filters genes from a count matrix based on the given threshold criteria.

# Arguments
- `data::AbstractMatrix`: Gene expression count matrix with genes as columns.
- `min_counts::Union{Int, Nothing}`: Minimum total counts per gene.
- `min_cells::Union{Int, Nothing}`: Minimum number of cells expressing the gene.
- `max_counts::Union{Int, Nothing}`: Maximum total counts per gene.
- `max_cells::Union{Int, Nothing}`: Maximum number of cells expressing the gene.

# Returns
- `gene_subset::BitVector`: A mask indicating genes that pass the filter.
- `number_per_gene::Vector`: A vector of counts or expressed cell numbers per gene.
"""
function filter_genes(
  data::AbstractMatrix;
  min_counts::Union{Int, Nothing}=nothing,
  min_cells::Union{Int, Nothing}=nothing,
  max_counts::Union{Int, Nothing}=nothing,
  max_cells::Union{Int, Nothing}=nothing,
)::Tuple{BitVector, Vector}

  # Check that only one filtering option is provided
  options = [min_cells, min_counts, max_cells, max_counts]
  n_given_options = sum(!isnothing(option) for option in options)
  if n_given_options != 1
    throw(
      ArgumentError(
        "Only provide one of the optional parameters `min_counts`, " *
        "`min_cells`, `max_counts`, `max_cells` per call.",
      ),
    )
  end

  min_number = isnothing(min_cells) ? min_counts : min_cells
  max_number = isnothing(max_cells) ? max_counts : max_cells

  # Calculate the number of counts or cells per gene
  number_per_gene = if isnothing(min_cells) && isnothing(max_cells)
    vec(sum(data, dims=1))
  else
    vec(sum(data .> 0, dims=1))
  end

  # Apply filtering conditions
  if !isnothing(min_number)
    gene_subset = number_per_gene .>= min_number
  end
  if !isnothing(max_number)
    gene_subset = number_per_gene .<= max_number
  end

  # Output filtering information
  s = sum(.!gene_subset)
  if s > 0
    msg = "Filtered out $s genes that have"
    if !isnothing(min_cells) || !isnothing(min_counts)
      msg *= " less than "
      msg *= min_counts === nothing ? "$min_cells cells expressed" : "$min_counts counts"
    end
    if !isnothing(max_cells) || !isnothing(max_counts)
      msg *= " more than "
      msg *= max_counts === nothing ? "$max_cells cells expressed" : "$max_counts counts"
    end
  else
    msg = "no genes passed the filtering threshold"
  end
  @info msg

  return gene_subset, number_per_gene
end

"""
    filter_genes!(
      data::Muon.AnnData;
      min_counts=nothing, min_cells=nothing,
      max_counts=nothing, max_cells=nothing
    ) -> Nothing

Filters genes in-place in a `Muon.AnnData` object.

# Arguments
- `data::Muon.AnnData`: The input single-cell data.
- `min_counts::Union{Int, Nothing}`: Minimum total counts per gene.
- `min_cells::Union{Int, Nothing}`: Minimum number of cells expressing the gene.
- `max_counts::Union{Int, Nothing}`: Maximum total counts per gene.
- `max_cells::Union{Int, Nothing}`: Maximum number of cells expressing the gene.

# Returns
- `Nothing`. The filtering is applied in-place.
"""
function filter_genes!(
  data::Muon.AnnData;
  min_counts::Union{Int, Nothing}=nothing,
  min_cells::Union{Int, Nothing}=nothing,
  max_counts::Union{Int, Nothing}=nothing,
  max_cells::Union{Int, Nothing}=nothing,
)::Nothing
  gene_subset, number = filter_genes(
    data;
    min_counts=min_counts,
    min_cells=min_cells,
    max_counts=max_counts,
    max_cells=max_cells,
  )
  # Update Muon.AnnData object
  if isnothing(min_cells) && isnothing(max_cells)
    data.var[!, :n_counts] = number
  else
    data.var[!, :n_cells] = number
  end

  # transform BitVector to something we can subset on
  if isa(gene_subset, BitVector)
    gene_subset = collect(gene_subset)
  end

  # filter adata in place 
  subset_adata!(data, gene_subset, :genes)

  return nothing
end
