# Tools Module

The `tools` module provides essential utility functions to support various analytical workflows in single-cell and multi-omics data processing. These functions serve as building blocks for higher-level analysis by offering efficient, reusable operations that streamline common computational tasks.

### Key Functionalities

#### 1. Highly Variable Genes (HVG)

The module includes methods for identifying highly variable genes, which are crucial for downstream analyses such as clustering and dimensionality reduction. These functions help select informative genes by computing variability metrics across cells.

#### 2. Principal Component Analysis (PCA)

PCA is a widely used dimensionality reduction technique that captures the most significant variations in the dataset. The `tools` module provides efficient implementations for computing PCA, enabling users to reduce data complexity while preserving important biological signals.

#### 3. Clustering

The module supports clustering techniques to group cells based on gene expression patterns. These methods help uncover underlying cellular heterogeneity and identify distinct cell populations, facilitating biological interpretation of single-cell data.

By providing these fundamental tools, the `tools` module enhances data processing workflows, making it easier to perform robust, reproducible analyses across different stages of research.

```@index
Pages = ["tools.md"]
```

## highly variable genes

```@docs
subset_to_hvg!
highly_variable_genes
highly_variable_genes!
```
