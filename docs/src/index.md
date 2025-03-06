Juscan.jl is a Julia implementation of Scanpy, tailored for single-cell data analysis. Currently in its development (dev) version, the library aims to deliver high-performance, flexible, and extensible tools for preprocessing, dimensionality reduction, clustering, and visualization of single-cell datasets.

## Features

- **Single-Cell Data Analysis**  
  Provides functionalities for data normalization, log transformation, dimensionality reduction (e.g., PCA), and graph-based clustering (e.g., Leiden algorithm).

- **High Performance Computing**  
  Leverages Julia's computational strengths to efficiently handle large-scale single-cell datasets.

- **Modular and Extensible**  
  Designed with a modular architecture, allowing users to easily customize and extend functionality to suit diverse analytical requirements.

## Development Status

Juscan.jl is under active development. Some features are still being refined and may change over time. We welcome feedback, bug reports, and contributions from the community.

## Installation

### Installing the Development Version via Julia Package Manager

You can install the development version using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/zehua0417/Juscan.jl")
```

### Cloning the Repository for Local Development

```bash
git clone https://github.com/zehua0417/Juscan.jl.git
cd Juscan.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start Example

Here’s a simple example demonstrating basic single-cell data preprocessing and analysis with Juscan.jl:

```julia
using Juscan

## Load single-cell data (assuming the data is in AnnData format)
adata = Juscan.read_h5ad("data/sample.h5ad")

## Data preprocessing: normalization and calculate quality control metrics
Juscan.normalize!(adata)
Juscan.calculate_qc_metrics!(adata)

## filter cells and genes
filter_cells(adata, min_counts = 100)
filter_genes(adata, min_cells = 10)

## calc highly variable genes
highly_variable_genes!(adata)

## TODO: Dimensionality reduction and clustering
Juscan.pca!(adata)
Juscan.neighbors!(adata)
Juscan.leiden!(adata)

## TODO: Visualize the results
Juscan.plot_umap(adata)
```

## Contributing

We welcome contributions of all kinds:

- **Bug Reports & Suggestions**: Please open an issue on GitHub.
- **Code or Documentation Updates**.
- **Discussions**: Share your feedback and experiences to help improve the project.

## License

Juscan.jl is licensed under the MIT License.

## Contact

If you have any questions or suggestions, please open an issue on GitHub or contact us at [my e-mail](mailto:zehuali0417@gmail.com).

## Acknowledgments

We would like to extend our sincere gratitude to the open source community for their continuous support and inspiration. In particular, we thank:

- [Scanpy](https://github.com/theislab/scanpy) for its groundbreaking approach to single-cell analysis,
- [AnnData](https://github.com/theislab/anndata) for providing a robust data structure for annotated data,
- [Muon.jl](https://github.com/scverse/Muon.jl.git) for its innovative multi-omic analysis tools, and
- [scVI.jl](https://github.com/maren-ha/scVI.jl.git) for offering valuable insights into variational inference methods.

Their contributions and ideas have been instrumental in shaping the development and direction of Juscan.jl.

---

Thank you for your interest in Juscan.jl. We hope it will empower your single-cell data analysis projects!
