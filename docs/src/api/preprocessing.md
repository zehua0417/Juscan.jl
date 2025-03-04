# Preprocessing Modules

Preprocessing is a crucial step in single-cell and multi-omics data analysis. This module provides functions for quality control, filtering, and normalization to ensure that datasets are clean and ready for downstream analyses.

Quality control (QC) metrics help identify low-quality cells and genes, such as those with extremely high or low counts, high dropout rates, or disproportionate expression levels. By computing these metrics, users can apply appropriate thresholds to retain only reliable data.

Filtering enables users to remove unwanted cells or genes based on QC criteria, ensuring that low-quality features do not affect downstream analyses.

Normalization methods adjust for differences in sequencing depth and technical variation, allowing meaningful comparisons across cells and conditions.

These preprocessing functions ensure that data is well-structured and suitable for clustering, dimensionality reduction, and other analytical tasks.

```@index
Pages = ["preprocessing.md"]
```

## quality control metrics

```@docs
calculate_qc_metrics
calculate_qc_metrics!
```

### sub functions

```@docs
describe_obs
describe_obs!
describe_var
describe_var!
```
