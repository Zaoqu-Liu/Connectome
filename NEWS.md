# Connectome 1.0.2

## New Features

* Added comprehensive vignettes:
  - Introduction to Connectome
  - Algorithm Principles and Mathematical Framework
  - Visualization Gallery
  - Differential Connectivity Analysis
  - Best Practices and Troubleshooting

* Built pkgdown documentation site

## Improvements

* Improved cross-platform compatibility using `Seurat::GetAssayData()` instead of direct slot access
* Replaced deprecated `Seurat::SubsetData()` with modern `subset()` approach
* Replaced `library()` calls with `requireNamespace()` checks for better CRAN compliance
* Fixed `Inf` value handling in differential analysis plots using `is.infinite()`
* Updated `ggplot2::element_rect()` from deprecated `size` to `linewidth`
* Replaced `1:length()` patterns with safer `seq_along()` for boundary condition handling
* Improved performance by replacing `rbind()` in loops with pre-allocation

## Bug Fixes

* Fixed scientific accuracy issue in DOR calculation by applying Haldane-Anscombe correction (pseudocount = 0.5)
* Fixed `PercentExpression` function to use `Seurat::Assays()` instead of internal function

## Documentation

* Updated author information (Zaoqu Liu as maintainer)
* Added GitHub repository URLs
* Cleaned up function documentation

---

# Connectome 1.0.1

* Initial release based on original Connectome package
* Multi-species ligand-receptor database support (human, mouse, rat, pig)
* Network construction and visualization functions
* Differential connectivity analysis
