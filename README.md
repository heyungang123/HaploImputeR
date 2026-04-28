# HaploImputeR

<!-- badges: start -->
[![R-CMD-check](https://github.com/heyungang123/HaploImputeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/heyungang123/HaploImputeR/actions/workflows/R-CMD-check.yaml)
[![Test Coverage](https://github.com/heyungang123/HaploImputeR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/heyungang123/HaploImputeR/actions/workflows/test-coverage.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/HaploImputeR)](https://CRAN.R-project.org/package=HaploImputeR)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Haplotypic Imputation and Generation Tools for R**

## Overview

HaploImputeR is an R package designed for haplotype imputation and synthetic haplotype generation. It uses Lasso regression and probabilistic methods to generate synthetic haplotypes based on reference population data while preserving allele frequency constraints.

### Key Features

- 🧬 **Synthetic Haplotype Generation** - Generate haplotypes matching target population allele frequencies
- 📊 **Dynamic Window Sizes** - Adaptive window sizes based on local LD structure
- ⚡ **Parallel Processing** - Efficient processing for large datasets
- 🔄 **Reproducibility** - Random seed control for reproducible results
- ✅ **Input Validation** - Comprehensive error checking and data validation

## Installation

### From CRAN (Coming Soon)

```r
install.packages("HaploImputeR")
```

### From GitHub

```r
# install.packages("devtools")
devtools::install_github("heyungang123/HaploImputeR")
```

### Dependencies

The package requires:
- R (>= 3.5.0)
- glmnet (required)

## Quick Start

### Basic Usage

```r
library(HaploImputeR)

# Create reference population (100 sites x 1000 chromosomes)
set.seed(42)
hap_ref <- matrix(sample(0:1, 100000, replace = TRUE), nrow = 100, ncol = 1000)

# Define target population allele counts (100 sites, 100 chromosomes)
count_obj <- matrix(c(50, 50), nrow = 100, ncol = 2)

# Generate synthetic haplotypes
simu_hap <- twinsGenerator(
  haplotype_ref = hap_ref,
  count_alleles_obj = count_obj,
  seed_size = 10,
  batch_size = 2,
  prior_window_size = 50,
  seed = 42
)

# Verify results
dim(simu_hap)  # 100 x 100
```

## Main Functions

| Function | Description |
|----------|-------------|
| `twinsGenerator()` | Generate haplotypes with automatic seed generation |
| `twinsGeneratorM()` | Generate haplotypes with external seeds |
| `twinsGeneratorRD()` | Generate haplotypes with dynamic window size |
| `haploSeeds()` | Generate initial haplotype seeds |
| `probComputationM()` | Compute probability matrix via Lasso |
| `imputedHaplo2ExactC()` | Impute haplotypes with exact allele counts |
| `winSizeR2()` | Compute dynamic window size based on R² |
| `calculateLD()` | Calculate LD r² between two SNPs |

## Documentation

- [Tutorial Vignette](vignettes/HaploImputeR-tutorial.Rmd)
- [Function Reference](man/)
- [NEWS](NEWS.md)

## Algorithm

The haplotype generation process:

1. **Seed Generation**: Initial haplotype patterns based on reference population frequencies
2. **Probability Estimation**: Lasso regression predicts haplotype probabilities
3. **Exact Matching**: Alleles assigned to match target counts exactly
4. **Extension**: Process continues until all sites are filled

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use HaploImputeR in your research, please cite:

```bibtex
@Manual{HaploImputeR,
  title = {HaploImputeR: Haplotype Imputation and Generation Tools},
  author = {YUQI CHEN and YUNGANG HE},
  year = {2024},
  note = {R package version 0.1.0},
  url = {https://github.com/heyungang123/HaploImputeR},
}
```

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## Contact

- **Issues**: [GitHub Issues](https://github.com/heyungang123/HaploImputeR/issues)
- **Email**: 625108562@qq.com