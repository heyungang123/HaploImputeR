# HaploImputeR 0.1.0

* Initial CRAN release

## New Features

* `twinsGenerator()` - Generate synthetic haplotypes with automatic seed generation
* `twinsGeneratorM()` - Generate synthetic haplotypes with external seeds
* `twinsGeneratorRD()` - Generate synthetic haplotypes with dynamic window size
* `haploSeeds()` - Generate initial haplotype seeds
* `probComputationM()` - Compute probability matrix via Lasso regression
* `imputedHaplo2ExactC()` - Impute haplotypes with exact allele counts
* `winSizeR2()` - Compute dynamic window size based on R-squared
* `calculateLD()` - Calculate LD r-squared between two SNPs

## Improvements

* Support for parallel processing in `winSizeR2()`
* Reproducible results through random seed control
* Comprehensive input validation
* Memory-efficient matrix pre-allocation
* Full English documentation