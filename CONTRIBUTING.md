# Contributing to HaploImputeR

Thank you for your interest in contributing to HaploImputeR! This document provides guidelines and instructions for contributing.

## Types of Contributions

### Bug Reports

If you find a bug, please open an issue with:
- A clear, descriptive title
- A minimal reproducible example
- Expected behavior
- Actual behavior
- R session information (`sessionInfo()`)

### Feature Requests

We welcome feature requests! Please open an issue with:
- A clear description of the feature
- Use case and motivation
- Possible implementation approach (if you have ideas)

### Pull Requests

We actively welcome pull requests!

## Development Setup

1. **Fork and Clone**
   ```bash
   git clone https://github.com/YOUR_USERNAME/HaploImputeR.git
   cd HaploImputeR
   ```

2. **Install Dependencies**
   ```r
   install.packages(c("glmnet", "testthat", "knitr", "rmarkdown"))
   ```

3. **Run Tests**
   ```r
   testthat::test_dir("tests/testthat")
   ```

## Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable names
- Add roxygen2 documentation for all functions
- Include examples in documentation
- Write tests for new functionality

## Pull Request Process

1. Create a feature branch
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes
   - Update documentation
   - Add/update tests
   - Update NEWS.md

3. Run checks
   ```r
   devtools::check()
   ```

4. Commit changes
   ```bash
   git commit -m "Description of changes"
   ```

5. Push to your fork
   ```bash
   git push origin feature/your-feature-name
   ```

6. Open a Pull Request

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## Questions?

Feel free to open an issue for any questions about contributing.