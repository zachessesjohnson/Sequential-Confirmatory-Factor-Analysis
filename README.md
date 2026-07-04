# Sequential Confirmatory Factor Analysis

> Sequential CFA is a novel approach to hierarchical Confirmatory Factor Analysis that estimates each stage independently, allowing for estimation of small-N samples.

---

## Table of Contents

1. [Overview](#overview)
2. [Background & Motivation](#background--motivation)
3. [Method Summary](#method-summary)
4. [Repository Structure](#repository-structure)
5. [Installation / Requirements](#installation--requirements)
6. [Quick Start](#quick-start)
7. [Usage Examples](#usage-examples)
8. [Reproducibility](#reproducibility)
9. [Results](#results)
10. [Data](#data)
11. [Citation](#citation)
12. [License](#license)
13. [Contributing](#contributing)
14. [Contact](#contact)

---

## Overview

**Sequential Confirmatory Factor Analysis (Sequential CFA)** is a novel statistical method for constructing hierarchical factor indices—such as national or sub-national composite measures—when sample sizes at upper levels of the hierarchy are small. Traditional hierarchical CFA attempts to estimate all levels of a factor model simultaneously, which can produce parameter bias and convergence failures under small-N conditions. Sequential CFA resolves this by estimating each level of the hierarchy independently, from the lowest to the highest, preserving the full sample size at each stage and eliminating the need to estimate cross-level covariances in a single model pass.

This repository provides the preprint paper and associated R package (`scfa`) implementing the Sequential CFA workflow.

**Preprint:** <https://osf.io/preprints/osf/akxtv_v2>

---

## Background & Motivation

Factor score estimation in small sample sizes often encounters parameter bias and convergence failures when constructing hierarchical national/sub-national indices. Many real-world composite indices—such as rule-of-law or governance indices—are built from survey data aggregated at the country or region level. At these upper levels, the effective sample size (number of countries/regions) can be very small (e.g., 30–120 units), making simultaneous multi-level CFA estimation unreliable or impossible.

This paper proposes **sequential Confirmatory Factor Analysis** as a principled solution. Instead of estimating multiple levels of factors simultaneously, this approach calculates factor scores sequentially from the lowest to highest levels. This sequential estimation:

- Keeps the original sample size in each estimation step.
- Removes the need to estimate cross-level covariances within a single model.
- Provides valid estimates in many cases where traditional or Bayesian CFA fail to converge.

### Abstract

Factor score estimation in small sample sizes often encounters parameter bias and convergence failures when constructing hierarchical national/sub-national indices. This paper proposes a novel method for hierarchical factor analysis called "sequential Confirmatory Factor Analysis". Instead of estimating multiple levels of factors at the same time, this approach calculates factor scores sequentially from the lowest to highest levels. This sequential estimation keeps the original sample size in each step and also removes cross-level covariance estimation. Using a series of Monte Carlo simulations, we isolate the difference between sequential Confirmatory Factor Analysis and traditional Confirmatory Factor Analysis by comparing their resulting factor scores to the true latent variables under varying conditions. We also estimate the WJP Rule of Law Index using traditional Confirmatory Factor Analysis, Bayesian Confirmatory Factor Analysis, and sequential Confirmatory Factor Analysis to test performance. Our findings demonstrate that sequential Confirmatory Factor Analysis significantly outperforms the traditional model for indices with simple/moderate complexity. Traditional Confirmatory Factor Analysis performs better where the data are skewed. Where the hierarchical model becomes complex, the two methods perform similarly. Finally, sequential Confirmatory Factor Analysis can provide valid estimates where traditional or Bayesian Confirmatory Factor Analysis fail to converge.

---

## Method Summary

Sequential CFA is a multi-stage hierarchical factor analysis approach. The high-level workflow is:

1. **Stage 1 – Lower-level CFA:** Fit a standard CFA model to the observed indicators at the lowest level of the hierarchy (e.g., individual survey items within a sub-factor). Extract and save the resulting factor scores.
2. **Stage 2 – Upper-level CFA:** Use the saved factor scores from Stage 1 as observed inputs to a new CFA model at the next level up. Repeat for each additional level.
3. **Final Stage – Index Construction:** Aggregate the highest-level factor scores to produce the composite index.

Because each stage is estimated independently:
- The full observed-data sample size is retained at each stage.
- Cross-level covariance structures do not need to be parameterized simultaneously.
- Standard CFA software can be used at each stage without modification.

The approach is benchmarked against:
- Traditional (simultaneous) Confirmatory Factor Analysis.
- Bayesian Confirmatory Factor Analysis.

Performance is evaluated using Monte Carlo simulations and a real-world application to the **WJP Rule of Law Index**.

---

## Repository Structure

```
Sequential-Confirmatory-Factor-Analysis/
├── README.md                                # This file
├── DESCRIPTION                              # R package metadata
├── NAMESPACE                                # R package exports
├── LICENSE                                  # Apache-2.0 license
├── CITATION.cff                             # Machine-readable citation metadata
├── Sequential_CFA_Paper___OSF_Preprint.pdf  # Preprint manuscript (PDF)
├── R/
│   ├── scfa_propagation.R    # Error-propagation diagnostics & corrections
│   └── run_scfa.R            # High-level orchestration function
├── man/                      # Roxygen2-generated documentation
├── vignettes/
│   └── sequential-cfa-workflow.Rmd  # End-to-end tutorial
└── tests/
    └── testthat/             # testthat unit & integration tests
```

---

## Installation / Requirements

This repository is structured as a standard R package. Install it directly from
GitHub with:

```r
# install.packages("remotes")
remotes::install_github("zachessesjohnson/Sequential-Confirmatory-Factor-Analysis")
```

Alternatively, clone the repository and install locally:

```r
# From the repository root
remotes::install_local(".")
```

**Requirements:**

- **R** ≥ 4.0.0
- R packages (installed automatically when using `remotes::install_github()`):
  - [`lavaan`](https://lavaan.ugent.be/) ≥ 0.6 – CFA model fitting at each stage
  - [`ggplot2`](https://ggplot2.tidyverse.org/) *(optional)* – enhanced diagnostic plots

For Bayesian CFA comparisons described in the paper:

```r
install.packages("blavaan")
```

---

## Quick Start

1. **Install the package:**
   ```r
   remotes::install_github("zachessesjohnson/Sequential-Confirmatory-Factor-Analysis")
   ```

2. **Load the package and run the full workflow with one call:**
   ```r
   library(lavaan)
   library(scfa)

   result <- run_scfa(
     stage_models = list(
       "subfactor1 =~ x1 + x2 + x3
        subfactor2 =~ x4 + x5 + x6",
       "index =~ subfactor1 + subfactor2"
     ),
     data      = your_data,
     method    = "bartlett",
     threshold = 0.70
   )

   print(result)
   head(result$index_scores)
   ```

3. **Inspect propagation diagnostics:**
   ```r
   print(result$diagnostics$stage_1)
   plot(result$diagnostics$stage_1)
   ```

4. **Read the vignette** for a complete walk-through:
   ```r
   vignette("sequential-cfa-workflow", package = "scfa")
   ```

---

## Usage Examples

### One-call orchestration with `run_scfa()`

```r
library(lavaan)
library(scfa)

result <- run_scfa(
  stage_models = list(
    "subfactor1 =~ item1 + item2 + item3
     subfactor2 =~ item4 + item5 + item6",
    "index =~ subfactor1 + subfactor2"
  ),
  data      = lower_level_data,
  method    = "bartlett",
  threshold = 0.70
)

# Summary
print(result)

# Stage-1 propagation diagnostics
print(result$diagnostics$stage_1)
plot(result$diagnostics$stage_1)

# Final index scores
head(result$index_scores)

# Residual variance correction (Bartlett workflow)
result$correction
```

### Error-propagation diagnostics (recommended pre-Stage-2 check)

Run the diagnostics on your fitted Stage-1 model to assess how much estimation
error will propagate forward before you ever run Stage 2.

```r
library(lavaan)
library(scfa)

# --- Stage 1: fit lower-level CFA ---
model_stage1 <- '
  subfactor1 =~ item1 + item2 + item3
  subfactor2 =~ item4 + item5 + item6
'
fit_stage1 <- cfa(model_stage1, data = lower_level_data)

# Full diagnostic table (one row per factor)
diag <- scfa_propagation_diagnostics(fit_stage1)
print(diag)
#   factor n_indicators  I_k psi_nu phi_k rho_k  flag
# 1     f1            3 8.42   0.12   1.0  0.894    ok
# 2     f2            3 3.17   0.32   1.0  0.760    ok

# Visualise reliabilities
plot(diag)

# Individual quantities
scfa_factor_information(fit_stage1)   # I_k
scfa_propagation_variance(fit_stage1) # psi_nu_k = 1 / I_k
scfa_factor_reliability(fit_stage1)   # rho_k
```

A factor with `rho_k < 0.70` (the default `threshold` in
`scfa_propagation_diagnostics()`) is flagged as a weak link where propagation
is material.

### Correcting Stage-2 loadings for regression-score attenuation

If you used **regression** factor scores as Stage-2 inputs, the loadings are
attenuated by a factor of `rho_k`.  Recover unbiased estimates with:

```r
# --- Stage 2: upper-level CFA using Stage-1 scores as inputs ---
upper_level_data <- as.data.frame(lavPredict(fit_stage1, method = "regression"))
model_stage2 <- '
  higher_factor =~ subfactor1 + subfactor2
'
fit_stage2 <- cfa(model_stage2, data = upper_level_data)

# De-attenuate the Stage-2 loadings
corrected_lambda <- scfa_correct_loadings(fit_stage2, fit_stage1)
```

### Correcting Stage-2 residual variances for Bartlett-score inflation

Under **Bartlett** scores the loadings are asymptotically unbiased, but the
propagated error variance `psi_nu_k = 1 / I_k` still inflates Stage-2 residual
variances. Recover the adjusted residuals with:

```r
scores_bart <- as.data.frame(lavPredict(fit_stage1, method = "bartlett"))
fit_stage2  <- cfa("higher_factor =~ subfactor1 + subfactor2",
                    data = scores_bart)
adj_theta   <- scfa_correct_residuals(fit_stage2, fit_stage1)
```

### Multi-stage (> 2 levels) propagation chain

For three or more levels, `scfa_propagate_chain()` accumulates propagation
variances stage by stage:

```r
chain <- scfa_propagate_chain(list(fit_stage1, fit_stage2, fit_stage3))
# Cumulative propagation variance entering Stage 3
chain[["stage_3"]]
```

### Complete sequential CFA workflow

```r
library(lavaan)
library(scfa)

# Stage 1
model_stage1 <- '
  subfactor1 =~ item1 + item2 + item3
  subfactor2 =~ item4 + item5 + item6
'
fit_stage1 <- cfa(model_stage1, data = lower_level_data)

# Pre-Stage-2 propagation check
diag <- scfa_propagation_diagnostics(fit_stage1, threshold = 0.70)
print(diag)

# Stage 2
scores_stage1    <- as.data.frame(lavPredict(fit_stage1))
model_stage2     <- 'higher_factor =~ subfactor1 + subfactor2'
fit_stage2       <- cfa(model_stage2, data = scores_stage1)

# Final index scores
index_scores <- lavPredict(fit_stage2)
```

### Monte Carlo Simulations

The paper includes Monte Carlo simulations comparing Sequential CFA against traditional and Bayesian CFA under varying:
- Sample sizes
- Model complexity levels
- Data skewness conditions

Simulation scripts will be added to this repository once the paper is accepted for publication.

### WJP Rule of Law Index Application

The paper demonstrates the method using the World Justice Project (WJP) Rule of Law Index dataset.

- **Obtain the data:** Download the dataset from <https://worldjusticeproject.org/rule-of-law-index/>.
- Application scripts will be added to this repository once the paper is accepted for publication.

---

## Reproducibility

To support reproducibility:

- Monte Carlo simulations use a fixed random seed (`set.seed(12345)` in R).
- To capture the exact package environment used when analysing results, run
  `renv::snapshot()` after installing all dependencies; commit the resulting
  `renv.lock` file.
- All analysis scripts should be run in the order documented in the Quick Start section.

---

## Results

The main findings from the paper are:

- **Sequential CFA significantly outperforms traditional CFA** for indices with simple or moderate hierarchical complexity.
- **Traditional CFA performs better** when the underlying data are substantially skewed.
- **For highly complex hierarchical models**, the two approaches perform similarly.
- **Sequential CFA provides valid estimates** in cases where traditional or Bayesian CFA fail to converge entirely.

Output files (factor score tables, simulation result CSVs, figures) will be added to an `output/` directory once scripts are published.

---

## Data

The empirical application in the paper uses the **WJP Rule of Law Index** dataset, produced by the [World Justice Project](https://worldjusticeproject.org/).

- **Availability:** The WJP dataset is publicly available at <https://worldjusticeproject.org/rule-of-law-index/>.
- **License:** Please refer to the WJP website for their data use terms.
- **Privacy:** The WJP data used in this analysis is aggregated at the country level and does not contain personally identifiable information.

Monte Carlo simulation data are generated synthetically within the analysis scripts and do not require external data sources.

---

## Citation

If you use this method or code in your work, please cite the preprint:

> zachessesjohnson (2024). *Sequential Confirmatory Factor Analysis*. OSF Preprints. <https://osf.io/preprints/osf/akxtv_v2>

A machine-readable [`CITATION.cff`](./CITATION.cff) file is included in this repository.

---

## License

This project is licensed under the **Apache License 2.0**. See the [`LICENSE`](./LICENSE) file for details.

---

## Contributing

Contributions, bug reports, and suggestions are welcome! Please open an issue or submit a pull request on [GitHub](https://github.com/zachessesjohnson/Sequential-Confirmatory-Factor-Analysis).

---

## Contact

For questions about this project, please open an issue on the [GitHub repository](https://github.com/zachessesjohnson/Sequential-Confirmatory-Factor-Analysis/issues) or contact the author via their GitHub profile: [@zachessesjohnson](https://github.com/zachessesjohnson).
