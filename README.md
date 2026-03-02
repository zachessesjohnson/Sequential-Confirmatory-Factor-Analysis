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

This repository provides the preprint paper and (<!-- TODO: add links once code/scripts are published -->) associated analysis code for the Sequential CFA project.

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
├── README.md                              # This file
├── LICENSE                                # Apache-2.0 license
├── CITATION.cff                           # Machine-readable citation metadata
└── Sequential_CFA_Paper___OSF_Preprint.pdf  # Preprint manuscript (PDF)
```

> **TODO:** Once analysis scripts and/or data are added to this repository, update this section to describe the folder structure (e.g., `R/`, `data/`, `output/`, `simulations/`).

---

## Installation / Requirements

> **TODO:** This section will be updated once analysis scripts are published. The following is a suggested setup based on the methods described in the paper.

The analysis is expected to require **R** (the standard environment for CFA via packages such as `lavaan`). Suggested requirements:

- **R** ≥ 4.0.0
- R packages:
  - [`lavaan`](https://lavaan.ugent.be/) – for CFA model fitting
  - [`blavaan`](https://ecmerkle.github.io/blavaan/) – for Bayesian CFA
  - Additional packages for data manipulation and visualization (e.g., `tidyverse`, `ggplot2`)

### Installing R packages

```r
install.packages(c("lavaan", "blavaan", "tidyverse", "ggplot2"))
```

> **TODO:** Confirm exact R version and package versions used in the analysis and add a `renv.lock` or `sessionInfo()` output for reproducibility.

---

## Quick Start

> **TODO:** Update with exact steps once scripts are available.

Suggested workflow based on the paper's methodology:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/zachessesjohnson/Sequential-Confirmatory-Factor-Analysis.git
   cd Sequential-Confirmatory-Factor-Analysis
   ```

2. **Install required R packages** (see [Installation / Requirements](#installation--requirements)).

3. **Prepare your data** (see [Data](#data)).

4. **Run the sequential CFA analysis:**
   ```r
   # TODO: Replace with actual script path once published
   source("R/run_sequential_cfa.R")
   ```

5. **Inspect outputs** in the `output/` directory (see [Results](#results)).

---

## Usage Examples

> **TODO:** Add concrete usage examples once analysis scripts are available.

The following is a conceptual example of the sequential CFA workflow using `lavaan` in R:

```r
library(lavaan)

# Stage 1: Lower-level CFA on observed indicators
model_stage1 <- '
  subfactor1 =~ item1 + item2 + item3
  subfactor2 =~ item4 + item5 + item6
'
fit_stage1 <- cfa(model_stage1, data = lower_level_data)

# Extract factor scores from Stage 1
scores_stage1 <- lavPredict(fit_stage1)

# Stage 2: Upper-level CFA using Stage 1 factor scores as inputs
upper_level_data <- as.data.frame(scores_stage1)
model_stage2 <- '
  higher_factor =~ subfactor1 + subfactor2
'
fit_stage2 <- cfa(model_stage2, data = upper_level_data)

# Extract final index scores
index_scores <- lavPredict(fit_stage2)
```

### Monte Carlo Simulations

The paper includes Monte Carlo simulations comparing Sequential CFA against traditional and Bayesian CFA under varying:
- Sample sizes
- Model complexity levels
- Data skewness conditions

> **TODO:** Link to simulation scripts once published.

### WJP Rule of Law Index Application

The paper demonstrates the method using the World Justice Project (WJP) Rule of Law Index dataset.

> **TODO:** Add instructions for obtaining the WJP data and running the empirical application scripts.

---

## Reproducibility

> **TODO:** Add specific seeds, R version, and package version information once scripts are published.

To support reproducibility:

- Monte Carlo simulations should be run with a fixed random seed (e.g., `set.seed(12345)` in R).
- A `renv.lock` file or `sessionInfo()` snapshot is recommended to capture the exact package environment.
- All analysis scripts should be run in the order documented in the Quick Start section.

---

## Results

The main findings from the paper are:

- **Sequential CFA significantly outperforms traditional CFA** for indices with simple or moderate hierarchical complexity.
- **Traditional CFA performs better** when the underlying data are substantially skewed.
- **For highly complex hierarchical models**, the two approaches perform similarly.
- **Sequential CFA provides valid estimates** in cases where traditional or Bayesian CFA fail to converge entirely.

> **TODO:** Once output files are added to the repository, describe the output directory structure and file names produced by the scripts (e.g., factor score tables, simulation result CSVs, figures).

---

## Data

The empirical application in the paper uses the **WJP Rule of Law Index** dataset, produced by the [World Justice Project](https://worldjusticeproject.org/).

- **Availability:** The WJP dataset is publicly available at <https://worldjusticeproject.org/rule-of-law-index/>.
- **License:** Please refer to the WJP website for their data use terms.
- **Privacy:** The WJP data used in this analysis is aggregated at the country level and does not contain personally identifiable information.

Monte Carlo simulation data are generated synthetically within the analysis scripts and do not require external data sources.

> **TODO:** Confirm exact dataset version/year used and add download/access instructions once scripts are published.

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
