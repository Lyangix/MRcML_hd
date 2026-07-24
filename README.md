# MRcML\_hd

**High-Dimensional Mendelian Randomization via Constrained Maximum Likelihood**

`MRcML_hd` implements the high-dimensional variant of the [MRcML](https://github.com/xue-hr/MRcML) method for Mendelian Randomization (MR). It estimates the causal effect of an exposure on an outcome using GWAS summary statistics, while robustly handling many weak invalid instrumental variables (IVs).

The key difference from the standard `MRcML` package is the **standard error estimator**: this package uses a sandwich variance formula (instead of the standard Fisher information), which provides better-calibrated inference when the number of genetic instruments is large.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("Lyangix/MRcML_hd")
```

---

## Background

Mendelian Randomization (MR) uses genetic variants (SNPs) as instrumental variables (IVs) to estimate causal effects from observational data. Traditional MR methods rely on a small set of strong IVs, which can make results sensitive to the choice of SNPs, statistically unstable, and underpowered — particularly in settings where GWAS yield few genome-wide significant hits.

Incorporating a large number of IVs can recover more genetic information and improve statistical power, but introduces two new challenges: **weak-IV bias** from the aggregate stochastic fluctuation of many individually weak SNPs, and a higher risk of **invalid IVs** through pleiotropic effects.

`MRcML\_hd` generalizes the robust MR method [MRcML](https://github.com/xue-hr/MRcML) to this high-dimensional regime. We show that invalid-IV selection remains consistent despite the weak-IV setting, and establish the consistency and asymptotic normality of the resulting estimator. By accounting for pleiotropic effects across a large SNP set, the method aggregates weak genetic signals to improve estimation stability and statistical power. 

---

## Functions

### `mr_cML_hd()`

| Estimator | Description |
|---|---|
| `cML-BIC` | Estimate from the single BIC-selected *K* |

```r
mr_cML_hd(
  b_exp,          # GWAS effect sizes for exposure
  b_out,          # GWAS effect sizes for outcome
  se_exp,         # Standard errors for exposure
  se_out,         # Standard errors for outcome
  n,              # Sample size
  K_vec = 0:(length(b_exp) - 2),  # Candidate K values
  random_start = 0,
  maxit = 100,
  random_seed = 0
)
```

**Returns** a list with `BIC_theta`, `BIC_se`, `BIC_p`, `BIC_invalid`, `BIC_vec`.

---

## Quick Example

```r
library(MRcMLhd)

set.seed(42)
p   <- 500      # number of SNPs
n   <- 50000   # sample size
m   <- 100      # number of invalid IVs

# Simulate summary statistics (toy example)
b_exp <- rnorm(p, 0, 0.0005)
b_out <- 0.3 * b_exp + runif(m, 0.2, 0.3)
se_exp <- rep(1/sqrt(n), p)
se_out <- rep(1/sqrt(n), p)

# cML-BIC (without data perturbation)
res <- mr_cML_hd(b_exp, b_out, se_exp, se_out, n = n, random_seed = 1)
cat("cML-BIC     theta:", res$BIC_theta,    " SE:", res$BIC_se,    " p:", res$BIC_p,    "\n")
cat("Invalid IVs selected by cML-BIC:", res$BIC_invalid, "\n")
```
