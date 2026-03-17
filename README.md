# MRcMLhd

**High-Dimensional Mendelian Randomization via Constrained Maximum Likelihood**

`MRcMLhd` is an R package that implements the high-dimensional variant of the [MRcML](https://github.com/xue-hr/MRcML) method for Mendelian Randomization (MR). It estimates the causal effect of an exposure on an outcome using GWAS summary statistics, while robustly handling invalid instrumental variables (IVs).

The key difference from the standard `MRcML` package is the **standard error estimator**: this package uses a RAPS-style sandwich variance formula (instead of the standard Fisher information), which provides better-calibrated inference when the number of genetic instruments is large.

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

`MRcMLhd` generalizes the robust MR method [MRcML](https://github.com/xue-hr/MRcML) to this high-dimensional regime. We show that invalid-IV selection remains consistent despite the weak-IV setting, and establish the consistency and asymptotic normality of the resulting estimator. By accounting for pleiotropic effects across a large SNP set, the method aggregates weak genetic signals to improve estimation stability and statistical power. 

---

## Functions

### `mr_cML_hd()`

Point estimation without data perturbation. Returns two estimators that differ in how the number of invalid IVs *K* is selected:

| Estimator | Description |
|---|---|
| `cML-MA-BIC` | Model-averaged estimate using BIC weights (recommended) |
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

**Returns** a list with `MA_BIC_theta`, `MA_BIC_se`, `MA_BIC_p`, `BIC_theta`, `BIC_se`, `BIC_p`, `BIC_invalid`, `BIC_vec`.

---

### `mr_cML_DP_hd()`

Estimation with **data perturbation** for additional robustness. Adds DP variants of the BIC estimators and two goodness-of-fit (GOF) tests.

```r
mr_cML_DP_hd(
  b_exp,               # GWAS effect sizes for exposure
  b_out,               # GWAS effect sizes for outcome
  se_exp,              # Standard errors for exposure
  se_out,              # Standard errors for outcome
  n,                   # Sample size
  K_vec = 0:(length(b_exp) - 2),
  random_start = 0,
  random_start_pert = 0,
  maxit = 100,
  num_pert = 200,      # Number of perturbation resamples
  random_seed = 0
)
```

**Returns** a list with `MA_BIC_theta/se/p`, `BIC_theta/se/p`, `BIC_invalid`, `MA_BIC_DP_theta/se/p`, `BIC_DP_theta/se/p`, `GOF1_p`, `GOF2_p`.

The two GOF p-values test whether the selected *K* is correctly specified. A small p-value suggests potential model misspecification.

---

## Quick Example

```r
library(MRcMLhd)

set.seed(42)
p   <- 50      # number of SNPs
n   <- 50000   # sample size

# Simulate summary statistics (toy example)
b_exp <- rnorm(p, 0, 0.05)
b_out <- 0.3 * b_exp + rnorm(p, 0, 0.01)
se_exp <- rep(0.01, p)
se_out <- rep(0.01, p)

# cML-BIC (without data perturbation)
res <- mr_cML_hd(b_exp, b_out, se_exp, se_out, n = n, random_seed = 1)
cat("cML-MA-BIC  theta:", res$MA_BIC_theta, " SE:", res$MA_BIC_se, " p:", res$MA_BIC_p, "\n")
cat("cML-BIC     theta:", res$BIC_theta,    " SE:", res$BIC_se,    " p:", res$BIC_p,    "\n")
cat("Invalid IVs selected by cML-BIC:", res$BIC_invalid, "\n")

# cML-BIC-DP (with data perturbation, slower but more robust)
res_dp <- mr_cML_DP_hd(b_exp, b_out, se_exp, se_out, n = n,
                        num_pert = 200, random_seed = 1)
cat("cML-MA-BIC-DP  theta:", res_dp$MA_BIC_DP_theta, " SE:", res_dp$MA_BIC_DP_se, " p:", res_dp$MA_BIC_DP_p, "\n")
cat("cML-BIC-DP     theta:", res_dp$BIC_DP_theta,    " SE:", res_dp$BIC_DP_se,    " p:", res_dp$BIC_DP_p,    "\n")
cat("GOF test p-values: GOF1 =", res_dp$GOF1_p, "  GOF2 =", res_dp$GOF2_p, "\n")
```


---

## License

MIT
