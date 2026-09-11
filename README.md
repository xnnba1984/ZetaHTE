# ZetaHTE

**Zero and Tail Treatment Effect Heterogeneity Test**

ZetaHTE implements the ZETA test for treatment effect heterogeneity with a
semicontinuous biomarker: a point mass at zero and a continuous positive tail.
It combines Jump, Tail linear, and Tail step scan channels into a global test
using a residual wild bootstrap. After a significant global test, it provides
exploratory channel evidence and descriptive effect summaries.

The package supports continuous and binary outcomes with two treatment groups.
Binary outcomes use the risk-difference scale. The primary target is a two-arm
randomized trial; the package does not adjust for observational confounding.

## Installation

Install the versioned release from GitHub:

```r
install.packages("remotes") # Run once if needed.
remotes::install_github("xnnba1984/ZetaHTE", ref = "v0.1.0")
library(ZetaHTE)
```

Alternatively, download `ZetaHTE_0.1.0.tar.gz` from the
[release page](https://github.com/xnnba1984/ZetaHTE/releases/tag/v0.1.0) and run:

```r
install.packages("ZetaHTE_0.1.0.tar.gz", repos = NULL, type = "source")
library(ZetaHTE)
```

R 4.1 or later is required. ZetaHTE uses only base/recommended R packages and
has no compiled components. Development and local testing use R 4.6.1 on macOS.

## Quick start

```r
set.seed(18)
x <- rep(c(0, seq(0.1, 1, length.out = 9)), each = 12)
a <- rep(0:1, length(x) / 2)
y <- rnorm(length(x), sd = 0.3) + a * (0.2 + 5 * x)

fit <- zeta_test(y, a, x, B = 999, seed = 42)
fit
fit$channels
plot(fit)

if (fit$p_global <= 0.05) {
  explanation <- zeta_interpret(fit, cutoff = 0.6)
  explanation
  explanation$groups
  plot(explanation, type = "contrasts")
}
```

`y` is the outcome, `a` is treatment (1) or control (0), and `x` is a
nonnegative numeric biomarker. Larger outcomes are treated as more favorable.
For binary outcomes, use `outcome = "binary"` and code a favorable response as 1.
Missing inputs cause an error by default; `na_action = "omit"` explicitly
removes incomplete rows and records their indices.

## Reading the output

- `fit$p_global` is the global ZETA test p value.
- `fit$channels` reports exploratory evidence for Jump, Tail linear, and
  Tail step scan. These p values do not establish separate confirmatory claims.
- `fit$scan` describes the three fixed positive-tail rank locations, 0.35,
  0.5, and 0.65. It does not estimate a precise clinical cutoff.
- `fit$combinations` contains Sparse and Dense combination diagnostics.
  They do not classify the number of signal sources.
- `zeta_interpret(fit)` reports descriptive effects after a significant test.
  A user-specified cutoff changes only the positive-tail descriptive groups,
  not the test or its p values. Independent confirmation is needed.
- An unavailable result is `NA` with a reason, not evidence of no signal.

Fit objects retain subject-level input data for interpretation. Keep them
private when working with sensitive data. The functions do not transmit data.
External real-data calibration, multiple-candidate adjustment, and patient
bootstrap analyses are separate workflows, not automatic package steps.

## Guides and plotting

The guides use synthetic data and are available as native R help, offline HTML,
and runnable R scripts:

| Guide | R help | Runnable source |
| --- | --- | --- |
| Continuous outcomes | `?zeta_continuous_guide` | [R script](inst/doc/zeta_continuous_guide.R) |
| Binary outcomes and a user cutoff | `?zeta_binary_guide` | [R script](inst/doc/zeta_binary_guide.R) |
| Input checks and multiple candidates | `?zeta_workflow_guide` | [R script](inst/doc/zeta_workflow_guide.R) |

```r
?zeta_test
?zeta_interpret
?zeta_plots
?zeta_reproduction
browseURL(system.file("doc", "index.html", package = "ZetaHTE"))
source(system.file("doc", "zeta_binary_guide.R", package = "ZetaHTE"))
```

`plot(fit)` displays channel evidence, and `plot(fit, type = "scan")` displays
observed scan statistics. For an interpretation result, `type = "groups"`
shows group means or response rates with sample counts, while
`type = "contrasts"` shows differences between subgroup treatment effects.
The Tail linear slope is retained in the effects table because it has different
units. Plot methods do not perform additional resampling. A 7 by 5 inch device
is suitable for most plots; use 7 by 6 inches for six group rows.

## Reproducibility and development

Record the package and R versions, `RNGkind()`, seed, bootstrap count `B`, and
input preprocessing. An explicit seed preserves the caller's random-number
state with supported built-in generators. Box-Muller and user-supplied
generators require external `set.seed()` and `seed = NULL`.

From the repository root, rebuild the offline documentation with
`Rscript tools/build-docs.R .`. From its parent directory, run:

```sh
R CMD build ZetaHTE
R CMD check --no-manual ZetaHTE_0.1.0.tar.gz
```

The check runs portable base-R tests for numerical behavior, input validation,
random-number handling, interpretation, and plotting. Examples and tests are
synthetic. This repository distributes the software, not a complete manuscript
reproduction bundle or clinical subject records.

## Citation, license, and contact

Use `citation("ZetaHTE")` for the software citation. ZetaHTE is released under
the [MIT license](LICENSE.md).

Author and maintainer: Miles Xi, [nxi@ucla.edu](mailto:nxi@ucla.edu).
Report reproducible software problems through
[GitHub Issues](https://github.com/xnnba1984/ZetaHTE/issues).
Do not post confidential or identifiable data in an issue.
