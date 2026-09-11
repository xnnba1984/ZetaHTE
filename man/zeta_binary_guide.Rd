\name{zeta_binary_guide}
\alias{zeta_binary_guide}
\title{Guide 2: Binary Outcomes and a User Cutoff}
\description{Analyze a binary response and describe a prespecified positive-tail
split. All data in this guide are synthetic.}
\details{
Use 1 for a favorable response and 0 otherwise. The test uses the same linear
working model as for continuous outcomes. Heterogeneity is on the risk-difference
scale, not the odds-ratio scale. Outcome coding and effect scale matter.

A user cutoff changes descriptive grouping only. It does not change the global
test, the three channel tests, or the scan grid. Zero is kept separate. With
\code{cutoff_inclusive = "upper"}, the positive groups are
0 < X < cutoff and X >= cutoff. The default \code{"lower"} instead uses
0 < X <= cutoff and X > cutoff. Choose the boundary convention before inspecting
results. The numerical value 10 here is an example, not a clinical recommendation.

In the result tables, binary means and differences remain on the probability
scale. Group plots multiply means by 100; contrast plots multiply differences
by 100. No additional resampling occurs. Groups are overlapping descriptions,
not independent analyses. A group lacking one treatment arm has no estimable
treatment difference; its reason is retained.
}
\examples{
library(ZetaHTE)
set.seed(202)
x <- rep(c(0, 2, 5, 8, 10, 12, 15, 20, 25), each = 40)
a <- rep(0:1, length(x) / 2)
probability <- 0.15 + a * (0.05 + 0.6 * (x >= 10))
y <- rbinom(length(x), size = 1, prob = probability)
fit <- zeta_test(y, a, x, outcome = "binary", B = 999, seed = 73)
print(fit)
plot(fit)

if (fit$p_global <= 0.05) {
  explanation <- zeta_interpret(fit, cutoff = 10, cutoff_inclusive = "upper")
  print(explanation$user_cutoff_rule)
  print(explanation$groups)
  print(explanation$effects)
  stopifnot(sum(explanation$groups$n_treatment[5:6]) == sum(a == 1 & x > 0))
  stopifnot(sum(explanation$groups$n_control[5:6]) == sum(a == 0 & x > 0))
  plot(explanation, type = "groups", cex = 0.85)
  plot(explanation, type = "contrasts")
}
}
\section{Study interpretation}{
The primary target is a two-arm randomized study. These functions do not add
confounder adjustment to an observational comparison. A small global p value
supports further investigation of treatment effect heterogeneity; the selected
split and effect summaries remain exploratory. Patient-level stability bootstrap,
external calibration, and the paper's family q values are separate analyses.}
\seealso{\code{\link{zeta_interpret}}, \code{\link{zeta_plots}},
\code{\link{zeta_reproduction}}}
