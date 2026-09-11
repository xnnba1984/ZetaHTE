\name{zeta_test}
\alias{zeta_test}
\title{Test Treatment Effect Heterogeneity with a Semicontinuous Biomarker}
\description{The global ZETA test with Jump, Tail linear,
and Tail step scan channels, jointly calibrated by a residual wild bootstrap.}
\usage{
zeta_test(y, a, x, outcome = c("continuous", "binary"), B = 999L,
  seed = NULL, zero_threshold = 0,
  na_action = c("fail", "omit"), keep_bootstrap = FALSE)
}
\arguments{
\item{y}{Numeric outcome vector. Larger values represent more favorable outcomes.
Binary outcomes must be coded 0/1.}
\item{a}{Numeric or logical treatment vector, coded 1 for treatment and 0 for
control. At least two complete subjects in each arm are required.}
\item{x}{Finite nonnegative numeric biomarker vector. Factors and matrices are
not accepted. Average ranks are used for tied positive values.}
\item{outcome}{Outcome type. Binary outcomes use the same linear working model,
with effects interpreted on the risk-difference scale.}
\item{B}{Number of internal bootstrap samples, an integer of at least two.
The default is 999; settings below 199 warn about coarse resolution.}
\item{seed}{Optional nonnegative integer random seed. An explicit seed restores
the previous random-number state after the call. With NULL, the usual R random
stream is advanced. Explicit seeds are rejected with the Box-Muller normal
generator or a user-supplied RNG because their hidden state cannot safely be
restored. In those cases, call set.seed() yourself and use seed = NULL.}
\item{zero_threshold}{Prespecified nonnegative threshold. Values at or below
this threshold are recoded to zero; larger values retain their original scale.}
\item{na_action}{Fail on missing values, or omit rows incomplete in any of the
three input vectors. Omitted row indices are recorded. Infinity is an error.}
\item{keep_bootstrap}{Retain the matrix of bootstrap channel statistics.}
}
\details{
The five basis directions are Jump, Tail linear, and three positive-tail steps
at the fixed ranks 0.35, 0.5, and 0.65.
Jump requires at least four zero and eight positive subjects; Tail
linear requires eight positive subjects. Each step requires at least five
positive subjects on each side. Constant residualized directions are omitted.
These are inherited computational availability rules, not sample-size guarantees.
If at least one channel is available, joint calibration uses the available
channels. All three channel rows are returned, with NA and a reason for missing
channels. An entirely unevaluable test returns NA, not a nonsignificant p value.

Scores use the sample variance of residual products and a variance floor of
1e-8. Bootstrap samples independently flip residual signs, regress the signed
residuals on the original null predictors, and recompute all available channels.
Sparse and Dense statistics are auxiliary combinations, not signal classifiers.
Their bootstrap upper-tail proportions and the final lower-tail minimum
calibration use inclusive comparisons and plus-one corrections. For each
bootstrap combination statistic, its reference includes itself, exactly as in
the research runners. If a channel bootstrap standard deviation is at most
1e-8, its scale is set to 1 and a diagnostic is emitted, preserving the research
implementation.

The residual-variation check uses the centered outcome scale. If the absolute
outcome mean exceeds 1e8 times the largest centered absolute outcome, this
common offset is removed before the null fit and recorded in diagnostics.
An intercept makes this a numerical stabilization, not a change in the model.
This cannot recover precision already lost in the input. Very small outcome
scales can activate the variance floor. Use numerically reasonable outcome units.

The target is one biomarker and two treatments. Causal interpretation depends
on the study design; this function does not adjust for observational confounding.
It does not perform the paper's additional real-data permutation/shuffling
calibration, candidate-family multiplicity adjustment, or patient bootstrap.
}
\value{An object of class \code{zeta_test} containing:
\item{p_global, status, reason}{Global p value and evaluation status.}
\item{channels}{Three rows with identifiers, labels, observed statistics,
standardized evidence, exploratory p values, availability, and reasons.}
\item{scan}{Three prespecified rank locations with observed statistics and the
largest lower-group and smallest upper-group biomarker values. These boundary
values describe an interval, not a uniquely estimated clinical cutoff.}
\item{combinations}{Sparse and Dense statistics and auxiliary p values.}
\item{settings, diagnostics}{Parameters, random-number kind, sample counts,
omissions, model rank, and numerical diagnostics.}
\item{data}{Complete analysis rows, including original and recoded biomarker
values, ranks, and original row indices. This object contains subject-level data;
do not publish it without permission. Nothing is transmitted or saved by the function.}
\item{bootstrap}{Bootstrap channel matrix when requested, otherwise NULL.}
}
\examples{
set.seed(12)
x <- c(rep(0, 40), runif(80))
a <- rep(0:1, 60)
y <- rnorm(120) + a * (0.3 + 3 * x)
fit <- zeta_test(y, a, x, B = 199, seed = 42)
fit
fit$channels
summary(fit)
}
\seealso{\code{\link{zeta_interpret}}}
