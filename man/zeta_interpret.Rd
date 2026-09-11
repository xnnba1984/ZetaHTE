\name{zeta_interpret}
\alias{zeta_interpret}
\title{Describe Heterogeneity Evidence after a ZETA Test}
\description{Report exploratory treatment effect summaries and the scan location.
This function performs no additional resampling and supplies no confirmatory intervals.}
\usage{
zeta_interpret(object, alpha = 0.05, cutoff = NULL,
  cutoff_inclusive = c("lower", "upper"),
  allow_nonsignificant = FALSE)
}
\arguments{
\item{object}{An evaluable result from \code{zeta_test}.}
\item{alpha}{Global significance threshold, strictly between zero and one.
A global p value no larger than alpha passes the gate.}
\item{cutoff}{Optional positive biomarker-scale cutoff used only for descriptive
comparisons within the positive tail. It does not change the global test or replace
the prespecified internal rank scan.}
\item{cutoff_inclusive}{Which group contains values exactly equal to a
user-specified cutoff. The default, lower, uses X <= cutoff and X > cutoff.
Upper uses X < cutoff and X >= cutoff. Zero values are excluded from both groups.}
\item{allow_nonsignificant}{Explicitly permit descriptive output after a
nonsignificant global result, with a warning. Default FALSE.}
}
\details{
All three summaries are returned when estimable. Channel p values guide
interpretation without classifying the number of signal sources. The Jump
contrast is the zero-group treatment effect minus the positive-group treatment
effect. The Tail linear slope is the coefficient of treatment times positive
rank in a positive-subject regression with intercept, treatment, rank, and their
interaction. It requires more than eight positive subjects and a full-rank model.
The Tail step contrast is upper-group minus lower-group treatment effect, with
both groups restricted to positive biomarker values. The selected rank maximizes
the observed step score. Exact ties select the first rank in the ordered grid.
The returned biomarker boundaries describe the split; no interpolation or precise
cutpoint estimation is performed. The same upper-minus-lower convention applies
to a user-specified cutoff. Binary effects are probabilities, not percentage points.
If a subgroup lacks a treatment arm, its effect and the associated contrast are
NA with a reason. Independent confirmation is needed for localized claims.
}
\value{An object of class \code{zeta_interpretation} with global p value,
significance status, all channel results, an \code{effects} table, subgroup
counts and means in \code{groups}, \code{selected_scan}, and the explicit user
cutoff rule. Summaries are descriptive even after a significant global test.}
\examples{
set.seed(18)
x <- rep(c(0, seq(0.1, 1, length.out = 9)), each = 12)
a <- rep(0:1, length(x) / 2)
y <- rnorm(length(x), sd = 0.3) + a * (0.2 + 5 * x)
fit <- zeta_test(y, a, x, B = 199, seed = 42)
if (fit$p_global <= 0.05) {
  zeta_interpret(fit, cutoff = 0.6, cutoff_inclusive = "upper")
}
}
\seealso{\code{\link{zeta_test}}}
