\name{zeta_continuous_guide}
\alias{zeta_continuous_guide}
\title{Guide 1: Continuous Outcomes}
\description{Fit one biomarker candidate, inspect evidence, and obtain exploratory
effect summaries. The synthetic example is an illustration, not a paper result.}
\details{
Supply one outcome, treatment indicator, and biomarker value per subject, in the
same row order. Code treatment as 1, control as 0, and use a nonnegative biomarker
with a point mass at zero and a positive tail. Higher outcomes are treated as
favorable when describing effect signs. ZETA tests whether treatment effects
vary with the biomarker; a prognostic association alone is not the target.

Start with \code{p_global}. Channel p values describe exploratory directions
and can be small in more than one channel. They do not classify subjects or
establish separate confirmatory discoveries. Sparse/Dense values in
\code{summary(fit)} are auxiliary combination diagnostics, not signal labels.
An unavailable channel is returned with NA and a reason, not p = 1.

The example deliberately uses a clear positive-tail signal to illustrate all
outputs. Its small p value is not an estimate of power. The internal bootstrap
count is 999; the smallest attainable p value is 1/(999 + 1). An explicit seed
preserves the caller's random state with the supported default RNG.

The interpretation function requires a significant global result by default.
Jump compares the treatment effect in zero versus positive subjects. The
Tail linear slope is the coefficient of treatment by positive-tail rank in a
linear regression among positive subjects. The scan chooses the largest
observed statistic from ranks 0.35, 0.5, and 0.65. These are coarse exploratory
locations. Treatment effects in each group are mean(treatment) minus mean(control).
}
\examples{
library(ZetaHTE)
set.seed(18)
x <- rep(c(0, seq(0.1, 1, length.out = 9)), each = 12)
a <- rep(0:1, length(x) / 2)
y <- rnorm(length(x), sd = 0.3) + a * (0.2 + 5 * x)
fit <- zeta_test(y, a, x, outcome = "continuous", B = 999, seed = 42)
print(fit)
stopifnot(fit$status == "ok", all(fit$channels$evaluable))
plot(fit, type = "channels")
plot(fit, type = "scan")

if (fit$p_global <= 0.05) {
  explanation <- zeta_interpret(fit)
  print(explanation)
  print(explanation$groups)
  plot(explanation, type = "groups")
  plot(explanation, type = "contrasts")
  print(explanation$effects[explanation$effects$summary == "Tail linear slope", ])
}

# A summary omits the subject-level data retained in fit.
public_summary <- summary(fit)
stopifnot(is.null(public_summary$data), is.null(public_summary$bootstrap))
}
\section{Saving results}{
To save numerical summaries, use \code{write.csv(fit$channels, "channels.csv",
row.names = FALSE)} in a chosen output directory. A full \code{saveRDS(fit, ...)}
also saves subject-level data; protect that file. Plot methods draw on the active
R device and can be exported through the usual R graphics tools.}
\seealso{\code{\link{zeta_test}}, \code{\link{zeta_interpret}},
\code{\link{zeta_plots}}, \code{\link{zeta_binary_guide}}}
