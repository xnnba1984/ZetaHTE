\name{zeta_plots}
\alias{zeta_plots}
\alias{plot.zeta_test}
\alias{plot.zeta_interpretation}
\title{Plot ZETA Test Results and Exploratory Summaries}
\description{Draw plots from an existing result without fitting or resampling.}
\usage{
\method{plot}{zeta_test}(x, type = c("channels", "scan"), alpha = 0.05,
  main = NULL, cex = 1, colors = "#5477C4", ...)
\method{plot}{zeta_interpretation}(x, type = c("groups", "contrasts"),
  main = NULL, cex = 1, colors = c("#2A9D8F", "#D5A445"), ...)
}
\arguments{
\item{x}{An object returned by \code{zeta_test} or \code{zeta_interpret}.}
\item{type}{Channel p values, observed scan statistics, treatment/control group
means, or contrasts between subgroup treatment effects.}
\item{alpha}{Exploratory reference line in the channel plot. This is not a
multiple-testing adjustment or a rule for selecting a unique channel.}
\item{main}{Optional plot title.}
\item{cex}{Positive text-size multiplier. Use a larger device if labels do not fit.}
\item{colors}{One channel/scan color, or two group colors in treatment/control
order. Group plots also distinguish arms by circle/square. Contrast points use
olive to distinguish descriptive differences from test evidence.}
\item{...}{Reserved. Extra arguments produce an error rather than being ignored.}
}
\details{
Plots draw on the current graphics device. A device of 7 by 5 inches is suitable
for most plots; use 7 by 6 inches for six group rows. All three channel rows and
three scan locations remain visible. Unavailable results are explicitly labeled.
Modified graphics parameters are restored even when drawing fails.

Group plots show means for continuous outcomes and response rates in percent
for binary outcomes. Treatment is coded 1 and control 0. Groups overlap: the
all-positive group contains both sides of each positive-tail split. The displayed
values are descriptive means, not adjusted predictions.

Contrasts compare treatment-minus-control differences across biomarker groups.
Jump is zero minus positive; step and user-cutoff contrasts are upper minus
lower. Binary contrasts are shown as percentage-point differences, with the
axis labeled response rate difference (percent). A value of 10 represents a
0.1 difference on the original probability scale, not a relative increase of
10 percent. The Tail linear slope has different units and remains in
\code{x$effects}; it is not placed on the contrast axis.

Plots add no confidence intervals. A descriptive override after a nonsignificant
global test remains labeled as such. Scan Q values are not p values and the
largest observed scan statistic does not identify a validated clinical cutoff.
}
\value{Invisibly returns the plotted table with display values and annotations.
The original result is unchanged. No subject records, files, or new random draws
are created by these methods.}
\examples{
set.seed(18)
x <- rep(c(0, seq(0.1, 1, length.out = 9)), each = 12)
a <- rep(0:1, length(x) / 2)
y <- rnorm(length(x), sd = 0.3) + a * (0.2 + 5 * x)
fit <- zeta_test(y, a, x, B = 199, seed = 42)
plot(fit)
plot(fit, type = "scan")
if (fit$p_global <= 0.05) {
  explanation <- zeta_interpret(fit)
  plot(explanation, type = "groups")
  plot(explanation, type = "contrasts")
}
}
\seealso{\code{\link{zeta_test}}, \code{\link{zeta_interpret}},
\code{\link{zeta_binary_guide}}}
