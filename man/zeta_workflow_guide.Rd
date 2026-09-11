\name{zeta_workflow_guide}
\alias{zeta_workflow_guide}
\title{Guide 3: Input Checks and Multiple Candidates}
\description{Handle missingness explicitly, preserve identifiers, and separate
single-candidate results from a multiple-testing analysis.}
\details{
Default input checks stop on missing values, infinite values, negative biomarker
values, invalid treatment coding, or too few subjects. Do not turn factors into
numbers with \code{as.numeric(factor)}: that returns level codes. Recode treatments
explicitly and verify the order against the subject identifiers. The biomarker's
zero definition and any small-value recoding threshold must be prespecified.

With \code{na_action = "omit"}, complete rows are used and original row indices
are retained. The user must judge the missing-data assumptions; omission is not
an imputation or a correction for missingness. An all-zero biomarker or a
nonvariable outcome can make the test unevaluable. This is distinct from a
nonsignificant result.

For multiple candidates, prespecify the testing family and preserve every row,
including failed or unavailable tests. The example applies ordinary BH to the
evaluable global p values, retaining NA elsewhere. Its scope is the evaluable
family. Error-control assumptions still apply to these p values and their
dependence. This example does not implement the paper's extra external real-data
calibration and must not be labeled as reproducing its family q values.
}
\examples{
library(ZetaHTE)
set.seed(37)
n <- 160
a <- rep(0:1, n / 2)
x <- c(rep(0, 60), rexp(100))
y <- rnorm(n) + a * (0.2 + 2 * (x > 0))
candidate_2 <- c(rep(0, 50), rexp(110))
candidate_2[7] <- NA_real_

# Fail first, then make and document an explicit missing-data choice.
failure <- tryCatch(zeta_test(y, a, candidate_2), error = conditionMessage)
stopifnot(is.character(failure))
fit_2 <- suppressWarnings(zeta_test(y, a, candidate_2, B = 199, seed = 6,
                                    na_action = "omit"))
stopifnot(identical(fit_2$diagnostics$omitted_rows, 7L))
print(fit_2$diagnostics[c("n_input", "n_used", "omitted_rows")])

candidates <- list(candidate_1 = x, candidate_2 = candidate_2,
                   constant_candidate = rep(0, n))
results <- do.call(rbind, lapply(seq_along(candidates), function(i) {
  fit <- suppressWarnings(zeta_test(y, a, candidates[[i]], B = 199,
    seed = 100 + i, na_action = "omit"))
  data.frame(candidate = names(candidates)[i], n_used = fit$diagnostics$n_used,
    status = fit$status, reason = fit$reason, p_global = fit$p_global)
}))
ok <- is.finite(results$p_global)
results$p_BH <- NA_real_
results$p_BH[ok] <- p.adjust(results$p_global[ok], method = "BH")
results$family_size_evaluable <- sum(ok)
print(results)
stopifnot(nrow(results) == 3, is.na(results$p_BH[3]))

# This writes a synthetic summary only, then removes the temporary export.
destination <- tempfile(fileext = ".csv")
write.csv(results, destination, row.names = FALSE)
round_trip <- read.csv(destination)
stopifnot(identical(round_trip$candidate, results$candidate))
unlink(destination)
}
\section{Troubleshooting}{
Inspect \code{status}, \code{reason}, and \code{channels$reason} before interpreting
NA results. Increase B for finer Monte Carlo resolution; repeating plotting or
interpretation does not increase B. Keep seeds, RNG kind, package version, input
coding, omitted rows, and the candidate family with the analysis record. Use a
larger graphics device if a plot reports that the panel is too small.}
\seealso{\code{\link{zeta_test}}, \code{\link{zeta_reproduction}}}
