### Name: zeta_workflow_guide
### Title: Guide 3: Input Checks and Multiple Candidates
### Aliases: zeta_workflow_guide

### ** Examples

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



