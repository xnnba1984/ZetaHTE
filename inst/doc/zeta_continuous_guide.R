### Name: zeta_continuous_guide
### Title: Guide 1: Continuous Outcomes
### Aliases: zeta_continuous_guide

### ** Examples

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



