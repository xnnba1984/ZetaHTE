### Name: zeta_binary_guide
### Title: Guide 2: Binary Outcomes and a User Cutoff
### Aliases: zeta_binary_guide

### ** Examples

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



