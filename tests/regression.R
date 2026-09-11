library(ZetaHTE)
n_checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  n_checks <<- n_checks + 1L
  cat("PASS:", label, "\n")
}
same <- function(x, y, tol = 1e-10) isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
quiet <- function(expr) suppressWarnings(expr)
fails <- function(expr, pattern) {
  msg <- tryCatch({force(expr); ""}, error = function(e) conditionMessage(e))
  nzchar(msg) && grepl(pattern, msg, fixed = TRUE)
}
internal <- function(name) getFromNamespace(name, "ZetaHTE")
set.seed(901)
x <- c(rep(0, 40), seq_len(80) / 80)
a <- rep(0:1, 60)
y <- rnorm(120) + a * x
fit <- quiet(zeta_test(y, a, x, B = 39, seed = 7, keep_bootstrap = TRUE))
check(identical(fit$settings$package_version, as.character(packageVersion("ZetaHTE"))), "recorded version matches installed package")
check(fit$diagnostics$outcome_offset == 0, "ordinary outcome is not numerically recentered")

tiny <- quiet(zeta_test(y, a, x * 1e-20, B = 39, seed = 7))
itiny <- quiet(zeta_interpret(tiny, cutoff = 5e-21, allow_nonsignificant = TRUE))
iorig <- quiet(zeta_interpret(fit, cutoff = 0.5, allow_nonsignificant = TRUE))
check(same(itiny$effects$estimate, iorig$effects$estimate), "tiny positive cutoff preserves descriptive groups and effects")
check(quiet(zeta_interpret(fit, alpha = 1e-20, allow_nonsignificant = TRUE))$alpha == 1e-20,
      "strictly positive alpha below machine epsilon accepted")
check(fails(zeta_interpret(fit, alpha = 1), "alpha must"), "alpha one rejected")
for (nm in c("y", "a", "x")) {
  d <- list(y = y, a = a, x = x, B = 39, seed = 7, na_action = "omit")
  other <- if (nm == "y") "x" else "y"
  d[[other]][1] <- NA_real_
  d[[nm]][1] <- Inf
  check(fails(do.call(zeta_test, d), "Infinite"), paste("infinite", nm, "not hidden by omission"))
}
offset <- quiet(zeta_test(y + 1e15, a, x, B = 39, seed = 7, keep_bootstrap = TRUE))
centered <- quiet(zeta_test((y + 1e15) - mean(y + 1e15), a, x, B = 39,
                           seed = 7, keep_bootstrap = TRUE))
check(offset$status == "ok", "large common offset does not imply no residual variation")
check(offset$diagnostics$outcome_offset != 0 && identical(offset$bootstrap, centered$bootstrap),
      "extreme-offset stabilization matches explicit centering of represented inputs")
for (off in c(0, 1e15)) {
  perfect <- quiet(zeta_test(off + a, a, x, B = 39, seed = 7))
  check(perfect$status == "not_evaluable" && is.na(perfect$p_global),
        paste("exact common-effect fit remains unavailable at offset", off))
}

old_kind <- RNGkind()
for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG", "Wichmann-Hill", "Super-Duper",
               "Knuth-TAOCP-2002")) {
  RNGkind(kind, "Inversion", "Rejection")
  set.seed(101)
  saved <- .Random.seed
  invisible(quiet(zeta_test(y, a, x, B = 19, seed = 42)))
  check(identical(saved, .Random.seed), paste("explicit seed restores", kind))
  f1 <- quiet(zeta_test(y, a, x, B = 19, seed = 42, keep_bootstrap = TRUE))
  set.seed(42)
  f2 <- quiet(zeta_test(y, a, x, B = 19, keep_bootstrap = TRUE))
  check(identical(f1$bootstrap, f2$bootstrap), paste("local and external seed agree for", kind))
}
RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
set.seed(29)
invisible(rnorm(1))
expected <- rnorm(1)
set.seed(29)
invisible(rnorm(1))
saved <- .Random.seed
check(fails(quiet(zeta_test(y, a, x, B = 19, seed = 1)), "hidden state"), "Box-Muller local seed is explicitly rejected")
check(identical(saved, .Random.seed) && identical(rnorm(1), expected), "rejection preserves seed and cached normal")
set.seed(29)
invisible(rnorm(1))
invisible(quiet(zeta_test(y, a, x, B = 19, seed = NULL)))
check(identical(rnorm(1), expected), "external-seed route preserves cached normal")
do.call(RNGkind, as.list(old_kind))

for (nz in c(0, 3, 4, 20)) {
  for (np in c(0, 7, 8, 9, 10, 14, 16, 30)) {
    if (nz + np < 4L) next
    xx <- c(rep(0, nz), seq_len(np))
    n <- length(xx)
    aa <- rep(0:1, length.out = n)
    dd <- internal("prepare_input")(seq_len(n), aa, xx, "continuous", 0, "fail")$data
    bb <- internal("make_basis")(dd)
    expected_jump <- nz >= 4 && np >= 8
    expected_linear <- np >= 8
    step_counts <- vapply(c(0.35, 0.5, 0.65), function(cut) {
      high <- dd$u[xx > 0] > cut
      sum(high) >= 5 && sum(!high) >= 5
    }, logical(1))
    expected <- c(if (expected_jump) "jump", if (expected_linear) "tail_linear",
                  if (any(step_counts)) paste0("step", which(step_counts)))
    check(identical(colnames(bb$basis), if (length(expected)) expected else NULL),
          paste("basis count boundaries: zero", nz, "positive", np))
  }
}

d <- internal("prepare_input")(y, a, x, "continuous", 0, "fail")$data
engine <- internal("prepare_engine")(d)
bases <- internal("make_basis")(d)$basis
x0 <- cbind(1, a, bases)
for (nm in names(engine$directions)) {
  cc <- engine$directions[[nm]]
  expected <- residuals(lm(I(a * bases[, nm]) ~ a + bases))
  check(same(cc[, 1], expected) && max(abs(crossprod(x0, cc))) < 1e-10,
        paste("residual direction and orthogonality:", nm))
}
set.seed(222)
signs <- sample(c(-1, 1), length(y), replace = TRUE)
signed <- engine$r * signs
rb <- qr.resid(engine$qr, signed)
pseudo_y <- (y - engine$r) + signed
check(same(rb, lm.fit(x0, pseudo_y)$residuals), "signed-residual refit equals pseudo-outcome null refit")
cc <- engine$directions[[1L]][, 1]
check(same(sum(signed * cc), sum(rb * cc)), "refit preserves score numerator against residualized direction")
check(abs(var(signed * cc) - var(rb * cc)) > 1e-8, "refit can change studentizing variance; not silently skipped")
set.seed(333)
order <- sample(seq_along(y))
de <- internal("prepare_engine")(d[order, ])
reb <- qr.resid(de$qr, de$r * signs[order])
check(same(internal("channel_scores")(rb, engine$directions),
           internal("channel_scores")(reb, de$directions)),
      "row reordering preserves scores when bootstrap signs follow the subjects")
for (n in c(60L, 120L, 480L)) {
  xx <- c(rep(0, n / 2), seq_len(n / 2))
  aa <- rep(0:1, n / 2)
  dd <- internal("prepare_input")(seq_len(n), aa, xx, "continuous", 0, "fail")$data
  bb <- internal("make_basis")(dd)$basis
  exact <- 1 + aa + as.numeric(bb %*% seq_len(ncol(bb)))
  check(quiet(zeta_test(exact, aa, xx, B = 19, seed = 1))$status == "not_evaluable",
        paste("exact null-basis fit is unavailable at n", n))
}
check(quiet(zeta_test(y, a, rep(1, length(x)), B = 19, seed = 1))$status == "not_evaluable",
      "all-positive constant biomarker is unavailable")
check(quiet(zeta_test(y, a, a, B = 19, seed = 1))$status == "not_evaluable",
      "biomarker perfectly confounded with treatment has no estimable direction")

# Direct enumeration checks inclusive tails, self-reference, and final minimum.
qb <- matrix(c(0, 1, 1, 3, 6, 0, 0, 2, 2, 5, 1, 1, 1, 4, 7), ncol = 3)
qo <- setNames(c(2, 3, 4), c("jump", "tail_linear", "tail_step_scan"))
wm <- (qo - colMeans(qb)) / apply(qb, 2, sd)
wb <- sweep(sweep(qb, 2, colMeans(qb), "-"), 2, apply(qb, 2, sd), "/")
obs_t <- c(max(wm), sum(pmax(wm, 0)) / sqrt(3))
boot_t <- cbind(apply(wb, 1, max), rowSums(pmax(wb, 0)) / sqrt(3))
obs_p <- vapply(1:2, function(j) (1 + sum(boot_t[, j] >= obs_t[j])) / 6, numeric(1))
self_p <- sapply(1:2, function(j) vapply(1:5, function(b) {
  (1 + sum(boot_t[, j] >= boot_t[b, j])) / 6
}, numeric(1)))
expected_p <- (1 + sum(apply(self_p, 1, min) <= min(obs_p))) / 6
combined <- internal("combine_channels")(qo, qb)
check(same(combined$combinations$p_value, obs_p), "combination proportions match direct enumeration")
check(same(combined$p_global, expected_p), "final minimum calibration matches direct enumeration")
all_ties <- internal("combine_channels")(qo * 0, qb * 0)
check(all_ties$p_global == 1 && length(all_ties$scale_fallback) == 3L,
      "all-tied bootstrap reference returns one with explicit scale fallback")

for (b in c(2L, 3L, 19L, 199L)) {
  f <- quiet(zeta_test(y, a, x, B = b, seed = 13))
  pp <- c(f$p_global, f$channels$p_value, f$combinations$p_value)
  check(all(pp > 0 & pp <= 1) && max(abs(pp * (b + 1) - round(pp * (b + 1)))) < 1e-10,
        paste("p bounds and resolution at B", b))
}

for (j in seq_len(20L)) {
  set.seed(6100 + j)
  n <- sample(c(60, 120, 240), 1)
  xx <- c(rep(0, floor(n * 0.4)), rexp(n - floor(n * 0.4)))
  aa <- sample(rep(0:1, length.out = n))
  yy <- rnorm(n) + aa * (0.4 + log1p(xx))
  f <- quiet(zeta_test(yy, aa, xx, B = 19, seed = j, keep_bootstrap = TRUE))
  rev <- quiet(zeta_test(yy, 1 - aa, xx, B = 19, seed = j, keep_bootstrap = TRUE))
  transformed <- quiet(zeta_test(-2 * yy + 3, aa, log1p(xx), B = 19, seed = j, keep_bootstrap = TRUE))
  check(same(f$bootstrap, rev$bootstrap) && same(f$channels$statistic, rev$channels$statistic) &&
          identical(f$p_global, rev$p_global), paste("treatment-label exchange property", j))
  check(same(f$bootstrap, transformed$bootstrap) && identical(f$p_global, transformed$p_global),
        paste("outcome/rank transformation property away from variance floor", j))
  ef <- quiet(zeta_interpret(f, allow_nonsignificant = TRUE))
  er <- quiet(zeta_interpret(rev, allow_nonsignificant = TRUE))
  check(same(ef$effects$estimate, -er$effects$estimate), paste("descriptive contrast signs under treatment swap", j))
}

tie_fit <- fit
tie_fit$scan$statistic[] <- 1
tie_summary <- quiet(zeta_interpret(tie_fit, allow_nonsignificant = TRUE))
check(tie_summary$selected_scan$rank_cutoff == 0.35, "exact scan ties select first ordered location")
wd <- getwd()
opts <- options()
temp <- tempfile(fileext = ".rds")
saveRDS(fit, temp)
restored <- readRDS(temp)
unlink(temp)
check(identical(fit, restored), "result survives RDS round trip")
check(identical(wd, getwd()) && identical(opts, options()), "calls preserve working directory and R options")
check(is.null(summary(restored)$data) && is.null(summary(restored)$call), "summary does not retain subject data or input call")
cat("\nTOTAL STEP 2 CHECKS:", n_checks, "\n")
