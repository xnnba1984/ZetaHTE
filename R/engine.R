zeta_cuts <- function() c(0.35, 0.5, 0.65)
zeta_channels <- function() c("jump", "tail_linear", "tail_step_scan")
zeta_labels <- function() c("Jump", "Tail linear", "Tail step scan")

drop_constant <- function(x) {
  if (!ncol(x)) return(x)
  keep <- apply(x, 2L, stats::sd) > 1e-10
  x[, keep, drop = FALSE]
}

make_basis <- function(d) {
  n <- nrow(d)
  z <- d$x == 0
  pos <- !z
  pi_zero <- mean(z)
  b <- matrix(numeric(0), nrow = n, ncol = 0L)
  add <- function(v, nm) {
    out <- cbind(b, v)
    colnames(out)[ncol(out)] <- nm
    out
  }
  if (sum(z) >= 4L && sum(pos) >= 8L) {
    b <- add((as.numeric(z) - pi_zero) / sqrt(pi_zero * (1 - pi_zero)), "jump")
  }
  if (sum(pos) >= 8L) {
    v <- numeric(n)
    v[pos] <- sqrt(3) * (2 * d$u[pos] - 1) / sqrt(1 - pi_zero)
    b <- add(v, "tail_linear")
  }
  scan <- data.frame(rank_cutoff = zeta_cuts(), statistic = NA_real_,
                     evaluable = FALSE, lower_max = NA_real_, upper_min = NA_real_,
                     reason = "Fewer than five positive subjects on one side.")
  for (j in seq_len(nrow(scan))) {
    high <- d$u[pos] > scan$rank_cutoff[j]
    if (any(!high)) scan$lower_max[j] <- max(d$x[pos][!high])
    if (any(high)) scan$upper_min[j] <- min(d$x[pos][high])
    if (sum(high) < 5L || sum(!high) < 5L) next
    ph <- mean(high)
    v <- numeric(n)
    v[pos] <- (as.numeric(high) - ph) / sqrt((1 - pi_zero) * ph * (1 - ph))
    b <- add(v, paste0("step", j))
    scan$reason[j] <- ""
  }
  list(basis = drop_constant(b), scan = scan)
}

prepare_engine <- function(d) {
  bases <- make_basis(d)
  x0 <- cbind(intercept = 1, treatment = d$a, bases$basis)
  if (any(!is.finite(x0))) stop("Nonfinite values in the null design.", call. = FALSE)
  centered_y <- d$y - mean(d$y)
  offset <- if (abs(mean(d$y)) > 1e8 * max(abs(centered_y))) mean(d$y) else 0
  # Remove an extreme common offset before fitting to avoid cancellation.
  working_y <- if (offset != 0) centered_y else d$y
  fit <- stats::lm.fit(x0, working_y)
  qrx <- qr(x0)
  r <- working_y - as.numeric(fit$fitted.values)
  if (any(!is.finite(r))) stop("Null model produced nonfinite residuals.", call. = FALSE)
  directions <- list()
  for (nm in colnames(bases$basis)) {
    cc <- qr.resid(qrx, bases$basis[, nm, drop = FALSE] * d$a)
    if (any(!is.finite(cc))) stop("Nonfinite interaction direction.", call. = FALSE)
    cc <- drop_constant(cc)
    if (ncol(cc)) directions[[nm]] <- cc
  }
  for (j in seq_len(nrow(bases$scan))) {
    available <- paste0("step", j) %in% names(directions)
    bases$scan$evaluable[j] <- available
    if (!available && bases$scan$reason[j] == "") {
      bases$scan$reason[j] <- "Residualized step direction is constant."
    }
  }
  list(r = r, qr = qrx, directions = directions, scan = bases$scan,
       rank = qrx$rank, df = nrow(d) - qrx$rank,
       basis_names = colnames(bases$basis), outcome_offset = offset)
}

direction_score <- function(r, cc) {
  products <- cc * as.numeric(r)
  score <- colSums(products) / sqrt(length(r))
  variance <- as.numeric(stats::cov(products))
  value <- as.numeric((score / sqrt(max(variance, 1e-8)))^2)
  if (!is.finite(value)) stop("Nonfinite score statistic; inspect the outcome scale.", call. = FALSE)
  value
}

channel_scores <- function(r, directions) {
  q <- vapply(directions, function(cc) direction_score(r, cc), numeric(1))
  out <- q[intersect(c("jump", "tail_linear"), names(q))]
  steps <- q[grepl("^step", names(q))]
  if (length(steps)) out["tail_step_scan"] <- max(steps)
  list(channels = out, steps = steps)
}

upper_p <- function(ref, observed) (1 + sum(ref >= observed)) / (length(ref) + 1)
lower_p <- function(ref, observed) (1 + sum(ref <= observed)) / (length(ref) + 1)

self_upper_p <- function(ref) {
  # Include the reference observation itself, matching the manuscript runners.
  n <- length(ref)
  (n - rank(ref, ties.method = "min") + 2) / (n + 1)
}

combine_channels <- function(observed, bootstrap) {
  center <- colMeans(bootstrap)
  scale <- apply(bootstrap, 2L, stats::sd)
  fallback <- !is.finite(scale) | scale <= 1e-8
  scale[fallback] <- 1
  w <- (as.numeric(observed) - center) / scale
  names(w) <- names(observed)
  wb <- sweep(sweep(bootstrap, 2L, center, "-"), 2L, scale, "/")
  ts <- max(w)
  td <- sum(pmax(w, 0)) / sqrt(length(w))
  tsb <- apply(wb, 1L, max)
  tdb <- rowSums(pmax(wb, 0)) / sqrt(ncol(wb))
  ps <- upper_p(tsb, ts)
  pd <- upper_p(tdb, td)
  pmin_boot <- pmin(self_upper_p(tsb), self_upper_p(tdb))
  list(p_global = lower_p(pmin_boot, min(ps, pd)), w = w,
       channel_p = vapply(seq_along(observed), function(j) {
         upper_p(bootstrap[, j], observed[j])
       }, numeric(1)),
       combinations = data.frame(combination = c("Sparse", "Dense"),
                                 statistic = c(ts, td), p_value = c(ps, pd)),
       center = center, scale = scale, scale_fallback = names(observed)[fallback],
       p_min = min(ps, pd))
}

zeta_test <- function(y, a, x, outcome = c("continuous", "binary"), B = 999L,
                      seed = NULL, zero_threshold = 0,
                      na_action = c("fail", "omit"), keep_bootstrap = FALSE) {
  call <- match.call()
  outcome <- match.arg(outcome)
  na_action <- match.arg(na_action)
  check_scalar(B, "B", 2, .Machine$integer.max, integer = TRUE)
  check_scalar(zero_threshold, "zero_threshold", 0)
  if (!is.null(seed)) check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  check_flag(keep_bootstrap, "keep_bootstrap")
  input <- prepare_input(y, a, x, outcome, zero_threshold, na_action)
  d <- input$data
  engine <- prepare_engine(d)
  reasons <- c("Requires at least four zero and eight positive subjects and an estimable direction.",
               "Requires at least eight positive subjects and an estimable rank direction.",
               "No evaluable step split; see scan results.")
  channels <- data.frame(channel = zeta_channels(), label = zeta_labels(),
                         statistic = NA_real_, standardized = NA_real_, p_value = NA_real_,
                         evaluable = FALSE, reason = reasons)
  diagnostics <- list(n_input = input$n_input, n_used = nrow(d),
                      omitted_rows = input$omitted, n_zero = sum(d$x == 0),
                      n_positive = sum(d$x > 0), n_recode = sum(d$x != d$x_original),
                      null_rank = engine$rank, residual_df = engine$df,
                      basis_names = engine$basis_names, outcome_offset = engine$outcome_offset,
                      scale_fallback = character())
  result <- list(call = call, p_global = NA_real_, status = "not_evaluable", reason = "",
                 channels = channels, scan = engine$scan, combinations = NULL,
                 settings = list(outcome = outcome, B = as.integer(B), seed = seed,
                                 zero_threshold = zero_threshold, scan_cuts = zeta_cuts(),
                                 ridge = 1e-8, na_action = na_action,
                                 rng_kind = RNGkind(), package_version = "0.1.0"),
                 diagnostics = diagnostics, data = d, bootstrap = NULL)
  class(result) <- "zeta_test"
  if (!length(engine$directions)) {
    result$reason <- "No evaluable interaction direction."
  } else if (engine$df <= 0L) {
    result$reason <- "No residual degrees of freedom in the null model."
  } else if (length(unique(d$y)) < 2L ||
             sqrt(mean((engine$r / max(abs(d$y - mean(d$y))))^2)) <= 100 * .Machine$double.eps) {
    result$reason <- "No usable outcome variation after fitting the null model."
  }
  if (nzchar(result$reason)) {
    result$channels$reason <- result$reason
    result$scan$evaluable <- FALSE
    result$scan$reason <- result$reason
    warning("ZETA test not evaluable: ", result$reason, call. = FALSE)
    return(result)
  }
  if (B < 199L) warning("B < 199 gives a coarse p value grid; use only for quick checks.", call. = FALSE)
  restore <- local_seed(seed)
  on.exit(restore(), add = TRUE)
  obs <- channel_scores(engine$r, engine$directions)
  qb <- matrix(NA_real_, B, length(obs$channels),
               dimnames = list(NULL, names(obs$channels)))
  for (b in seq_len(B)) {
    signed <- engine$r * sample(c(-1, 1), nrow(d), replace = TRUE)
    rb <- as.numeric(qr.resid(engine$qr, signed))
    if (any(!is.finite(rb))) stop("Nonfinite bootstrap residuals.", call. = FALSE)
    qb[b, ] <- channel_scores(rb, engine$directions)$channels
  }
  combined <- combine_channels(obs$channels, qb)
  idx <- match(names(obs$channels), result$channels$channel)
  result$channels$statistic[idx] <- as.numeric(obs$channels)
  result$channels$standardized[idx] <- as.numeric(combined$w)
  result$channels$p_value[idx] <- combined$channel_p
  result$channels$evaluable[idx] <- TRUE
  result$channels$reason[idx] <- ""
  for (j in seq_len(nrow(result$scan))) {
    nm <- paste0("step", j)
    if (nm %in% names(obs$steps)) result$scan$statistic[j] <- obs$steps[[nm]]
  }
  result$p_global <- combined$p_global
  result$status <- "ok"
  result$combinations <- combined$combinations
  result$diagnostics$bootstrap_center <- combined$center
  result$diagnostics$bootstrap_scale <- combined$scale
  result$diagnostics$scale_fallback <- combined$scale_fallback
  result$diagnostics$p_min <- combined$p_min
  if (keep_bootstrap) result$bootstrap <- qb
  if (any(!result$channels$evaluable)) {
    warning("Some channels are not evaluable; the global test uses the remaining channels. See channels$reason.", call. = FALSE)
  }
  if (length(combined$scale_fallback)) {
    warning("Near-zero bootstrap standard deviation: scale set to 1 as in the research implementation. See diagnostics.", call. = FALSE)
  }
  result
}
