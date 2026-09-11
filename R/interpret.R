group_effect <- function(d, index, name) {
  treated <- d$y[index & d$a == 1]
  control <- d$y[index & d$a == 0]
  valid <- length(treated) > 0L && length(control) > 0L
  mt <- if (length(treated)) mean(treated) else NA_real_
  mc <- if (length(control)) mean(control) else NA_real_
  data.frame(group = name, n_treatment = length(treated), n_control = length(control),
             mean_treatment = mt, mean_control = mc, effect = mt - mc,
             reason = if (valid) "" else "One treatment group has no subjects.")
}

zeta_interpret <- function(object, alpha = 0.05, cutoff = NULL,
                           cutoff_inclusive = c("lower", "upper"),
                           allow_nonsignificant = FALSE) {
  if (!inherits(object, "zeta_test")) stop("object must be returned by zeta_test().", call. = FALSE)
  check_scalar(alpha, "alpha", 0, 1)
  if (alpha <= 0 || alpha >= 1) stop("alpha must be strictly between zero and one.", call. = FALSE)
  check_flag(allow_nonsignificant, "allow_nonsignificant")
  cutoff_inclusive <- match.arg(cutoff_inclusive)
  if (!is.null(cutoff)) {
    check_scalar(cutoff, "cutoff", 0)
    if (cutoff <= 0) stop("cutoff must be strictly positive.", call. = FALSE)
  }
  if (!identical(object$status, "ok")) stop("The global test is not evaluable.", call. = FALSE)
  significant <- object$p_global <= alpha
  if (!significant && !allow_nonsignificant) {
    stop("Global p value exceeds alpha. Set allow_nonsignificant = TRUE for explicitly descriptive output.", call. = FALSE)
  }
  if (!significant) warning("Global test is not significant; these summaries do not indicate a discovery.", call. = FALSE)
  d <- object$data
  pos <- d$x > 0
  groups <- rbind(group_effect(d, !pos, "zero"), group_effect(d, pos, "positive"))
  scale <- if (object$settings$outcome == "binary") "risk difference" else "mean difference"
  effects <- data.frame(summary = c("Jump contrast", "Tail linear slope", "Tail step contrast"),
                        estimate = NA_real_, scale = c(scale, paste(scale, "per unit rank"), scale),
                        reason = "")
  effects$estimate[1] <- groups$effect[1] - groups$effect[2]
  if (!is.finite(effects$estimate[1])) effects$reason[1] <- "Zero or positive group lacks one treatment arm."
  tail <- d[pos, , drop = FALSE]
  if (nrow(tail) > 8L && length(unique(tail$a)) == 2L && stats::sd(tail$u) > 0) {
    xx <- stats::model.matrix(~ a * u, data = tail)
    fit <- stats::lm.fit(xx, tail$y)
    cc <- fit$coefficients["a:u"]
    if (is.finite(cc) && fit$rank == ncol(xx)) effects$estimate[2] <- unname(cc)
  }
  if (!is.finite(effects$estimate[2])) {
    effects$reason[2] <- "Requires more than eight positive subjects and an estimable treatment by rank slope."
  }
  scan <- object$scan
  selected <- scan[FALSE, , drop = FALSE]
  if (any(scan$evaluable)) {
    selected <- scan[which.max(replace(scan$statistic, !scan$evaluable, -Inf)), , drop = FALSE]
    low <- pos & d$u <= selected$rank_cutoff
    high <- pos & d$u > selected$rank_cutoff
    step <- rbind(group_effect(d, low, "step_lower"), group_effect(d, high, "step_upper"))
    effects$estimate[3] <- step$effect[2] - step$effect[1]
    groups <- rbind(groups, step)
    if (!is.finite(effects$estimate[3])) effects$reason[3] <- "One side lacks one treatment arm."
  } else {
    effects$reason[3] <- "No evaluable scan location."
  }
  user_rule <- NULL
  if (!is.null(cutoff)) {
    if (cutoff_inclusive == "lower") {
      low <- pos & d$x <= cutoff
      high <- pos & d$x > cutoff
      user_rule <- "positive X <= cutoff versus X > cutoff"
    } else {
      low <- pos & d$x < cutoff
      high <- pos & d$x >= cutoff
      user_rule <- "positive X < cutoff versus X >= cutoff"
    }
    user <- rbind(group_effect(d, low, "user_lower"), group_effect(d, high, "user_upper"))
    est <- user$effect[2] - user$effect[1]
    effects <- rbind(effects, data.frame(summary = "User cutoff contrast", estimate = est,
                                         scale = scale, reason = if (is.finite(est)) "" else "One side lacks one treatment arm."))
    groups <- rbind(groups, user)
  }
  if (any(!is.finite(effects$estimate))) warning("Some descriptive summaries are not estimable; see effects$reason.", call. = FALSE)
  structure(list(p_global = object$p_global, alpha = alpha, global_significant = significant,
                  channels = object$channels, effects = effects, groups = groups,
                  selected_scan = selected, user_cutoff = cutoff,
                  user_cutoff_rule = user_rule,
                  note = "Exploratory summaries; localized effects require independent confirmation."),
             class = "zeta_interpretation")
}
