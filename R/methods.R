format_zeta_number <- function(x, digits = 4L) {
  vapply(x, function(v) if (is.na(v)) "NA" else format(signif(v, digits), trim = TRUE), character(1))
}

print_zeta_table <- function(x, digits) {
  for (nm in names(x)) if (is.numeric(x[[nm]])) x[[nm]] <- format_zeta_number(x[[nm]], digits)
  print(x, row.names = FALSE, right = FALSE)
}

print.zeta_test <- function(x, digits = 4L, ...) {
  cat("ZETA test\n")
  cat("Subjects:", x$diagnostics$n_used, " | Zero:", x$diagnostics$n_zero,
      " | Positive:", x$diagnostics$n_positive, " | Bootstrap:", x$settings$B, "\n")
  cat("Global p value:", format_zeta_number(x$p_global, digits), "\n")
  if (x$status != "ok") cat("Not evaluable:", x$reason, "\n")
  print_zeta_table(x$channels[, c("label", "standardized", "p_value", "evaluable")], digits)
  if (any(!x$channels$evaluable)) cat("See channels$reason for unavailable channels.\n")
  cat("Channel evidence is exploratory.\n")
  invisible(x)
}

summary.zeta_test <- function(object, ...) {
  object$data <- NULL
  object$bootstrap <- NULL
  object$call <- NULL
  class(object) <- "summary.zeta_test"
  object
}

print.summary.zeta_test <- function(x, digits = 4L, ...) {
  print.zeta_test(x, digits = digits)
  if (!is.null(x$combinations)) {
    cat("\nCombination diagnostics (not signal classifications):\n")
    print_zeta_table(x$combinations, digits)
  }
  invisible(x)
}

print.zeta_interpretation <- function(x, digits = 4L, ...) {
  cat("ZETA exploratory interpretation\n")
  cat("Global p value:", format_zeta_number(x$p_global, digits),
      " | Significant at", format_zeta_number(x$alpha), ":", x$global_significant, "\n")
  print_zeta_table(x$effects, digits)
  if (nrow(x$selected_scan)) cat("Exploratory rank location:", x$selected_scan$rank_cutoff, "\n")
  cat(x$note, "\n")
  invisible(x)
}
