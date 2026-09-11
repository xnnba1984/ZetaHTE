check_zeta_plot <- function(main, cex, colors, dots) {
  if (!is.null(main) && (!is.character(main) || length(main) != 1L || is.na(main))) {
    stop("main must be NULL or one character string.", call. = FALSE)
  }
  check_scalar(cex, "cex", 0)
  if (cex <= 0) stop("cex must be positive.", call. = FALSE)
  if (length(dots)) stop("Unused plot arguments; use main, cex, or colors.", call. = FALSE)
  if (!is.character(colors) || anyNA(colors)) stop("colors must contain color names or hex codes.", call. = FALSE)
  tryCatch(grDevices::col2rgb(colors), error = function(e) {
    stop("Invalid plot color.", call. = FALSE)
  })
}

zeta_plot_frame <- function(labels, limits, title, subtitle, xlab, cex, annotation) {
  graphics::par(cex = cex, mar = c(4.6, 2, 3.7, 0.8), mgp = c(2.2, 0.6, 0), las = 1)
  left <- max(graphics::strwidth(labels, units = "inches")) + 0.2
  graphics::par(mai = c(1.03, left, 0.72, 0.14))
  pin <- graphics::par("pin")
  ann <- max(graphics::strwidth(annotation, units = "inches")) + 0.2
  if (pin[1] - ann < 1.5 || pin[2] < max(0.9, length(labels) * 0.27)) {
    stop("Plot panel is too small. Enlarge the device or reduce cex.", call. = FALSE)
  }
  right <- 1 - ann / pin[1]
  map <- function(v) 0.035 + (v - limits[1]) / diff(limits) * (right - 0.08)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0.45, length(labels) + 0.55),
                        xaxs = "i", yaxs = "i")
  yy <- rev(seq_along(labels))
  graphics::axis(2, at = yy, labels = labels, tick = FALSE, line = -0.3)
  ticks <- pretty(limits, n = 4)
  ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  graphics::axis(1, at = map(ticks), labels = format_zeta_number(ticks),
                 col = "#4B4B4B", col.axis = "#000000")
  graphics::mtext(xlab, side = 1, line = 2.1, at = right / 2, cex = cex)
  graphics::title(main = title, cex.main = 1.05, adj = 0, line = 2)
  graphics::mtext(subtitle, side = 3, adj = 0, line = 0.6, cex = 0.82 * cex)
  graphics::segments(map(limits[1]), yy, map(limits[2]), yy, col = "#E6E8F0")
  list(y = yy, map = map, annotation_x = right + 0.03)
}

zeta_plot_limits <- function(values, include_zero = TRUE) {
  values <- values[is.finite(values)]
  if (include_zero) values <- c(0, values)
  if (!length(values)) return(c(0, 1))
  lim <- range(values)
  if (diff(lim) == 0) lim <- lim + c(-0.5, 0.5)
  lim
}

plot.zeta_test <- function(x, type = c("channels", "scan"), alpha = 0.05,
                           main = NULL, cex = 1, colors = "#5477C4", ...) {
  type <- match.arg(type)
  check_zeta_plot(main, cex, colors, list(...))
  if (length(colors) != 1L) stop("Supply one channel color.", call. = FALSE)
  check_scalar(alpha, "alpha", 0, 1)
  if (alpha <= 0 || alpha >= 1) stop("alpha must be strictly between zero and one.", call. = FALSE)
  old <- graphics::par(c("cex", "mar", "mgp", "las"))
  on.exit(graphics::par(old), add = TRUE)
  if (type == "channels") {
    d <- x$channels
    d$display_value <- -log10(d$p_value)
    d$annotation <- ifelse(d$evaluable, paste0("p = ", format_zeta_number(d$p_value)), "Not evaluable")
    lim <- c(0, max(c(1, -log10(alpha), d$display_value[is.finite(d$display_value)])))
    f <- zeta_plot_frame(d$label, lim, if (is.null(main)) "Channel evidence" else main,
      paste0("Global p = ", format_zeta_number(x$p_global), "  |  Exploratory channel p values"),
      expression(-log[10](p)), cex, d$annotation)
    graphics::segments(f$map(-log10(alpha)), 0.65, f$map(-log10(alpha)), 3.35,
                       col = "#7A7A7A", lty = 2)
    graphics::mtext(paste0("Dashed line: p = ", format_zeta_number(alpha), " (exploratory reference)"),
                    side = 1, line = 3.6, cex = 0.72 * cex, adj = 0)
  } else {
    d <- x$scan
    d$display_value <- d$statistic
    d$annotation <- ifelse(d$evaluable, paste0("Q = ", format_zeta_number(d$statistic)), "Not evaluable")
    f <- zeta_plot_frame(paste("Rank", format_zeta_number(d$rank_cutoff)),
      c(0, max(c(1, d$statistic[is.finite(d$statistic)]))),
      if (is.null(main)) "Tail step scan" else main,
      "Three prespecified locations in the positive tail", "Q statistic", cex, d$annotation)
    graphics::mtext("Locations describe coarse exploratory splits.", side = 1, line = 3.6,
                    cex = 0.72 * cex, adj = 0)
  }
  valid <- d$evaluable & is.finite(d$display_value)
  graphics::points(f$map(d$display_value[valid]), f$y[valid], pch = 21,
                   bg = colors, col = "#4B4B4B", cex = 1.35)
  graphics::text(f$annotation_x, f$y, labels = d$annotation, adj = 0, cex = 0.88)
  invisible(d)
}

plot.zeta_interpretation <- function(x, type = c("groups", "contrasts"),
                                     main = NULL, cex = 1,
                                     colors = c("#2A9D8F", "#D5A445"), ...) {
  type <- match.arg(type)
  check_zeta_plot(main, cex, colors, list(...))
  if (length(colors) != 2L) stop("Supply treatment and control colors, in that order.", call. = FALSE)
  old <- graphics::par(c("cex", "mar", "mgp", "las"))
  on.exit(graphics::par(old), add = TRUE)
  binary <- identical(x$effects$scale[1], "risk difference")
  mult <- if (binary) 100 else 1
  subtitle <- if (x$global_significant) "Exploratory summaries after a significant global test" else
    "Descriptive only: global test is not significant"
  if (type == "groups") {
    d <- x$groups
    labels <- c(zero = "Zero", positive = "All positive", step_lower = "Step: lower",
      step_upper = "Step: upper", user_lower = "User: lower", user_upper = "User: upper")
    if (nrow(x$selected_scan)) {
      cut <- format_zeta_number(x$selected_scan$rank_cutoff)
      labels[c("step_lower", "step_upper")] <- paste0("Step: rank ", c("<= ", "> "), cut)
    }
    if (!is.null(x$user_cutoff)) {
      cut <- format_zeta_number(x$user_cutoff)
      upper <- identical(x$user_cutoff_rule, "positive X < cutoff versus X >= cutoff")
      labels[c("user_lower", "user_upper")] <- c(
        paste0("User: 0 < X ", if (upper) "< " else "<= ", cut),
        paste0("User: X ", if (upper) ">= " else "> ", cut))
    }
    d$display_group <- unname(labels[d$group])
    d$treatment_display <- d$mean_treatment * mult
    d$control_display <- d$mean_control * mult
    annotate <- function(v, n) ifelse(is.finite(v),
      paste0(format_zeta_number(v), if (binary) "%" else "", " (n = ", n, ")"),
      paste0("NA (n = ", n, ")"))
    d$treatment_label <- annotate(d$treatment_display, d$n_treatment)
    d$control_label <- annotate(d$control_display, d$n_control)
    lim <- if (binary) c(0, 100) else zeta_plot_limits(c(d$treatment_display, d$control_display))
    f <- zeta_plot_frame(d$display_group, lim,
      if (is.null(main)) "Group summaries" else main, subtitle,
      if (binary) "Response rate (%)" else "Outcome mean", cex,
      c(d$treatment_label, d$control_label))
    both <- is.finite(d$treatment_display) & is.finite(d$control_display)
    graphics::segments(f$map(d$treatment_display[both]), f$y[both] + 0.13,
      f$map(d$control_display[both]), f$y[both] - 0.13, col = "#C5CAD3")
    for (k in 1:2) {
      val <- if (k == 1) d$treatment_display else d$control_display
      lab <- if (k == 1) d$treatment_label else d$control_label
      yy <- f$y + if (k == 1) 0.13 else -0.13
      ok <- is.finite(val)
      graphics::points(f$map(val[ok]), yy[ok], pch = if (k == 1) 21 else 22,
                       bg = colors[k], col = "#4B4B4B", cex = 1.2)
      graphics::text(f$annotation_x, yy, lab, adj = 0, col = colors[k], cex = 0.78)
    }
    graphics::mtext("Circle: treatment; square: control. Groups overlap.",
                    side = 1, line = 3.6, cex = 0.72 * cex, adj = 0)
  } else {
    d <- x$effects[x$effects$summary != "Tail linear slope", , drop = FALSE]
    d$display_value <- d$estimate * mult
    d$annotation <- ifelse(is.finite(d$estimate), format_zeta_number(d$display_value), "Not estimable")
    f <- zeta_plot_frame(d$summary, zeta_plot_limits(d$display_value),
      if (is.null(main)) "Effect contrasts" else main, subtitle,
      if (binary) "Response rate difference (%)" else "Treatment effect difference", cex, d$annotation)
    graphics::segments(f$map(0), 0.65, f$map(0), nrow(d) + 0.35, lty = 2, col = "#7A7A7A")
    ok <- is.finite(d$display_value)
    graphics::points(f$map(d$display_value[ok]), f$y[ok], pch = 21,
                     bg = "#7B8F45", col = "#4B4B4B", cex = 1.35)
    graphics::text(f$annotation_x, f$y, d$annotation, adj = 0, cex = 0.88)
    graphics::mtext("Jump: zero - positive; step/user: upper - lower.",
                    side = 1, line = 3.6, cex = 0.68 * cex, adj = 0)
  }
  invisible(d)
}
