library(ZetaHTE)
ncheck <- 0L
check <- function(x) {
  ncheck <<- ncheck + 1L
  if (!isTRUE(x)) stop("Plot check failed: ", ncheck)
}
error <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
set.seed(18)
x <- rep(c(0, seq(0.1, 1, length.out = 9)), each = 12)
a <- rep(0:1, length(x) / 2)
y <- rnorm(length(x), sd = 0.3) + a * (0.2 + 5 * x)
fit <- zeta_test(y, a, x, B = 199, seed = 42)
explanation <- zeta_interpret(fit, cutoff = 0.6)
file <- tempfile(fileext = ".pdf")
pdf(file, width = 8, height = 6, title = "", author = "")
before <- serialize(fit, NULL)
before_e <- serialize(explanation, NULL)
rng <- .Random.seed
gp <- par(c("cex", "mar", "mgp", "las"))
for (kind in c("channels", "scan")) {
  d <- plot(fit, type = kind)
  check(is.data.frame(d))
  check(nrow(d) == 3L)
  check(identical(par(names(gp)), gp))
  check(identical(.Random.seed, rng))
  check(identical(serialize(fit, NULL), before))
  expected <- if (kind == "channels") -log10(fit$channels$p_value) else fit$scan$statistic
  check(identical(d$display_value, expected))
}
for (kind in c("groups", "contrasts")) {
  d <- plot(explanation, type = kind)
  check(is.data.frame(d))
  check(identical(par(names(gp)), gp))
  check(identical(.Random.seed, rng))
  check(identical(serialize(explanation, NULL), before_e))
  if (kind == "groups") {
    check(identical(d$treatment_display, explanation$groups$mean_treatment))
    check(identical(d$control_display, explanation$groups$mean_control))
    check(nrow(d) == 6L)
  } else {
    check(!any(d$summary == "Tail linear slope"))
    check(identical(d$display_value, d$estimate))
    check(nrow(d) == 3L)
  }
}
set.seed(202)
xb <- rep(c(0, 2, 5, 8, 10, 12, 15, 20, 25), each = 40)
ab <- rep(0:1, length(xb) / 2)
yb <- rbinom(length(xb), 1, 0.15 + ab * (0.05 + 0.6 * (xb >= 10)))
fb <- zeta_test(yb, ab, xb, outcome = "binary", B = 199, seed = 73)
eb <- zeta_interpret(fb, cutoff = 10, cutoff_inclusive = "upper")
gb <- plot(eb, type = "groups")
cb <- plot(eb, type = "contrasts")
check(identical(gb$treatment_display, 100 * eb$groups$mean_treatment))
check(identical(gb$control_display, 100 * eb$groups$mean_control))
check(identical(cb$display_value, 100 * cb$estimate))
check(all(gb$treatment_display >= 0 & gb$treatment_display <= 100))
check(all(grepl("% \\(n = ", gb$treatment_label)))
check(all(grepl("% \\(n = ", gb$control_label)))
check(all(cb$scale == "risk difference"))
check(identical(gb$display_group[5:6], c("User: 0 < X < 10", "User: X >= 10")))
check(gb$display_group[3] == paste("Step: rank <=", eb$selected_scan$rank_cutoff))
check(gb$display_group[4] == paste("Step: rank >", eb$selected_scan$rank_cutoff))
lower <- plot(zeta_interpret(fb, cutoff = 10), type = "groups")
check(identical(lower$display_group[5:6], c("User: 0 < X <= 10", "User: X > 10")))
check(sum(gb$n_treatment[5:6]) == sum(lower$n_treatment[5:6]))

# Genuine input-generated unavailable channels and empty descriptive groups.
nozero <- suppressWarnings(zeta_test(y, a, x + 1, B = 199, seed = 42))
partial <- plot(nozero)
check(is.na(partial$display_value[1]))
check(partial$annotation[1] == "Not evaluable")
flat <- suppressWarnings(zeta_test(y, a, rep(0, length(x)), B = 199, seed = 42))
empty <- plot(flat)
check(all(is.na(empty$display_value)))
check(all(empty$annotation == "Not evaluable"))
empty_scan <- plot(flat, type = "scan")
check(all(empty_scan$annotation == "Not evaluable"))
outside <- suppressWarnings(zeta_interpret(fit, cutoff = 100))
eg <- plot(outside, type = "groups")
ec <- plot(outside, type = "contrasts")
check(eg$treatment_label[6] == "NA (n = 0)")
check(eg$control_label[6] == "NA (n = 0)")
check(is.na(ec$display_value[3]))
check(ec$annotation[3] == "Not estimable")

for (thunk in list(function() plot(fit, type = "wrong"),
                   function() plot(fit, alpha = 0),
                   function() plot(fit, alpha = 1),
                   function() plot(fit, cex = 0),
                   function() plot(fit, cex = Inf),
                   function() plot(fit, main = NA_character_),
                   function() plot(fit, colors = "invalid_color"),
                   function() plot(fit, colors = c("red", "blue")),
                   function() plot(explanation, colors = "red"),
                   function() plot(fit, xlim = c(0, 2)))) {
  check(error(thunk()))
  check(identical(par(names(gp)), gp))
}
check(is.data.frame(plot(fit, colors = "black", main = "Custom title", cex = 0.9)))
check(is.data.frame(plot(eb, type = "groups", colors = c("black", "gray40"), cex = 0.9)))
dev.off()
unlink(file)

tiny <- tempfile(fileext = ".pdf")
pdf(tiny, width = 2, height = 2, title = "", author = "")
small_gp <- par(c("cex", "mar", "mgp", "las"))
check(error(plot(fit)))
check(identical(par(names(small_gp)), small_gp))
dev.off()
unlink(tiny)
cat("PASS:", ncheck, "plot assertions.\n")
