args <- commandArgs(trailingOnly = TRUE)
pkg <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
version <- read.dcf(file.path(pkg, "DESCRIPTION"))[[1, "Version"]]
out <- file.path(pkg, "inst/doc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- list.files(file.path(pkg, "man"), pattern = "\\.Rd$", full.names = TRUE)
topics <- sub("\\.Rd$", "", basename(files))
links <- setNames(paste0(topics, ".html"), topics)
links <- c(links, "plot.zeta_test" = "zeta_plots.html",
  "plot.zeta_interpretation" = "zeta_plots.html")
Sys.setenv(`_R_HELP_ENABLE_ENHANCED_HTML_` = "false")
for (i in seq_along(files)) {
  dest <- file.path(out, paste0(topics[i], ".html"))
  tools::Rd2HTML(files[i], out = dest, package = pkg, Links = links,
    stylesheet = "guide.css", texmath = "", standalone = TRUE)
  html <- readLines(dest, warn = FALSE)
  html <- html[!grepl('<meta name="(generator|author)"', html, ignore.case = TRUE)]
  html <- gsub("00Index.html", "index.html", html, fixed = TRUE)
  writeLines(html, dest)
}
guides <- c("zeta_continuous_guide", "zeta_binary_guide", "zeta_workflow_guide")
for (nm in guides) tools::Rd2ex(file.path(pkg, "man", paste0(nm, ".Rd")),
                               out = file.path(out, paste0(nm, ".R")))
writeLines(c(
  "html { background: #fff; color: #191919; font: 16px/1.55 Arial, sans-serif; }",
  "body { max-width: 960px; margin: 32px auto; padding: 0 24px; }",
  "h1 { font-size: 26px; } h2 { font-size: 22px; } h3 { font-size: 19px; }",
  "h1, h2, h3 { line-height: 1.25; letter-spacing: 0; }",
  "a { color: #365b9d; } code, pre { font-family: monospace; font-size: 14px; }",
  "pre { padding: 14px; background: #f3f5f6; white-space: pre-wrap; overflow-wrap: anywhere; }",
  "p, td, code { overflow-wrap: anywhere; } table { border-collapse: collapse; max-width: 100%; }",
  "td, th { padding: 6px 12px 6px 0; text-align: left; vertical-align: top; }",
  "hr { border: 0; border-top: 1px solid #c5cad3; }",
  "@media (max-width: 600px) { body { margin: 16px auto; padding: 0 14px; } }"
), file.path(out, "guide.css"))
writeLines(c(
  '<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">',
  '<meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>ZetaHTE documentation</title><link rel="stylesheet" href="guide.css"></head><body>',
  paste0('<h1>ZetaHTE</h1><p>Version ', version, '</p>'),
  '<p>One semicontinuous biomarker candidate, two treatment groups, and a continuous or binary outcome.</p>',
  '<h2>Start here</h2><ol>',
  '<li><a href="zeta_continuous_guide.html">Continuous outcomes</a> | <a href="zeta_continuous_guide.R">R script</a></li>',
  '<li><a href="zeta_binary_guide.html">Binary outcomes and a user cutoff</a> | <a href="zeta_binary_guide.R">R script</a></li>',
  '<li><a href="zeta_workflow_guide.html">Input checks and multiple candidates</a> | <a href="zeta_workflow_guide.R">R script</a></li></ol>',
  '<p>Scripts run after installing and loading ZetaHTE. Examples use synthetic data only. Plot scripts use the active graphics device.</p>',
  '<h2>Reference</h2><ul>',
  '<li><a href="zeta_test.html">zeta_test: global test</a></li>',
  '<li><a href="zeta_interpret.html">zeta_interpret: exploratory summaries</a></li>',
  '<li><a href="zeta_methods.html">Printing and summaries</a></li>',
  '<li><a href="zeta_plots.html">Plot methods</a></li>',
  '<li><a href="zeta_reproduction.html">Reproducible analyses and scope</a></li></ul>',
  '<h2>Scope</h2><p>Channel results and effect summaries are exploratory. External real-data calibration, multiple-candidate adjustment, and patient bootstrap are separate workflows.</p>',
  '<p>Examples are synthetic. No clinical subject records are included.</p>',
  '</body></html>'
), file.path(out, "index.html"))
cat("Built", length(files), "HTML reference pages and", length(guides), "executable guides.\n")
