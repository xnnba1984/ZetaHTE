\name{zeta_reproduction}
\alias{zeta_reproduction}
\title{Reproducible Package Analyses and Scope}
\description{Record an analysis configuration and distinguish the package
workflow from additional study-level analyses.}
\section{Reproducing a test}{
Retain the input vectors and their preprocessing, including missing-value
handling and any zero threshold. Record the outcome type, treatment coding,
bootstrap count, seed, package version, R version, and random-number generators.
The test object stores the analysis settings and complete analysis rows.
Keep objects containing confidential subject data private.

Use the same inputs, settings, and software environment to reproduce a result.
An explicit seed preserves the caller's random-number state with supported
built-in generators. When using Box-Muller or user-supplied generators, call
\code{set.seed()} externally and use \code{seed = NULL}.
The internal bootstrap count controls p-value resolution; the default is 999.
Independent simulated datasets are needed to study type I error or power.
}
\section{Examples and tests}{
The three installed guides contain synthetic examples with fixed seeds.
They can be read as help pages or run from the matching scripts in
\code{system.file("doc", package = "ZetaHTE")}.
The source package includes automated tests for numerical behavior,
input validation, random-number handling, interpretation, and plotting.
}
\section{Scope}{
Each call tests one biomarker candidate with two treatment groups.
Channel results and effect summaries are exploratory. User cutoffs affect
descriptive summaries only. External calibration across shuffled datasets,
multiple-candidate adjustment, and patient bootstrap analyses are separate
workflows. The package does not run these analyses automatically.

This software release includes neither clinical subject records nor a complete
manuscript reproduction bundle. Synthetic examples demonstrate package use;
they are not clinical results or evidence of performance across all settings.
}
\seealso{\code{\link{zeta_continuous_guide}}, \code{\link{zeta_binary_guide}},
\code{\link{zeta_workflow_guide}}}
