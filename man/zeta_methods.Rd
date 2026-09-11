\name{zeta_methods}
\alias{print.zeta_test}
\alias{summary.zeta_test}
\alias{print.summary.zeta_test}
\alias{print.zeta_interpretation}
\title{Print and Summarize ZETA Results}
\description{Compact displays of the global test, channel evidence, and
descriptive summaries. Summary objects omit subject-level data and bootstrap matrices.}
\usage{
\method{print}{zeta_test}(x, digits = 4L, ...)
\method{summary}{zeta_test}(object, ...)
\method{print}{summary.zeta_test}(x, digits = 4L, ...)
\method{print}{zeta_interpretation}(x, digits = 4L, ...)
}
\arguments{
\item{x, object}{The corresponding ZETA result object.}
\item{digits}{Number of significant digits in printed numbers.}
\item{...}{Reserved for compatibility with the standard generic functions.}
}
\value{Print methods return their argument invisibly. The summary method returns
a compact object including parameters and diagnostic tables.}
