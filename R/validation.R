check_scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " must be a finite ", if (integer) "integer " else "number ",
         "between ", lower, " and ", upper, ".", call. = FALSE)
  }
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

prepare_input <- function(y, a, x, outcome, zero_threshold, na_action) {
  for (nm in c("y", "a", "x")) {
    v <- get(nm)
    if (!(is.numeric(v) || (nm != "x" && is.logical(v))) ||
        !is.null(dim(v)) || is.object(v)) {
      stop(nm, " must be a plain numeric vector (0/1 for treatment).",
           call. = FALSE)
    }
  }
  if (length(y) != length(a) || length(y) != length(x)) {
    stop("y, a, and x must have the same length.", call. = FALSE)
  }
  if (any(is.infinite(y)) || any(is.infinite(a)) || any(is.infinite(x))) {
    stop("Infinite values are not supported, including on incomplete rows.", call. = FALSE)
  }
  missing <- is.na(y) | is.na(a) | is.na(x)
  if (any(missing) && na_action == "fail") {
    stop("Missing values found. Set na_action = 'omit' to remove incomplete rows.",
         call. = FALSE)
  }
  d <- data.frame(y = as.numeric(y), a = as.numeric(a), x = as.numeric(x),
                  row_id = seq_along(y))
  d <- d[!missing, , drop = FALSE]
  rownames(d) <- NULL
  if (nrow(d) < 4L) stop("At least four complete subjects are required.", call. = FALSE)
  if (any(!is.finite(as.matrix(d)))) {
    stop("Infinite values are not supported.", call. = FALSE)
  }
  if (any(d$x < 0)) stop("x must be nonnegative; recode the biomarker explicitly first.", call. = FALSE)
  if (!all(d$a %in% c(0, 1)) || any(tabulate(d$a + 1L, nbins = 2L) < 2L)) {
    stop("a must contain 0 and 1 with at least two complete subjects in each group.", call. = FALSE)
  }
  if (outcome == "binary" && !all(d$y %in% c(0, 1))) {
    stop("For outcome = 'binary', y must be coded 0/1.", call. = FALSE)
  }
  d$x_original <- d$x
  d$x[d$x <= zero_threshold] <- 0
  d$u <- 0
  pos <- d$x > 0
  d$u[pos] <- rank(d$x[pos], ties.method = "average") / (sum(pos) + 1)
  if (any(missing)) warning(sum(missing), " incomplete rows omitted.", call. = FALSE)
  list(data = d, omitted = which(missing), n_input = length(y))
}

local_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  if (RNGkind()[2L] == "Box-Muller" || any(RNGkind() == "user-supplied")) {
    stop("An explicit seed cannot safely restore this RNG's hidden state. ",
         "Use set.seed() before the call and seed = NULL, or a standard RNG ",
         "with the Inversion normal generator.", call. = FALSE)
  }
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (present) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (present) {
      assign(".Random.seed", previous, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}
