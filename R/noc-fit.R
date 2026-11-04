noc_fit <- function(y, method = c("EB0", "RAW", "NOC")) {
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) < 8L || any(!is.finite(y))) {
    stop("y must be a finite numeric vector with at least eight observations")
  }
  method <- match.arg(method, choices = noc_methods, several.ok = TRUE)
  y <- as.numeric(y)
  n <- length(y)
  sigma <- ms_sigma(y)
  if (sigma == 0) {
    if (!all(y == y[1L])) stop("Zero estimated noise on nonconstant data is unsupported")
    return(list(prediction = matrix(y, n, length(method), dimnames = list(NULL, method)),
                strength = setNames(numeric(length(method)), method), sigma = sigma))
  }

  # Preserve any pre-existing value in the historical ASH trace slot.
  trace_slot <- ".noc_native"
  had_slot <- exists(trace_slot, envir = .GlobalEnv, inherits = FALSE)
  if (had_slot) saved_slot <- get(trace_slot, envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_slot) {
      assign(trace_slot, saved_slot, envir = .GlobalEnv)
    } else if (exists(trace_slot, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = trace_slot, envir = .GlobalEnv)
    }
  })

  if (n <= 256L && n == 2^floor(log2(n))) {
    design <- noc_design(n)
    engine <- environment(noc_fit)
  } else {
    impulse <- n3a_impulse(n)
    levels <- lapply(seq_len(impulse$J) - 1L, function(j) n3a_level(y, impulse, j))
    design <- ncp_descriptor(y, impulse, levels)
    noc_assert(design$analysis, noc_details(noc_wd(y[design$index])))
    engine <- ncp_functions()
    engine$noc_base <- ncp_rewrite(engine$noc_base,
      list(list(quote(c(y, rev(y))), quote(y[D$index]))), engine)
    engine$noc_inverse <- ncp_rewrite(noc_inverse,
      list(list(quote(N <- 2L * n), quote(N <- length(theta) %/% w$nlevels))), engine)
  }
  base <- engine$noc_base(y, sigma, design)
  engine$niv_check_base(base, design)
  fits <- lapply(method, function(member) {
    result <- engine$noc_output(base, design, member)
    engine$noc_validate_output(result, base, design, member)
    result[c("prediction", "omega")]
  })
  list(prediction = vapply(setNames(fits, method), `[[`, numeric(n), "prediction"),
       strength = vapply(setNames(fits, method), `[[`, numeric(1L), "omega"),
       sigma = sigma)
}
