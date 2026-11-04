library(nocurve)
set.seed(51)
x <- (seq_len(32) - 0.5) / 32
y <- sin(2 * pi * x) + rnorm(length(x), sd = 0.15)
z <- noc_fit(y)
stopifnot(identical(dim(z$prediction), c(32L, 3L)), all(is.finite(z$prediction)),
          all(z$strength >= 0 & z$strength <= 1), z$strength[["EB0"]] == 0)
one <- noc_fit(rep(2, 29L), "NOC")
stopifnot(identical(as.numeric(one$prediction), rep(2, 29L)), one$sigma == 0)
bad <- try(noc_fit(c(y, NA_real_)), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
