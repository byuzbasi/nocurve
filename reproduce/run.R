args <- commandArgs(TRUE)
if (!length(args) %in% 2:3 || !args[1L] %in% c("simulation", "raman") ||
    (args[1L] == "raman" && length(args) != 2L)) {
  stop("Usage: Rscript reproduce/run.R simulation|raman NEW_DIRECTORY [SIMULATION_CASE]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- dirname(normalizePath(script, mustWork = TRUE))
out <- args[2L]
if (file.exists(out) || !dir.exists(dirname(out))) stop("Choose a new directory with an existing parent")

pins <- jsonlite::fromJSON(file.path(root, "dependency_pins.json"))
for (i in seq_len(nrow(pins))) {
  name <- pins$name[i]
  if (!requireNamespace(name, quietly = TRUE) ||
      as.character(packageVersion(name)) != pins$version[i]) {
    stop("Install the pinned dependency ", name, " ", pins$version[i])
  }
}
stopifnot(requireNamespace("nocurve", quietly = TRUE))
python <- Sys.which("python3")
if (!nzchar(python)) stop("Python 3 is required for the numerical certificate and simulation summaries")
source(file.path(root, "lib/methods.R"), local = TRUE)
for (name in c("ms_sigma", "noc_design", "noc_base", "noc_output", "noc_validate_output",
               "noc_assert", "niv_check_base", "niv_check_native", "niv_pm", "n3a_impulse",
               "n3a_level", "ncp_descriptor", "ncp_functions", "ncp_rewrite", "ncp_analysis",
               "noc_inverse", "noc_details", "noc_wd")) {
  assign(name, getFromNamespace(name, "nocurve"))
}

# The spline's numerical controls are identical for the two study parts.
cp_cfg <- list(ss = list(low = -1.5, high = 1.5, tol = 1e-4, eps = 2e-8, maxit = 500L),
               ss_prediction_tol = 1e-8, ss_boundary_tol = 1e-3,
               ss_stable = list(probe_spar = 0.5, search_points = 129L, matrix_tol = 1e-10))
nr_raw <- function(...) invisible(NULL)
ap4_certificate <- function(sys, c, z) {
  input <- jsonlite::toJSON(ap4_payload(sys, c, z), auto_unbox = TRUE)
  text <- system2(python, c("-B", shQuote(file.path(root, "lib/certify_ap_residual_v04.py"))),
                  input = input, stdout = TRUE, stderr = TRUE)
  status <- attr(text, "status")
  stopifnot(is.null(status) || status == 0L)
  cert <- jsonlite::fromJSON(paste(text, collapse = "\n"))
  stopifnot(cert$valid, cert$normal_pass_exact, cert$readout_pass_exact,
            cert$normal_equation_residual_upper <= 1e-8, cert$tolerance == 1e-8,
            !cert$optimization_executed)
  cert
}
cache <- tempfile("nocurve-build-")
dir.create(cache)
Rcpp::sourceCpp(file.path(root, "lib/ap_strong_penalty_v04.cpp"), cacheDir = cache)
if (args[1L] == "raman") {
  Rcpp::sourceCpp(file.path(root, "lib/noc_three_real_kernel_v05.cpp"), cacheDir = cache)
}

# The simulation uses the manuscript's amended REML acceptance gate.
if (args[1L] == "simulation") {
  cp_gam_check <- function(f, d, k, method) mc3_check(f, d, k, method, allow_solve = TRUE)
}
spline_fit <- function(d, method) {
  z <- cp_gam(d, 96L, method)
  if (args[1L] == "simulation") {
    mc3_validate_anchor(z, d, 96L, method)
  } else {
    ncv6_validate_anchor(z, d, 96L, method)
  }
  z
}

raman_geometry <- function(y) {
  impulse <- n3a_impulse(length(y))
  levels <- lapply(0:(impulse$J - 1L), function(j) n3a_level(y, impulse, j))
  design <- ncp_descriptor(y, impulse, levels)
  noc_assert(design$analysis, noc_details(noc_wd(y[design$index])))
  design
}
raman_engine <- function() {
  f <- ncp_functions()
  f$noc_base <- ncp_rewrite(f$noc_base,
    list(list(quote(c(y, rev(y))), quote(y[D$index]))), f)
  f$noc_inverse <- ncp_rewrite(noc_inverse,
    list(list(quote(N <- 2L * n), quote(N <- length(theta) %/% w$nlevels))), f)
  f
}
fit_family <- function(d, D, engine) {
  base <- engine$noc_base(d$y, d$sigma_hat, D)
  engine$niv_check_base(base, D)
  lapply(setNames(c("NOC", "RAW", "EB0"), c("NOC", "RAW", "EB0")), function(m) {
    z <- engine$noc_output(base, D, m)
    engine$noc_validate_output(z, base, D, m)
    z[c("prediction", "omega")]
  })
}
fit_study <- function(d, design, engine, raman = FALSE) {
  family <- fit_family(d, design, engine)
  if (raman) {
    smash <- nsm_smash(d$y[design$pad_index], d$sigma_hat)
    smash$prediction <- smash$prediction[seq_len(design$n)]
    check_env <- new.env(parent = globalenv())
    check_env$ncp_synthesis <- ncp_functions()$ncp_synthesis
    check <- ncp_rewrite(niv_check_smash, list(
      list(quote(drop(D$A %*% y)), quote(ncp_analysis(D, y))),
      list(quote(drop(D$synthesis %*% theta + D$scaling %*% y)),
           quote(ncp_synthesis(D, theta, y)))), check_env)
    check(smash, d$y, d$sigma_hat, design)
    ss <- cvp_ss(d, data.frame(candidate = "native_GCV", value = NA_real_))
    kernel <- n3r5_kernel(d, data.frame(candidate = "native_dpill", value = NA_real_))
    ss <- list(prediction = ss$prediction[, 1L], spar = ss$table$value)
    kernel <- list(prediction = kernel$prediction[, 1L], bandwidth = kernel$table$value)
  } else {
    smash <- niv_smash(d$y, d$sigma_hat)
    niv_check_smash(smash, d$y, d$sigma_hat, design)
    ss <- pc_ss_stable_fit(d, d$x, cp_cfg)
    ss$spar <- ss$solution$spar
    kernel <- nr_kernel(d, nr_dpill(d))
  }
  ps <- spline_fit(d, "PS")
  ap <- spline_fit(d, "AP")
  prediction <- do.call(cbind, c(lapply(family, `[[`, "prediction"),
    list(SMASH4 = smash$prediction, PS = ps$prediction, AP = ap$prediction,
         SS = ss$prediction, KERNEL = kernel$prediction)))
  stopifnot(nrow(prediction) == length(d$y), ncol(prediction) == 8L, all(is.finite(prediction)))
  tuning <- c(lapply(family, `[[`, "omega"),
              list(SMASH4 = NA_real_, PS = ps$fit$sp, AP = ap$fit$sp,
                   SS = ss$spar, KERNEL = kernel$bandwidth))
  list(prediction = prediction, tuning = tuning)
}
tuning_table <- function(id, fit) {
  data.frame(dataset = id, method = names(fit$tuning),
    selected_parameter = vapply(fit$tuning, function(x) paste(format(x, digits = 17),
      collapse = ";"), ""))
}

if (args[1L] == "simulation") {
  spec <- jsonlite::fromJSON(file.path(root, "data/simulation.json"))
  grid <- expand.grid(curve = spec$curves, replication = seq_len(spec$replications_per_curve),
                      stringsAsFactors = FALSE)
  grid$id <- paste0("NMC100_", grid$curve, "_R", grid$replication)
  if (length(args) == 3L) {
    grid <- grid[grid$id == args[3L], ]
    if (nrow(grid) != 1L) stop("Unknown manuscript simulation case")
  }
  stopifnot(dir.create(out))
  geometry <- noc_design(spec$n)
  x <- (seq_len(spec$n) - 0.5) / spec$n
  output <- lapply(seq_len(nrow(grid)), function(i) {
    case <- grid[i, ]
    truth <- nt_eval_mean(spec$mean_functions[[case$curve]], x)
    seed <- spec$seed_base + 4L * (case$replication - 1L) + match(case$curve, spec$curves)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    y <- truth + spec$generating_sigma * rnorm(spec$n)
    d <- list(id = case$id, x = x, y = y, sigma_hat = ms_sigma(y))
    fit <- fit_study(d, geometry, globalenv())
    message("Completed ", case$id)
    list(curve = data.frame(dataset = case$id, curve = case$curve, replication = case$replication,
                            x = x, y = y, truth = truth, fit$prediction),
         tuning = tuning_table(case$id, fit))
  })
  write.csv(do.call(rbind, lapply(output, `[[`, "curve")),
            file.path(out, "predictions.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(output, `[[`, "tuning")),
            file.path(out, "tuning.csv"), row.names = FALSE)
  if (length(args) == 2L) {
    status <- system2(python, c("-B", shQuote(file.path(root, "lib/simulation_metrics.py")),
      "--directory", shQuote(normalizePath(out)), "--expected-cases", "400", "--repetitions", "100"))
    stopifnot(status == 0L)
  }
} else {
  data <- file.path(root, "data/raman")
  split <- jsonlite::fromJSON(file.path(data, "split.json"))
  read_split <- function(name) {
    path <- file.path(data, name)
    stopifnot(identical(digest::digest(file = path, algo = "sha256"), split$sha256[[name]]))
    read.delim(path)
  }
  stopifnot(split$source_row_zero_based == 440, split$training_acquisition_one_based == 1481,
            identical(split$test_reference_acquisitions_one_based, 1485:1500))
  train <- read_split("train.tsv")
  stopifnot(identical(train$channel, seq_len(split$channels)))
  d <- list(id = "Raman_440", x = train$x, y = train$observed, sigma_hat = ms_sigma(train$observed))
  stopifnot(dir.create(out))
  fit <- fit_study(d, raman_geometry(d$y), raman_engine(), raman = TRUE)
  # Read the test reference only after fitting and tuning all eight methods.
  test <- read_split("test.tsv")
  stopifnot(identical(train[c("channel", "x", "raman_shift")],
                      test[c("channel", "x", "raman_shift")]))
  errors <- abs(fit$prediction - test$reference_mean)
  metrics <- data.frame(method = colnames(errors), RMSE = sqrt(colMeans(errors^2)),
                        MAE = colMeans(errors), MAXAE = apply(errors, 2L, max))
  write.csv(cbind(train, reference_mean = test$reference_mean, fit$prediction),
            file.path(out, "predictions.csv"), row.names = FALSE)
  write.csv(tuning_table(d$id, fit), file.path(out, "tuning.csv"), row.names = FALSE)
  write.csv(metrics, file.path(out, "test_errors.csv"), row.names = FALSE)
}
