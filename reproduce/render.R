# Render the manuscript's empirical tables and figures from completed results.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript reproduce/render.R SIMULATION_DIR RAMAN_DIR NEW_OUTPUT_DIR",
       call. = FALSE)
}
simulation <- normalizePath(args[1L], mustWork = TRUE)
raman <- normalizePath(args[2L], mustWork = TRUE)
output <- args[3L]
if (file.exists(output) || !dir.exists(dirname(output))) {
  stop("Choose a new output directory with an existing parent.", call. = FALSE)
}
for (package in c("knitr", "kableExtra")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Rendering requires the R package ", package, ".", call. = FALSE)
  }
}

methods <- c("NOC", "RAW", "EB0", "SMASH4", "PS", "AP", "SS", "KERNEL")
labels <- c("NOC", "RAW", "EB0", "SMASH4", "P-spline", "Adaptive P-spline",
            "Smoothing spline", "Kernel")
colors <- c("#D62728", "#228833", "#A66A00", "#7B3294", "#00A6D6",
            "#D45087", "#233B8B", "#7A4E2D")
line_types <- c("solid", "44", "13", "4313", "2262", "3313", "122212", "51")
scenarios <- paste0("C", 1:4)
losses <- c("MSE", "MAE", "MAXAE")
test_losses <- c("RMSE", "MAE", "MAXAE")
simulation_columns <- c(losses, paste0("MCSE_", losses))

read_results <- function(path, columns, rows, comma = FALSE) {
  reader <- if (comma) read.csv else read.delim
  result <- reader(path, check.names = FALSE)
  stopifnot(nrow(result) == rows, !anyDuplicated(names(result)),
            all(columns %in% names(result)), !anyNA(result[columns]))
  result
}

ordered_methods <- function(block) {
  stopifnot(length(block$method) == length(methods), setequal(block$method, methods))
  block[match(methods, block$method), , drop = FALSE]
}

scenario_scores <- read_results(file.path(simulation, "scenario_metrics.tsv"),
  c("curve", "method", "R", simulation_columns), 32L)
overall_scores <- read_results(file.path(simulation, "balanced_metrics.tsv"),
  c("method", "R_per_curve", "fixed_curves", simulation_columns), 8L)
mean_curves <- read_results(file.path(simulation, "mean_curves.tsv"),
  c("curve", "method", "R", "x", "truth", "prediction"), 8192L)
stopifnot(all(scenario_scores$R == 100L), all(overall_scores$R_per_curve == 100L),
          all(overall_scores$fixed_curves == 4L), all(mean_curves$R == 100L),
          all(is.finite(as.matrix(scenario_scores[simulation_columns]))),
          all(is.finite(as.matrix(overall_scores[simulation_columns]))),
          all(is.finite(as.matrix(mean_curves[c("x", "truth", "prediction")]))))
panels <- lapply(scenarios, function(scenario) {
  block <- mean_curves[mean_curves$curve == scenario, ]
  stopifnot(nrow(block) == 256L * length(methods))
  baseline <- block[block$method == "NOC", ]
  stopifnot(identical(baseline$x, (seq_len(256L) - 0.5) / 256))
  estimates <- vapply(methods, function(method) {
    curve <- block[block$method == method, ]
    stopifnot(identical(curve$x, baseline$x), identical(curve$truth, baseline$truth))
    curve$prediction
  }, numeric(256L))
  list(x = baseline$x, truth = baseline$truth, estimates = estimates, title = paste("Scenario", scenario))
})

spectra <- read_results(file.path(raman, "predictions.csv"),
  c("channel", "x", "raman_shift", "observed", "reference_mean", methods), 800L, comma = TRUE)
raman_scores <- ordered_methods(read_results(file.path(raman, "test_errors.csv"),
  c("method", test_losses), 8L, comma = TRUE))
stopifnot(all(is.finite(as.matrix(spectra[c("x", "raman_shift", "observed", "reference_mean", methods)]))),
          identical(spectra$channel, seq_len(800L)), all(diff(spectra$raman_shift) > 0))
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
data_dir <- file.path(dirname(normalizePath(script, mustWork = TRUE)), "data", "raman")
for (name in c("train", "test")) {
  split <- read.delim(file.path(data_dir, paste0(name, ".tsv")))
  stopifnot(isTRUE(all.equal(spectra[names(split)], split, tolerance = 0)))
}
errors <- abs(as.matrix(spectra[methods]) - spectra$reference_mean)
expected <- cbind(RMSE = sqrt(colMeans(errors^2)), MAE = colMeans(errors),
                  MAXAE = apply(errors, 2L, max))
stopifnot(isTRUE(all.equal(unname(as.matrix(raman_scores[test_losses])), unname(expected),
                          tolerance = 1e-10)))

simulation_table <- function(block, group, heading) {
  out <- data.frame(Group = group, Method = labels, check.names = FALSE)
  names(out)[1L] <- heading
  for (loss in losses) {
    scale <- if (loss == "MSE") 1000 else 1
    digits <- if (loss == "MSE") 4L else 5L
    value <- block[[loss]]
    number <- formatC(scale * value, digits = digits, format = "f", decimal.mark = ".")
    best <- value == min(value)
    number[best] <- paste0("\\textbf{", number[best], "}")
    se <- formatC(scale * block[[paste0("MCSE_", loss)]], digits = digits,
                  format = "f", decimal.mark = ".")
    out[[if (loss == "MSE") "MSE $\\times10^3$" else loss]] <- paste0(number, " (", se, ")")
  }
  out
}
by_scenario <- do.call(rbind, lapply(scenarios, function(scenario) {
  block <- ordered_methods(scenario_scores[scenario_scores$curve == scenario, ])
  simulation_table(block, scenario, "Scenario")
}))
overall <- simulation_table(ordered_methods(overall_scores),
  rep(c("Matched family", "Published reference"), c(3L, 5L)), "Group")
raman_table <- data.frame(Method = labels)
for (loss in test_losses) {
  value <- raman_scores[[loss]]
  number <- formatC(value, digits = 4L, format = "f", decimal.mark = ".")
  best <- value == min(value)
  number[best] <- paste0("\\textbf{", number[best], "}")
  raman_table[[loss]] <- number
}

save_table <- function(data, name, alignment, grouped = FALSE) {
  latex <- knitr::kable(data, format = "latex", booktabs = TRUE, escape = FALSE,
                       row.names = FALSE, linesep = "", align = alignment)
  if (grouped) {
    latex <- kableExtra::collapse_rows(latex, columns = 1, valign = "top", latex_hline = "major")
  } else {
    latex <- sub("^\n", "", as.character(latex))
    latex <- gsub(" +(?=\\\\\\\\)", "", latex, perl = TRUE)
  }
  writeLines(as.character(latex), file.path(output, paste0(name, ".tex")))
}

draw_legend <- function(first_label, first_color, first_width) {
  par(mar = rep(0, 4), cex = 1)
  plot.new()
  keys <- c(first_label, labels)
  settings <- list(x = "center", legend = keys, col = c(first_color, colors),
    lty = c("solid", line_types), lwd = c(first_width, rep(1.05, 8)),
    horiz = TRUE, bty = "n", cex = .78, seg.len = 1.2, x.intersp = .3,
    text.width = strwidth(keys, cex = .78))
  test <- do.call(legend, c(settings, list(plot = FALSE)))
  bounds <- par("usr")
  stopifnot(length(unique(test$text$y)) == 1L, test$rect$left >= bounds[1],
            test$rect$left + test$rect$w <= bounds[2], test$rect$top <= bounds[4],
            test$rect$top - test$rect$h >= bounds[3])
  do.call(legend, settings)
}

draw_estimates <- function(x, estimates) {
  for (j in seq_along(methods)) {
    lines(x, estimates[, j], col = colors[j], lty = line_types[j], lwd = 1.05)
  }
}

save_figure <- function(name, height, draw) {
  postscript(file.path(output, paste0(name, ".eps")), onefile = FALSE,
             horizontal = FALSE, paper = "special", width = 7, height = height,
             pointsize = 10, family = "Helvetica")
  on.exit(dev.off())
  draw()
}

stopifnot(dir.create(output))
save_table(by_scenario, "scenario_losses", "llrrr", grouped = TRUE)
save_table(overall, "overall_losses", "llrrr", grouped = TRUE)
save_table(raman_table, "raman_case_metrics", "lrrr")
save_figure("simulation_mean_curves", 6.7, function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, byrow = TRUE), heights = c(1, 1, .14))
  par(mar = c(3.3, 3.7, 2.1, .7), mgp = c(2, .6, 0), tcl = -.25, cex.axis = .85, cex.lab = .9)
  for (panel in panels) {
    plot(panel$x, panel$truth, type = "n", xlab = "Design coordinate", ylab = "Response",
         main = panel$title, cex.main = .95, ylim = range(c(panel$truth, panel$estimates)))
    lines(panel$x, panel$truth, col = "#8D99A6", lwd = .8)
    draw_estimates(panel$x, panel$estimates)
  }
  draw_legend("True mean", "#8D99A6", .8)
})
save_figure("raman_case", 4.1, function() {
  layout(matrix(1:2, ncol = 1), heights = c(1, .15))
  par(mar = c(3.3, 3.7, .8, .7), mgp = c(2, .6, 0), tcl = -.25, cex.axis = .85, cex.lab = .9)
  estimates <- as.matrix(spectra[methods])
  plot(spectra$raman_shift, spectra$observed, type = "n", xlab = "Stored Raman coordinate",
       ylab = "Intensity", ylim = range(c(spectra$observed, estimates)))
  lines(spectra$raman_shift, spectra$observed, col = "#EBEBEB", lwd = .55)
  draw_estimates(spectra$raman_shift, estimates)
  draw_legend("Input", "#EBEBEB", .55)
})
message("Wrote three LaTeX tables and two EPS figures to ", output, ".")
