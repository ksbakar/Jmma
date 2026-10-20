
################################################################################
# Print method for Jmma objects
################################################################################

print.Jmma <- function(x, digits = 3, ...) {

  if (!inherits(x, "Jmma")) {
    stop("Object must inherit from class 'Jmma'.")
  }
  settings <- x$settings
  cat("\nJoint Mediation Analysis (JMMA)\n")
  cat("----------------------------------------\n")
  cat("Outcome:             ", settings$outcome, "\n", sep = "")
  cat("Exposure:            ", settings$exposure, "\n", sep = "")
  cat("Reference category:  ", settings$reference, "\n", sep = "")
  cat("Analysis sample size: ", settings$analysis_sample_size, "\n", sep = "")
  cat("Mediators:           ",
      paste(settings$mediators, collapse = ", "), "\n", sep = "")
  cat("Covariates:          ",
      paste(settings$covariates, collapse = ", "), "\n", sep = "")
  if (isTRUE(settings$use_pseudo)) {
    cat("Outcome type:        Competing-risks RMST pseudo-values\n")
    cat("Restriction time:    ", settings$tau, "\n", sep = "")
    cat("Event code:          ", settings$event_code, "\n", sep = "")
  } else {
    cat("Outcome type:        Observed outcome\n")
  }
  cat("Bootstrap replicates requested: ",
      settings$n_boot, "\n", sep = "")
  cat("Confidence level:     ",
      100 * settings$conf_level, "%\n", sep = "")
  cat("Comparisons:          ",
      length(settings$comparisons), "\n", sep = "")
  invisible(x)
}


################################################################################
# Summary method for Jmma objects
################################################################################

summary.Jmma <- function(object, cov_matrix = FALSE, ...) {

  if (!inherits(object, "Jmma")) {
    stop("Object must inherit from class 'Jmma'.")
  }
  if (!is.logical(cov_matrix) ||
      length(cov_matrix) != 1L ||
      is.na(cov_matrix)) {
    stop("'cov_matrix' must be TRUE or FALSE.")
  }
  # Collect diagnostics for each comparison.
  diagnostics <- do.call(
    rbind,
    lapply(names(object$inference), function(comp) {
      inf <- object$inference[[comp]]
      dec <- inf$decomposition_diagnostic
      data.frame(
        Comparison = comp,
        Bootstrap_replicates = inf$n_boot,
        Sample_size = inf$n,
        Covariance_rank = inf$estimated_covariance_rank,
        TE_minus_NDE_minus_NIE = dec$point_estimate_error,
        Max_abs_bootstrap_decomposition_error =
          dec$max_absolute_bootstrap_error,
        row.names = NULL
      )
    })
  )
  rownames(diagnostics) <- NULL
  # Store covariance matrices by comparison.
  bootstrap_covariance <- lapply(
    object$inference,
    function(inf) inf$bootstrap_covariance
  )
  covariance_NDE_NIE <- lapply(
    object$inference,
    function(inf) inf$covariance_NDE_NIE
  )
  Sigma_hat <- lapply(
    object$inference,
    function(inf) inf$Sigma_hat
  )
  V_gamma_hat <- lapply(
    object$inference,
    function(inf) inf$V_gamma_hat
  )
  out <- list(
    settings = object$settings,
    results = object$results,
    joint_tests = object$joint_tests,
    diagnostics = diagnostics,
    bootstrap_covariance = bootstrap_covariance,
    covariance_NDE_NIE = covariance_NDE_NIE,
    Sigma_hat = Sigma_hat,
    V_gamma_hat = V_gamma_hat,
    inference = object$inference,
    cov_matrix = cov_matrix
  )
  class(out) <- "summary.Jmma"
  out
}


################################################################################
# Print method for summary.Jmma objects
################################################################################

print.summary.Jmma <- function(x, digits = 4, ...) {

  cat("\nSummary of Joint Mediation Analysis\n")
  cat("========================================\n")
  cat("\nModel settings\n")
  cat("----------------------------------------\n")
  cat("Exposure:             ", x$settings$exposure, "\n", sep = "")
  cat("Reference:            ", x$settings$reference, "\n", sep = "")
  cat("Outcome:              ", x$settings$outcome, "\n", sep = "")
  cat("Analysis sample size:  ",
      x$settings$analysis_sample_size, "\n", sep = "")
  cat("Confidence level:      ",
      100 * x$settings$conf_level, "%\n", sep = "")
  cat("\nMarginal effect estimates\n")
  cat("----------------------------------------\n")
  print(x$results, row.names = FALSE, digits = digits)
  cat("\nJoint Wald tests\n")
  cat("----------------------------------------\n")
  cat("Null hypothesis (H0): NDE = 0 and NIE = 0\n")
  cat("Alternative (H1):     NDE != 0 or NIE != 0\n\n")
  print(x$joint_tests, row.names = FALSE, digits = digits)
  cat("\nBootstrap and decomposition diagnostics\n")
  cat("----------------------------------------\n")
  print(x$diagnostics, row.names = FALSE, digits = digits)
  # Optionally print covariance matrices.
  if (isTRUE(x$cov_matrix)) {
    for (comp in names(x$bootstrap_covariance)) {
      cat("\nBootstrap covariance matrix: ", comp, "\n", sep = "")
      cat("----------------------------------------\n")
      print(x$bootstrap_covariance[[comp]], digits = digits)
      cat("\nNDE-NIE bootstrap covariance matrix: ",
          comp, "\n", sep = "")
      cat("----------------------------------------\n")
      print(x$covariance_NDE_NIE[[comp]], digits = digits)
      cat("\nAsymptotic covariance matrix Sigma_hat: ",
          comp, "\n", sep = "")
      cat("----------------------------------------\n")
      print(x$Sigma_hat[[comp]], digits = digits)
      cat("\nAsymptotic NDE-NIE covariance matrix: ",
          comp, "\n", sep = "")
      cat("----------------------------------------\n")
      print(x$V_gamma_hat[[comp]], digits = digits)
    }
  }
  cat("\nEffect definitions:\n")
  cat("  TE  = Total effect\n")
  cat("  NDE = Natural direct effect\n")
  cat("  NIE = Natural indirect effect\n")
  invisible(x)
}



################################################################################
# Plot method for Jmma objects
#
# type = "forest": marginal effect estimates and confidence intervals
# type = "joint":  joint NDE-NIE confidence ellipse(s)
################################################################################

plot.Jmma <- function(x,
                      type = c("forest", "joint"),
                      comparison = NULL,
                      conf_level = NULL,
                      digits = 3,
                      ...) {

  if (!inherits(x, "Jmma")) {
    stop("Object must inherit from class 'Jmma'.")
  }
  type <- match.arg(type)
  if (is.null(conf_level)) {
    conf_level <- x$settings$conf_level
  }
  if (length(conf_level) != 1L ||
      !is.numeric(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.")
  }
  # Forest plot: TE, NDE, and NIE
  if (type == "forest") {
    results <- x$results
    if (!is.null(comparison)) {
      if (any(!comparison %in% results$Comparison)) {
        stop("Unknown comparison supplied to the forest plot.")
      }
      results <- results[
        results$Comparison %in% comparison,
        ,
        drop = FALSE
      ]
    }
    if (nrow(results) == 0L) {
      stop("No results available to plot.")
    }
    labels <- paste(
      results$Comparison,
      results$Effect,
      sep = " : "
    )
    n_effects <- nrow(results)
    ypos <- rev(seq_len(n_effects))
    x_range <- range(
      results$Lower_CI,
      results$Upper_CI,
      0,
      finite = TRUE
    )
    if (diff(x_range) == 0) {
      x_range <- x_range + c(-1, 1)
    }
    padding <- 0.08 * diff(x_range)
    x_range <- x_range + c(-padding, padding)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 11, 4, 2) + 0.1)
    plot(
      results$Estimate,
      ypos,
      xlim = x_range,
      ylim = c(0.5, n_effects + 0.5),
      yaxt = "n",
      ylab = "",
      xlab = "Effect estimate",
      main = paste0(
        "Mediation effects (",
        100 * conf_level,
        "% confidence intervals)"
      ),
      pch = 19,
      ...
    )
    axis(
      side = 2,
      at = ypos,
      labels = labels,
      las = 1,
      cex.axis = 0.85
    )
    abline(v = 0, lty = 2, col = "gray50")
    segments(
      x0 = results$Lower_CI,
      y0 = ypos,
      x1 = results$Upper_CI,
      y1 = ypos,
      lwd = 1.5
    )
    points(
      results$Estimate,
      ypos,
      pch = 19
    )
    return(invisible(x))
  }
  # Joint confidence ellipse(s) for NDE and NIE
  if (type == "joint") {
    available <- names(x$inference)
    if (length(available) == 0L) {
      stop("No fitted comparisons are available.")
    }
    # If comparison is NULL, plot all comparisons.
    if (is.null(comparison)) {
      comparisons <- available
    } else {
      if (any(!comparison %in% available)) {
        stop(
          "Unknown comparison(s): ",
          paste(
            setdiff(comparison, available),
            collapse = ", "
          )
        )
      }
      # Preserve the order requested by the user.
      comparisons <- unique(comparison)
    }
    # Helper function to draw one joint confidence ellipse.
    draw_joint_ellipse <- function(comp) {
      inf <- x$inference[[comp]]
      center <- inf$gamma_hat[c("NDE", "NIE")]
      if (length(center) != 2L ||
          any(!is.finite(center))) {
        stop(
          "Invalid NDE/NIE estimates for comparison: ",
          comp
        )
      }
      center <- as.numeric(center)
      # Bootstrap covariance on the estimator scale.
      V <- inf$joint_confidence_region$
        covariance_estimator_scale
      if (!is.matrix(V) ||
          !all(dim(V) == c(2L, 2L)) ||
          any(!is.finite(V))) {
        stop(
          "Invalid NDE-NIE covariance matrix for comparison: ",
          comp
        )
      }
      V <- (V + t(V)) / 2
      eig <- eigen(V, symmetric = TRUE)
      eig_values <- eig$values
      if (any(!is.finite(eig_values)) ||
          min(eig_values) <= 0) {
        stop(
          "The NDE-NIE covariance matrix is not positive definite for: ",
          comp
        )
      }
      # Chi-square critical value for a two-dimensional region.
      chi2_crit <- qchisq(conf_level, df = 2)
      # Generate the joint confidence ellipse.
      angle <- seq(0, 2 * pi, length.out = 361)
      circle <- rbind(
        cos(angle),
        sin(angle)
      )
      transform <- eig$vectors %*%
        diag(sqrt(eig_values * chi2_crit), nrow = 2)
      ellipse <- sweep(
        transform %*% circle,
        MARGIN = 1,
        STATS = center,
        FUN = "+"
      )
      # Include zero in both axes for interpretation.
      x_range <- range(ellipse[1, ], 0)
      y_range <- range(ellipse[2, ], 0)
      x_padding <- max(diff(x_range) * 0.08, 0.01)
      y_padding <- max(diff(y_range) * 0.08, 0.01)
      x_range <- x_range + c(-x_padding, x_padding)
      y_range <- y_range + c(-y_padding, y_padding)
      plot(
        ellipse[1, ],
        ellipse[2, ],
        type = "n",
        xlim = x_range,
        ylim = y_range,
        xlab = "Natural direct effect (NDE)",
        ylab = "Natural indirect effect (NIE)",
        main = paste0(
          comp,
          " vs reference: ",
          100 * conf_level,
          "% joint confidence region"
        ),
        ...
      )
      abline(
        h = 0,
        v = 0,
        lty = 2,
        col = "gray50"
      )
      lines(
        ellipse[1, ],
        ellipse[2, ],
        lwd = 2
      )
      points(
        center[1],
        center[2],
        pch = 19
      )
      invisible(ellipse)
    }
    # Save and restore the user's graphics settings.
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    if (length(comparisons) == 1L) {
      # Single comparison: one large plot.
      par(mar = c(5, 5, 4, 2) + 0.1)
      draw_joint_ellipse(comparisons[1])
    } else {
      # Multiple comparisons: one panel per comparison.
      n_panels <- length(comparisons)
      n_rows <- ceiling(sqrt(n_panels))
      n_cols <- ceiling(n_panels / n_rows)
      par(
        mfrow = c(n_rows, n_cols),
        mar = c(4, 4, 3, 1) + 0.1
      )
      for (comp in comparisons) {
        draw_joint_ellipse(comp)
      }
    }
    return(invisible(x))
  }
}

################################################################################
################################################################################
