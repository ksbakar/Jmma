
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
    #cat("Event code:          ", settings$event_code, "\n", sep = "")
  } else {
    cat("Outcome type:        Observed outcome\n")
  }
  #cat("Bootstrap replicates requested: ", settings$n_boot, "\n", sep = "")
  cat("Confidence level:     ", 100 * settings$conf_level, "%\n", sep = "")
  cat("Comparisons:          ", length(settings$comparisons), "\n", sep = "")
  invisible(x)
}



################################################################################
# Summary method for Jmma objects
################################################################################

summary.Jmma <- function(object,
                         cov_matrix = FALSE,
                         diagnostics = FALSE,
                         wald_test = FALSE,
                         ...) {

  if (!inherits(object, "Jmma")) {
    stop("Object must inherit from class 'Jmma'.")
  }
  # Validate options
  if (!is.logical(cov_matrix) ||
      length(cov_matrix) != 1L ||
      is.na(cov_matrix)) {
    stop("'cov_matrix' must be TRUE or FALSE.")
  }
  if (!is.logical(diagnostics) ||
      length(diagnostics) != 1L ||
      is.na(diagnostics)) {
    stop("'diagnostics' must be TRUE or FALSE.")
  }
  if (!is.logical(wald_test) ||
      length(wald_test) != 1L ||
      is.na(wald_test)) {
    stop("'wald_test' must be TRUE or FALSE.")
  }
  # Collect diagnostics for each comparison
  diagnostics_list <- lapply(names(object$inference), function(comp) {
    inf <- object$inference[[comp]]
    dec <- inf$decomposition_diagnostic
    data.frame(
      Comparison = comp,
      Bootstrap_replicates = inf$n_boot,
      Sample_size = inf$n,
      Covariance_rank = inf$estimated_covariance_rank,
      TE_minus_NDE_minus_NIE =
        if (!is.null(dec)) dec$point_estimate_error else NA_real_,
      Max_abs_bootstrap_decomposition_error =
        if (!is.null(dec)) {
          dec$max_absolute_bootstrap_error
        } else {
          NA_real_
        },
      row.names = NULL
    )
  })
  diagnostics_table <- if (length(diagnostics_list)) {
    do.call(rbind, diagnostics_list)
  } else {
    data.frame()
  }
  rownames(diagnostics_table) <- NULL
  # Mediator-specific indirect effects
  mediator_results <- if (!is.null(object$mediator_results)) {
    object$mediator_results
  } else if (
    !is.null(object$results) &&
    "Effect" %in% names(object$results)
  ) {
    object$results[
      grepl("^Indirect: ", object$results$Effect),
      ,
      drop = FALSE
    ]
  } else {
    data.frame(
      Comparison = character(),
      Effect = character(),
      Estimate = numeric(),
      SE = numeric(),
      Lower_CI = numeric(),
      Upper_CI = numeric()
    )
  }
  # Helper: safely retrieve a field from the Jmma object
  get_field <- function(name, default = NULL) {
    value <- object[[name]]
    if (is.null(value)) default else value
  }
  # Existing covariance matrices
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
  # All three joint Wald tests
  joint_tests <- get_field("joint_tests")
  mediator_joint_tests <- get_field("mediator_joint_tests")
  nde_mediator_joint_tests <- get_field("nde_mediator_joint_tests")
  # Joint confidence regions and additional covariance outputs
  joint_confidence_regions <-
    get_field("joint_confidence_regions")
  mediator_joint_confidence_regions <-
    get_field("mediator_joint_confidence_regions")
  nde_mediator_joint_confidence_regions <-
    get_field("nde_mediator_joint_confidence_regions")
  # Return summary object
  out <- list(
    settings = object$settings,
    results = object$results,
    mediator_results = mediator_results,
    # Joint tests
    joint_tests = joint_tests,
    mediator_joint_tests = mediator_joint_tests,
    nde_mediator_joint_tests = nde_mediator_joint_tests,
    # Joint confidence regions
    joint_confidence_regions = joint_confidence_regions,
    mediator_joint_confidence_regions =
      mediator_joint_confidence_regions,
    nde_mediator_joint_confidence_regions =
      nde_mediator_joint_confidence_regions,
    # Diagnostics and covariance matrices
    diagnostics = diagnostics_table,
    bootstrap_covariance = bootstrap_covariance,
    covariance_NDE_NIE = covariance_NDE_NIE,
    Sigma_hat = Sigma_hat,
    V_gamma_hat = V_gamma_hat,
    # Full inference results
    inference = object$inference,
    cov_matrix = cov_matrix,
    show_diagnostics = diagnostics,
    show_wald_test = wald_test
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
  # Model settings
  cat("\nModel settings\n")
  cat("----------------------------------------\n")
  cat("Exposure:             ", x$settings$exposure, "\n", sep = "")
  cat("Reference:            ", x$settings$reference, "\n", sep = "")
  cat("Outcome:              ", x$settings$outcome, "\n", sep = "")
  cat("Analysis sample size:  ",
      x$settings$analysis_sample_size, "\n", sep = "")
  cat("Confidence level:      ",
      100 * x$settings$conf_level, "%\n", sep = "")
  # Overall effects
  cat("\nTotal, direct, and total indirect effects\n")
  cat("----------------------------------------\n")
  main_results <- x$results[
    x$results$Effect %in% c("TE", "NDE", "NIE"),
    ,
    drop = FALSE
  ]
  print(main_results, row.names = FALSE, digits = digits)
  # Mediator-specific indirect effects
  cat("\nMediator-specific indirect effects\n")
  cat("----------------------------------------\n")
  if (!is.null(x$mediator_results) &&
      nrow(x$mediator_results) > 0L) {
    print(
      x$mediator_results,
      row.names = FALSE,
      digits = digits
    )
  } else {
    cat("No mediator-specific indirect effects available.\n")
  }
  # Helper for safely printing joint tests
  print_joint_test <- function(test, title, null_text, alt_text) {
    cat("\n", title, "\n", sep = "")
    cat("----------------------------------------\n")
    cat("Null hypothesis (H0): ", null_text, "\n", sep = "")
    cat("Alternative (H1):     ", alt_text, "\n\n", sep = "")
    if (!is.null(test) &&
        is.data.frame(test) &&
        nrow(test) > 0L) {
      print(test, row.names = FALSE, digits = digits)
    } else {
      cat("Joint test results are not available.\n")
    }
  }
  # Optional wald_test
  if (isTRUE(x$show_wald_test)) {
    # Joint test 1: NDE and total NIE
    print_joint_test(
      test = x$joint_tests,
      title = "Overall joint Wald test",
      null_text = "NDE = 0 and total NIE = 0",
      alt_text = "NDE != 0 or total NIE != 0"
    )
    # Joint test 2: all mediator-specific effects
    print_joint_test(
      test = x$mediator_joint_tests,
      title = "Mediator-specific joint Wald test",
      null_text = "All mediator-specific indirect effects = 0",
      alt_text = "At least one mediator-specific indirect effect != 0"
    )
    # Joint test 3: NDE and all mediator-specific effects
    print_joint_test(
      test = x$nde_mediator_joint_tests,
      title = "Combined direct and mediator-specific joint Wald test",
      null_text = paste(
        "NDE = 0 and all mediator-specific",
        "indirect effects = 0"
      ),
      alt_text = paste(
        "NDE != 0 or at least one mediator-specific",
        "indirect effect != 0"
      )
    )
  }
  # Optional diagnostics
  if (isTRUE(x$show_diagnostics)) {
    cat("\nBootstrap and decomposition diagnostics\n")
    cat("----------------------------------------\n")
    if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
      print(
        x$diagnostics,
        row.names = FALSE,
        digits = digits
      )
    } else {
      cat("Diagnostics are not available.\n")
    }
  }
  # Optional covariance matrices
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
      # Print additional covariance matrices when present
      inf <- x$inference[[comp]]
      additional_covariances <- inf[
        grepl("covariance|Sigma|V_gamma", names(inf),
              ignore.case = TRUE)
      ]
      if (length(additional_covariances) > 0L) {
        for (nm in names(additional_covariances)) {
          # Avoid printing matrices already displayed above
          if (nm %in% c(
            "bootstrap_covariance",
            "covariance_NDE_NIE",
            "Sigma_hat",
            "V_gamma_hat"
          )) {
            next
          }
          cat("\n", nm, ": ", comp, "\n", sep = "")
          cat("----------------------------------------\n")
          print(additional_covariances[[nm]], digits = digits)
        }
      }
    }
  }
  # Effect definitions
  cat("\nEffect definitions:\n")
  cat("  TE  = Total effect\n")
  cat("  NDE = Natural direct effect\n")
  cat("  NIE = Total natural indirect effect\n")
  invisible(x)
}


################################################################################
# Plot method for Jmma objects using ggplot2
#
# type = "forest":   TE, NDE, and total NIE
# type = "indirect": mediator-specific indirect effects
# type = "joint":    joint NDE-NIE confidence ellipse(s)
################################################################################

plot.Jmma <- function(x,
                      type = c("forest", "indirect", "joint"),
                      comparison = NULL,
                      conf_level = NULL,
                      digits = 3,
                      ...) {

  if (!inherits(x, "Jmma")) {
    stop("Object must inherit from class 'Jmma'.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
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
  # Common ggplot2 theme
  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )
  # Helper for forest plots
  forest_plot <- function(results, xlab, title) {
    labels <- paste(
      results$Comparison,
      results$Effect,
      sep = " : "
    )
    results$Label <- factor(
      labels,
      levels = rev(unique(labels))
    )
    ggplot2::ggplot(
      results,
      ggplot2::aes(y = Label)
    ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray50",
        linewidth = 0.6
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = Lower_CI,
          xend = Upper_CI,
          yend = Label
        ),
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = Estimate),
        size = 3
      ) +
      ggplot2::labs(
        title = title,
        x = xlab,
        y = NULL
      ) +
      base_theme
  }
  # Overall effects forest plot
  if (type == "forest") {
    results <- x$results[
      x$results$Effect %in% c("TE", "NDE", "NIE"),
      ,
      drop = FALSE
    ]
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
      stop(
        "No total, direct, or total indirect effects available to plot."
      )
    }
    return(
      forest_plot(
        results = results,
        xlab = "Effect estimate",
        title = paste0(
          "Overall mediation effects (",
          100 * conf_level,
          "% CI)"
        )
      )
    )
  }
  # Mediator-specific indirect effects forest plot
  if (type == "indirect") {
    results <- x$mediator_results
    # Backward compatibility
    if (is.null(results)) {
      if ("Effect" %in% names(x$results)) {
        results <- x$results[
          grepl("^Indirect: ", x$results$Effect),
          ,
          drop = FALSE
        ]
      } else {
        stop(
          "Mediator-specific indirect effects are not available."
        )
      }
    }
    if (nrow(results) == 0L) {
      stop(
        "No mediator-specific indirect effects available to plot."
      )
    }
    required_cols <- c(
      "Comparison",
      "Effect",
      "Estimate",
      "Lower_CI",
      "Upper_CI"
    )
    missing_cols <- setdiff(required_cols, names(results))
    if (length(missing_cols) > 0L) {
      stop(
        "Mediator results are missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    # Filter comparisons
    if (!is.null(comparison)) {
      available <- unique(as.character(results$Comparison))
      if (any(!comparison %in% available)) {
        stop(
          "Unknown comparison(s): ",
          paste(
            setdiff(comparison, available),
            collapse = ", "
          )
        )
      }
      results <- results[
        results$Comparison %in% comparison,
        ,
        drop = FALSE
      ]
    }
    if (nrow(results) == 0L) {
      stop(
        "No mediator-specific indirect effects for selected comparisons."
      )
    }
    # Exclude the joint mediator effect, if present
    results <- results[
      !grepl("^Indirect: j1$", results$Effect) &
        !grepl("^j1$", results$Effect),
      ,
      drop = FALSE
    ]
    if (nrow(results) == 0L) {
      stop(
        "No individual mediator-specific effects available to plot."
      )
    }
    # Use the Mediator column when available
    mediator_label <- if ("Mediator" %in% names(results)) {
      as.character(results$Mediator)
    } else {
      sub(
        "^Indirect: ",
        "",
        as.character(results$Effect)
      )
    }
    results$Effect <- mediator_label
    return(
      forest_plot(
        results = results,
        xlab = "Mediator-specific indirect effect",
        title = paste0(
          "Mediator-specific indirect effects (",
          100 * conf_level,
          "% CI)"
        )
      )
    )
  }
  # Joint confidence ellipse(s) for NDE and NIE
  if (type == "joint") {
    available <- names(x$inference)
    if (length(available) == 0L) {
      stop("No fitted comparisons are available.")
    }
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
      comparisons <- unique(comparison)
    }
    # Critical value for a two-dimensional confidence region
    chi2_crit <- stats::qchisq(conf_level, df = 2)
    angle <- seq(0, 2 * pi, length.out = 361)
    circle <- rbind(
      cos(angle),
      sin(angle)
    )
    ellipse_data <- vector("list", length(comparisons))
    center_data <- vector("list", length(comparisons))
    for (i in seq_along(comparisons)) {
      comp <- comparisons[i]
      inf <- x$inference[[comp]]
      # Center of the joint confidence ellipse
      center <- inf$gamma_hat[c("NDE", "NIE")]
      if (length(center) != 2L ||
          any(!is.finite(center))) {
        stop(
          "Invalid NDE/NIE estimates for comparison: ",
          comp
        )
      }
      center <- as.numeric(center)
      # Covariance matrix on the estimator scale
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
      if (any(!is.finite(eig$values)) ||
          min(eig$values) <= 0) {
        stop(
          "The NDE-NIE covariance matrix is not positive definite for: ",
          comp
        )
      }
      # Construct the ellipse
      transform <- eig$vectors %*%
        diag(
          sqrt(eig$values * chi2_crit),
          nrow = 2
        )
      ellipse <- sweep(
        transform %*% circle,
        MARGIN = 1,
        STATS = center,
        FUN = "+"
      )
      ellipse_data[[i]] <- data.frame(
        NDE = ellipse[1, ],
        NIE = ellipse[2, ],
        Comparison = comp
      )
      center_data[[i]] <- data.frame(
        NDE = center[1],
        NIE = center[2],
        Comparison = comp
      )
    }
    ellipse_df <- do.call(rbind, ellipse_data)
    center_df <- do.call(rbind, center_data)
    ellipse_df$Comparison <- factor(
      ellipse_df$Comparison,
      levels = comparisons
    )
    center_df$Comparison <- factor(
      center_df$Comparison,
      levels = comparisons
    )
    p <- ggplot2::ggplot(
      ellipse_df,
      ggplot2::aes(
        x = NDE,
        y = NIE,
        group = Comparison
      )
    ) +
      ggplot2::geom_polygon(
        fill = "steelblue",
        alpha = 0.20,
        color = NA
      ) +
      ggplot2::geom_path(
        color = "steelblue",
        linewidth = 1
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "gray50",
        linewidth = 0.6
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray50",
        linewidth = 0.6
      ) +
      ggplot2::geom_point(
        data = center_df,
        ggplot2::aes(x = NDE, y = NIE),
        inherit.aes = FALSE,
        size = 2.5
      ) +
      ggplot2::labs(
        title = paste0(
          "Joint NDE-NIE confidence region (",
          100 * conf_level,
          "% CI)"
        ),
        x = "Natural direct effect (NDE)",
        y = "Natural indirect effect (NIE)"
      ) +
      base_theme
    # Separate panels when multiple comparisons are selected
    if (length(comparisons) > 1L) {
      p <- p +
        ggplot2::facet_wrap(
          ~Comparison,
          scales = "free"
        )
    } else {
      p <- p +
        ggplot2::labs(
          subtitle = paste0(
            comparisons[1],
            " vs reference"
          )
        )
    }
    return(p)
  }
}

################################################################################
################################################################################
