
################################################################################
# Run joint MMA mediation analysis with optional competing-risks RMST
# pseudo-values
################################################################################

jmma <- function(data,
                 exposure,
                 comparison = NULL,
                 reference,
                 covariates,
                 mediators,
                 outcome,
                 n_boot = 10,
                 conf_level = 0.95,
                 use_pseudo = TRUE,
                 tau = NULL,
                 event_code = 1,
                 time_var,
                 status_var,
                 ...) {

  # Validate inputs
  if (length(conf_level) != 1L ||
      !is.numeric(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop("conf_level must be between 0 and 1.")
  }
  if (length(n_boot) != 1L ||
      !is.numeric(n_boot) ||
      !is.finite(n_boot) ||
      n_boot < 2 ||
      n_boot != as.integer(n_boot)) {
    stop("n_boot must be an integer greater than or equal to 2.")
  }
  if (use_pseudo) {
    if (missing(time_var) || missing(status_var)) {
      stop(
        "Specify both time_var and status_var ",
        "when use_pseudo = TRUE."
      )
    }
    if (is.null(tau) ||
        !is.numeric(tau) ||
        length(tau) != 1L ||
        !is.finite(tau) ||
        tau <= 0) {
      stop(
        "Specify a positive, finite numeric tau ",
        "when use_pseudo = TRUE."
      )
    }
    if (length(event_code) != 1L ||
        is.na(event_code) ||
        event_code == 0) {
      stop("event_code must be one nonzero event-status code.")
    }
  }
  if (length(outcome) != 1L ||
      !is.character(outcome) ||
      is.na(outcome) ||
      !nzchar(outcome)) {
    stop("outcome must be a single nonempty column name.")
  }
  # Prepare complete-case analysis dataset
  if (use_pseudo) {
    dd <- prepare_mma_data(
      data = data,
      covariates = covariates,
      mediators = mediators,
      exposure = exposure,
      outcome = NULL,
      time_var = time_var,
      event_var = status_var
    )
    if (nrow(dd) < 2L) {
      stop(
        "Insufficient complete observations ",
        "to calculate RMST pseudo-values."
      )
    }
    # Calculate competing-risks RMST pseudo-values
    pseudo <- get_pseudo_rmst(
      data = dd,
      tau = tau,
      event_code = event_code,
      time_var = time_var,
      status_var = status_var
    )
    if (length(pseudo) != nrow(dd)) {
      stop(
        "The number of RMST pseudo-values does not ",
        "match the analysis sample size."
      )
    }
    # Assign pseudo-values as the MMA outcome
    dd[[outcome]] <- pseudo
    # Remove non-finite pseudo-values
    valid <- is.finite(dd[[outcome]])
    if (any(!valid)) {
      warning(
        sum(!valid),
        " observations removed because their pseudo-values ",
        "were not finite."
      )
    }
    dd <- dd[valid, , drop = FALSE]
    dd <- droplevels(dd)
  } else {
    dd <- prepare_mma_data(
      data = data,
      covariates = covariates,
      mediators = mediators,
      exposure = exposure,
      outcome = outcome
    )
    tau <- NULL
    event_code <- NULL
    time_var <- NULL
    status_var <- NULL
  }
  if (nrow(dd) < 2L) {
    stop("Insufficient complete observations for analysis.")
  }
  if (length(unique(dd[[exposure]])) < 2L) {
    stop("At least two exposure categories are required.")
  }
  if (!reference %in% as.character(dd[[exposure]])) {
    stop("Reference category not found in the analysis dataset.")
  }
  # Run MMA for each exposure comparison
  fits <- run_mma_multiple(
    data = dd,
    exposure = exposure,
    comparison = comparison,
    reference = reference,
    covariates = covariates,
    mediators = mediators,
    outcome = outcome,
    n_boot = as.integer(n_boot),
    ...
  )
  if (length(fits) == 0L) {
    stop("No MMA models were fitted.")
  }
  # Calculate marginal and joint inference for each comparison
  inference <- lapply(fits, function(fit) {
    joint_mediation_inference(
      fit = fit,
      conf_level = conf_level
    )
  })
  # Combine TE, NDE, and total NIE results
  results_main <- dplyr::bind_rows(
    lapply(names(inference), function(comp) {
      out <- inference[[comp]]$results
      out$Comparison <- comp
      out
    })
  ) %>%
    dplyr::select(
      Comparison,
      Effect,
      Estimate,
      SE,
      Lower_CI,
      Upper_CI
    )
  # Combine mediator-specific indirect-effect results
  mediator_results <- dplyr::bind_rows(
    lapply(names(inference), function(comp) {
      out <- inference[[comp]]$mediator_results
      if (is.null(out) || nrow(out) == 0L) {
        return(NULL)
      }
      out$Comparison <- comp
      out
    })
  )
  # Append mediator-specific effects to the main results
  if (nrow(mediator_results) > 0L) {
    results <- dplyr::bind_rows(
      results_main,
      mediator_results %>%
        dplyr::select(
          Comparison,
          Effect,
          Estimate,
          SE,
          Lower_CI,
          Upper_CI
        )
    )
  } else {
    results <- results_main
  }
  # Combine joint NDE-NIE Wald tests
  joint_tests <- dplyr::bind_rows(
    lapply(names(inference), function(comp) {
      out <- inference[[comp]]$joint_wald_test
      if (is.null(out) || nrow(out) == 0L) {
        return(NULL)
      }
      out$Comparison <- comp
      out
    })
  )
  if (nrow(joint_tests) > 0L) {
    joint_tests <- joint_tests %>%
      dplyr::select(
        Comparison,
        Test,
        Statistic,
        DF,
        P_value
      )
  }
  # Combine joint mediator-specific indirect-effect tests
  mediator_joint_tests <- dplyr::bind_rows(
    lapply(names(inference), function(comp) {
      out <- inference[[comp]]$mediator_joint_wald_test
      if (is.null(out) || nrow(out) == 0L) {
        return(NULL)
      }
      out$Comparison <- comp
      out
    })
  )
  if (nrow(mediator_joint_tests) > 0L) {
    mediator_joint_tests <- mediator_joint_tests %>%
      dplyr::select(
        Comparison,
        Test,
        Statistic,
        DF,
        P_value
      )
  }
  # Combine joint NDE and all mediator-specific effect tests
  nde_mediator_joint_tests <- dplyr::bind_rows(
    lapply(names(inference), function(comp) {
      out <- inference[[comp]]$nde_mediator_joint_wald_test
      if (is.null(out) || nrow(out) == 0L) {
        return(NULL)
      }
      out$Comparison <- comp
      out
    })
  )
  if (nrow(nde_mediator_joint_tests) > 0L) {
    nde_mediator_joint_tests <- nde_mediator_joint_tests %>%
      dplyr::select(
        Comparison,
        Test,
        Statistic,
        DF,
        P_value
      )
  }
  # Collect joint confidence regions by comparison
  joint_confidence_regions <- lapply(
    inference,
    function(x) x$joint_confidence_region
  )
  mediator_joint_confidence_regions <- lapply(
    inference,
    function(x) x$mediator_joint_confidence_region
  )
  nde_mediator_joint_confidence_regions <- lapply(
    inference,
    function(x) x$nde_mediator_joint_confidence_region
  )
  # Collect covariance matrices and diagnostics
  covariance_matrices <- lapply(inference, function(x) {
    list(
      overall_effects = x$bootstrap_covariance,
      NDE_NIE = x$covariance_NDE_NIE,
      mediator_effects = x$covariance_mediator_effects,
      NDE_mediator_effects =
        x$covariance_NDE_mediator_effects
    )
  })
  decomposition_diagnostics <- lapply(
    inference,
    function(x) x$decomposition_diagnostic
  )
  # Store settings and return the Jmma object
  out <- structure(
    list(
      # Analysis data and fitted models
      data = dd,
      fits = fits,
      inference = inference,
      # Marginal mediation estimates
      results = results,
      mediator_results = mediator_results,
      # Joint Wald test tables
      joint_tests = joint_tests,
      mediator_joint_tests = mediator_joint_tests,
      nde_mediator_joint_tests = nde_mediator_joint_tests,
      # Joint confidence regions
      joint_confidence_regions = joint_confidence_regions,
      mediator_joint_confidence_regions =
        mediator_joint_confidence_regions,
      nde_mediator_joint_confidence_regions =
        nde_mediator_joint_confidence_regions,
      # Covariance matrices and diagnostics
      covariance_matrices = covariance_matrices,
      decomposition_diagnostics = decomposition_diagnostics,
      # Analysis settings
      settings = list(
        exposure = exposure,
        reference = reference,
        comparisons = names(fits),
        covariates = covariates,
        mediators = mediators,
        outcome = outcome,
        use_pseudo = use_pseudo,
        tau = tau,
        event_code = event_code,
        time_var = time_var,
        status_var = status_var,
        n_boot = as.integer(n_boot),
        conf_level = conf_level,
        analysis_sample_size = nrow(dd)
      )
    ),
    class = "Jmma"
  )
  return(out)
}

################################################################################
################################################################################
