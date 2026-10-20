##
## run function
##

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

  # Validate pseudo-value settings
  if (use_pseudo && is.null(tau)) {
    stop("Specify tau when use_pseudo = TRUE.")
  }
  if (use_pseudo && !is.numeric(tau)) {
    stop("tau must be numeric.")
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
    # Calculate competing-risks RMST pseudo-values
    pseudo <- get_pseudo_rmst(
      data = dd,
      tau = tau,
      event_code = event_code,
      time_var = time_var,
      status_var = status_var
    )
    # Add pseudo-values as the MMA outcome
    dd[[outcome]] <- pseudo
    # Remove invalid pseudo-values, if any
    valid <- is.finite(dd[[outcome]])
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
  }
  if (nrow(dd) < 2L) {
    stop("Insufficient complete observations for analysis.")
  }
  if (length(unique(dd[[exposure]])) < 2L) {
    stop("At least two exposure categories are required.")
  }
  # Run MMA for each comparison against the reference
  fits <- run_mma_multiple(
    data = dd,
    exposure = exposure,
    comparison = comparison,
    reference = reference,
    covariates = covariates,
    mediators = mediators,
    outcome = outcome,
    n_boot = n_boot,
    ...
  )
  # Joint and marginal mediation inference for each comparison
  inference <- lapply(fits, function(fit) {
    joint_mediation_inference(
      fit = fit,
      conf_level = conf_level
    )
  })
  # Combine marginal estimates and confidence intervals
  results <- dplyr::bind_rows(
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
  # Combine joint NDE-NIE Wald tests
  joint_tests <- dplyr::bind_rows(
    lapply(names(inference), function(comp) {
      out <- inference[[comp]]$joint_wald_test
      out$Comparison <- comp
      out
    })
  ) %>%
    dplyr::select(
      Comparison,
      Test,
      Statistic,
      DF,
      P_value
    )
  # Return fitted models, inference, and analysis dataset
  out <- list(
    data = dd,
    fits = fits,
    inference = inference,
    results = results,
    joint_tests = joint_tests,
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
      n_boot = n_boot,
      conf_level = conf_level,
      analysis_sample_size = nrow(dd)
    )
  )
  class(out) <- "Jmma"
  out
}

