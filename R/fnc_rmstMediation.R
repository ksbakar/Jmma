

################################################################################
# Calculate restricted mean survival time from competing risks
# data: data.frame containing time and status variables
# tau: restriction time
# time_var: name of follow-up time variable
# status_var: name of competing risk status variable
# event_code: event of interest
################################################################################
get_rmst_cr <- function(data,
                        tau,
                        time_var = "GF_time_cr",
                        status_var = "GF_status_cr",
                        event_code = 1) {

  cif_obj <- cmprsk::cuminc(
    ftime   = data[[time_var]],
    fstatus = data[[status_var]],
    cencode = 0
  )
  curve_name <- names(cif_obj)[
    grepl(paste0(" ", event_code, "$"), names(cif_obj))
  ]
  if(length(curve_name) == 0){
    stop("Event code not found in cumulative incidence object")
  }
  cif_event <- cif_obj[[curve_name]]
  time <- c(0, cif_event$time)
  cif  <- c(0, cif_event$est)
  keep <- time <= tau
  time <- time[keep]
  cif  <- cif[keep]
  if(max(time) < tau){

    time <- c(time, tau)
    cif <- c(cif, tail(cif,1))

  }
  surv <- 1 - cif
  rmst <- sum(
    diff(time) * surv[-length(surv)]
  )
  return(rmst)
}
################################################################################
# Calculate pseudo values for RMST
# data: dataset
# tau: restriction time
# event_code: event of interest
################################################################################
get_pseudo_rmst <- function(data,
                            tau,
                            event_code = 1,
                            time_var = "GF_time_cr",
                            status_var = "GF_status_cr"){
  n <- nrow(data)
  rmst_full <- get_rmst_cr(
    data,
    tau,
    time_var,
    status_var,
    event_code
  )
  rmst_loo <- sapply(seq_len(n), function(i){
    get_rmst_cr(
      data[-i, ],
      tau,
      time_var,
      status_var,
      event_code
    )
  })
  pseudo <- n * rmst_full -
    (n-1) * rmst_loo
  return(pseudo)
}
################################################################################
# Prepare data for MMA mediation analysis
################################################################################

prepare_mma_data <- function(data,
                             covariates,
                             mediators,
                             exposure,
                             outcome = NULL,
                             time_var = NULL,
                             event_var = NULL) {
  required_vars <- unique(c(
    exposure,
    covariates,
    mediators,
    outcome,
    time_var,
    event_var
  ))
  required_vars <- required_vars[
    !is.na(required_vars) & nzchar(required_vars)
  ]
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "Variables not found in data: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  data %>%
    dplyr::select(dplyr::all_of(required_vars)) %>%
    tidyr::drop_na() %>%
    droplevels()
}
################################################################################
# Run MMA for binary exposure comparison
################################################################################
run_mma_single <- function(data,
                           exposure,
                           comparison,
                           reference,
                           covariates,
                           mediators,
                           outcome,
                           n_boot,
                           ...){
  dd <- data %>%
    dplyr::mutate(
      pred_bin = as.integer(.data[[exposure]] %in% comparison)
    )
  x <- dd %>%
    dplyr::select(
      all_of(c(covariates, mediators))
    )
  # Remove variables with only one level
  keep_vars <- sapply(x, function(z){
    if(is.factor(z)){
      nlevels(droplevels(z)) > 1
    } else {
      length(unique(z)) > 1
    }
  })
  removed <- names(x)[!keep_vars]
  if(length(removed)>0){
    warning(
      "Removed variables with only one level: ",
      paste(removed, collapse=", ")
    )
    x <- x[,keep_vars,drop=FALSE]
    mediators <- intersect(
      mediators,
      names(x)
    )
  }
  y <- dd[[outcome]]
  pred <- dd$pred_bin
  set.seed(123)
  mma::mma(
    x=x,
    y=y,
    pred=pred,
    mediator=mediators,
    binmed=mediators,
    jointm=list(n=1,j1=mediators),
    binref=rep(0,length(mediators)),
    predref="0",
    alpha = 0.99,
    alpha2 = 0.99,
    n2=n_boot,
    ...
  )
}
################################################################################
# Run MMA for multiple exposure categories
################################################################################

run_mma_multiple <- function(data,
                             exposure,
                             comparison = NULL,
                             reference,
                             covariates,
                             mediators,
                             outcome,
                             n_boot,
                             ...) {

  exposure_values <- as.character(data[[exposure]])
  if (!reference %in% exposure_values) {
    stop("Reference category not found in exposure.")
  }
  if (is.null(comparison)) {
    comparisons <- setdiff(unique(exposure_values), reference)
  } else {
    comparisons <- setdiff(
      unique(as.character(comparison)),
      reference
    )
  }
  if (length(comparisons) == 0L) {
    stop("No comparison categories available.")
  }
  if (any(!comparisons %in% exposure_values)) {
    stop("One or more comparison categories are absent from the data.")
  }
  fits <- lapply(seq_along(comparisons), function(i) {
    comp <- comparisons[i]
    cat("\n")
    cat("==================================================\n")
    cat("JMMA comparison ", i, " of ", length(comparisons), "\n",
        sep = "")
    cat("Comparison category: ", comp, "\n", sep = "")
    cat("Reference category:  ", reference, "\n", sep = "")
    keep <- exposure_values %in% c(reference, comp)
    dd <- data[keep, , drop = FALSE]
    dd[[exposure]] <- droplevels(
      factor(as.character(dd[[exposure]]))
    )
    cat("Sample size:         ", nrow(dd), "\n", sep = "")
    cat("Outcome:             ", outcome, "\n", sep = "")
    cat("Mediators:           ",
        paste(mediators, collapse = ", "), "\n", sep = "")
    cat("Starting JMMA model...\n")
    cat("==================================================\n")
    fit <- run_mma_single(
      data = dd,
      exposure = exposure,
      comparison = comp,
      reference = reference,
      covariates = covariates,
      mediators = mediators,
      outcome = outcome,
      n_boot = n_boot,
      ...
    )
    cat("Finished JMMA model for: ", comp, "\n", sep = "")
    fit
  })
  names(fits) <- comparisons
  fits
}
################################################################################
# R function for joint adjustment and standard errors
################################################################################

# Joint multivariate inference for mediation effects
# TE, NDE, and NIE
#
# - Paired bootstrap covariance
# - Marginal standard errors and Wald confidence intervals
# - Joint Wald inference for NDE and NIE
# - No covariance-weighted adjustment

joint_mediation_inference <- function(fit, conf_level = 0.95) {

  # Extract original estimates and bootstrap results
  bootsresults <- fit$a.binx$bootsresults
  nie_col <- "y1.all"
  theta_hat <- c(
    TE  = as.numeric(fit$a.binx$estimation$te),
    NDE = as.numeric(fit$a.binx$estimation$de),
    NIE = as.numeric(fit$a.binx$estimation$ie[1])
  )
  if (length(theta_hat) != 3L || any(!is.finite(theta_hat))) {
    stop("TE, NDE, and NIE point estimates must all be finite.")
  }
  # Extract bootstrap estimates.
  te_boot  <- as.numeric(bootsresults$te)
  nde_boot <- as.numeric(bootsresults$de)
  nie_boot <- as.numeric(bootsresults$ie$pred[, nie_col])
  if (length(te_boot) != length(nde_boot) ||
      length(te_boot) != length(nie_boot)) {
    stop("Bootstrap estimates must have the same number of replicates.")
  }
  boot_mat <- cbind(
    TE  = te_boot,
    NDE = nde_boot,
    NIE = nie_boot
  )
  # Retain paired, complete bootstrap replicates.
  boot_mat <- boot_mat[
    complete.cases(boot_mat) & apply(boot_mat, 1, function(x) all(is.finite(x))),
    , drop = FALSE]
  B <- nrow(boot_mat)
  if (B < 2L) {
    stop("At least two complete bootstrap replicates are required.")
  }
  # Sample size and bootstrap covariance
  n <- nrow(fit$a.binx$data$x)
  if (length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("Could not determine a valid original sample size.")
  }
  # Covariance of the bootstrap estimators, on the estimator scale.
  cov_boot <- cov(boot_mat)
  # Estimate of the asymptotic covariance of sqrt(n) * theta_hat.
  Sigma_hat <- n * cov_boot
  # Numerical symmetrization.
  cov_boot <- (cov_boot + t(cov_boot)) / 2
  Sigma_hat <- (Sigma_hat + t(Sigma_hat)) / 2
  effect_names <- c("TE", "NDE", "NIE")
  dimnames(cov_boot) <- list(effect_names, effect_names)
  dimnames(Sigma_hat) <- list(effect_names, effect_names)
  # Marginal standard errors and Wald confidence intervals
  SE <- sqrt(pmax(diag(cov_boot), 0))
  names(SE) <- effect_names
  alpha <- 1 - conf_level
  z_crit <- qnorm(1 - alpha / 2)
  lower <- theta_hat - z_crit * SE
  upper <- theta_hat + z_crit * SE
  results <- data.frame(
    Effect = effect_names,
    Estimate = unname(theta_hat),
    SE = unname(SE),
    Lower_CI = unname(lower),
    Upper_CI = unname(upper),
    row.names = NULL
  )
  # Joint multivariate inference for NDE and NIE
  # NDE and NIE are the two independent effects.
  gamma_hat <- theta_hat[c("NDE", "NIE")]
  gamma_boot <- boot_mat[, c("NDE", "NIE"), drop = FALSE]
  # Covariance of the two effect estimators.
  cov_gamma_boot <- cov(gamma_boot)
  # Asymptotic covariance of sqrt(n) * (NDE, NIE).
  V_gamma_hat <- n * cov_gamma_boot
  cov_gamma_boot <- (cov_gamma_boot + t(cov_gamma_boot)) / 2
  V_gamma_hat <- (V_gamma_hat + t(V_gamma_hat)) / 2
  dimnames(cov_gamma_boot) <- list(
    c("NDE", "NIE"), c("NDE", "NIE")
  )
  dimnames(V_gamma_hat) <- list(
    c("NDE", "NIE"), c("NDE", "NIE")
  )
  # Check whether the covariance is positive definite.
  eig_gamma <- eigen(
    cov_gamma_boot,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (any(!is.finite(eig_gamma)) || min(eig_gamma) <= 0) {
    stop(
      paste(
        "The NDE-NIE bootstrap covariance matrix is not positive definite.",
        "Check the bootstrap estimates and numerical conditioning."
      )
    )
  }
  # Joint Wald test:
  # H0: NDE = 0 and NIE = 0.
  # The test has 2 degrees of freedom.
  joint_wald_stat <- drop(
    t(gamma_hat) %*%
      solve(cov_gamma_boot, gamma_hat)
  )
  joint_p_value <- pchisq(
    joint_wald_stat,
    df = 2,
    lower.tail = FALSE
  )
  joint_test <- data.frame(
    Test = "Joint Wald test: NDE = 0 and NIE = 0",
    Statistic = joint_wald_stat,
    DF = 2,
    P_value = joint_p_value
  )
  # Joint Wald confidence region for NDE and NIE
  # Region:
  # (gamma_hat - gamma)' cov_gamma_boot^{-1}
  # (gamma_hat - gamma) <= chi-square critical value.
  #
  # This defines an ellipse in the NDE-NIE plane.
  chi2_crit <- qchisq(conf_level, df = 2)
  joint_region <- list(
    confidence_level = conf_level,
    chi_square_critical_value = chi2_crit,
    center = gamma_hat,
    covariance_estimator_scale = cov_gamma_boot,
    covariance_asymptotic_scale = V_gamma_hat,
    equation = paste0(
      "(gamma_hat - gamma)' inverse(cov_gamma_boot) ",
      "(gamma_hat - gamma) <= ", signif(chi2_crit, 6)
    )
  )
  # Check the mediation decomposition
  # Point-estimate discrepancy.
  point_decomposition_error <-
    unname(theta_hat["TE"] - theta_hat["NDE"] - theta_hat["NIE"])
  # Bootstrap decomposition discrepancies.
  bootstrap_decomposition_error <-
    boot_mat[, "TE"] -  boot_mat[, "NDE"] -  boot_mat[, "NIE"]
  decomposition_diagnostic <- list(
    point_estimate_error = point_decomposition_error,
    max_absolute_bootstrap_error =
      max(abs(bootstrap_decomposition_error)),
    bootstrap_errors = bootstrap_decomposition_error
  )
  # Rank diagnostic for the full effect covariance matrix.
  eig_full <- eigen(
    cov_boot,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  tolerance <- max(abs(eig_full)) * 1e-8
  estimated_rank <- sum(eig_full > tolerance)
  # Return results
  list(
    results = results,
    theta_hat = theta_hat,
    # Covariance of the estimators and asymptotic covariance.
    bootstrap_covariance = cov_boot,
    Sigma_hat = Sigma_hat,
    # Joint inference for NDE and NIE.
    gamma_hat = gamma_hat,
    covariance_NDE_NIE = cov_gamma_boot,
    V_gamma_hat = V_gamma_hat,
    joint_wald_test = joint_test,
    joint_confidence_region = joint_region,
    # Diagnostics.
    decomposition_diagnostic = decomposition_diagnostic,
    estimated_covariance_rank = estimated_rank,
    bootstrap_matrix = boot_mat,
    n = n,
    n_boot = B,
    conf_level = conf_level
  )
}
################################################################################
################################################################################

