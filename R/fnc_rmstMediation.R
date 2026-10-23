

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

  # Validate inputs and extract MMA estimates and bootstraps
  if (length(conf_level) != 1L ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop("conf_level must be a single number between 0 and 1.")
  }
  bootsresults <- fit$a.binx$bootsresults
  estimation <- fit$a.binx$estimation
  if (is.null(bootsresults) || is.null(estimation)) {
    stop("Missing MMA estimation or bootstrap results.")
  }
  # Extract the 'est' row from an MMA effect table
  extract_mma_estimates <- function(x, effect_name = "effect") {
    if (is.null(x)) {
      stop("MMA is missing the ", effect_name, " estimates.")
    }
    if (is.data.frame(x) || is.matrix(x)) {
      rn <- rownames(x)
      if (!is.null(rn) && "est" %in% rn) {
        ans <- as.numeric(x["est", , drop = TRUE])
        names(ans) <- colnames(x)
        return(ans)
      }
      # Also allow a one-row matrix/data frame
      if (nrow(x) == 1L) {
        ans <- as.numeric(x[1, , drop = TRUE])
        names(ans) <- colnames(x)
        return(ans)
      }
      stop("Could not identify the 'est' row in ", effect_name, ".")
    }
    ans <- as.numeric(x)
    names(ans) <- names(x)
    ans
  }
  # Point estimates for TE, NDE, and total NIE
  te_hat <- as.numeric(estimation$te)[1]
  nde_hat <- as.numeric(estimation$de)[1]
  ie_hat <- extract_mma_estimates(
    estimation$ie,
    "indirect effects"
  )
  if (is.null(names(ie_hat)) ||
      any(is.na(names(ie_hat)) | names(ie_hat) == "")) {
    stop("MMA indirect-effect estimates must have named columns.")
  }
  if (!"y1.all" %in% names(ie_hat)) {
    stop("Could not find 'y1.all' in MMA indirect-effect estimates.")
  }
  theta_hat <- c(
    TE = te_hat,
    NDE = nde_hat,
    NIE = unname(ie_hat["y1.all"])
  )
  if (any(!is.finite(theta_hat))) {
    stop("TE, NDE, and NIE point estimates must all be finite.")
  }
  # Extract paired bootstrap estimates
  te_boot <- as.numeric(bootsresults$te)
  nde_boot <- as.numeric(bootsresults$de)
  ie_boot <- bootsresults$ie$pred
  if (is.null(dim(ie_boot))) {
    stop("MMA bootstrap indirect effects must be a matrix.")
  }
  ie_boot <- as.matrix(ie_boot)
  if (!"y1.all" %in% colnames(ie_boot)) {
    stop("Could not find 'y1.all' in bootstrap indirect effects.")
  }
  if (length(te_boot) != length(nde_boot) ||
      length(te_boot) != nrow(ie_boot)) {
    stop("Bootstrap estimates must have the same number of replicates.")
  }
  # Keep paired replicates for TE, NDE, and total NIE
  boot_mat <- cbind(
    TE = te_boot,
    NDE = nde_boot,
    NIE = as.numeric(ie_boot[, "y1.all"])
  )
  valid <- complete.cases(boot_mat) &
    apply(boot_mat, 1, function(z) all(is.finite(z)))
  boot_mat <- boot_mat[valid, , drop = FALSE]
  B <- nrow(boot_mat)
  if (B < 2L) {
    stop("At least two complete bootstrap replicates are required.")
  }
  n <- nrow(fit$a.binx$data$x)
  if (length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("Could not determine a valid original sample size.")
  }
  # Joint bootstrap covariance for TE, NDE, and total NIE
  effect_names <- c("TE", "NDE", "NIE")
  cov_boot <- cov(boot_mat)
  cov_boot <- (cov_boot + t(cov_boot)) / 2
  Sigma_hat <- n * cov_boot
  Sigma_hat <- (Sigma_hat + t(Sigma_hat)) / 2
  dimnames(cov_boot) <- list(effect_names, effect_names)
  dimnames(Sigma_hat) <- list(effect_names, effect_names)
  # Marginal standard errors and Wald confidence intervals
  SE <- sqrt(pmax(diag(cov_boot), 0))
  names(SE) <- effect_names
  alpha <- 1 - conf_level
  z_crit <- qnorm(1 - alpha / 2)
  results <- data.frame(
    Effect = effect_names,
    Estimate = unname(theta_hat),
    SE = unname(SE),
    Lower_CI = unname(theta_hat - z_crit * SE),
    Upper_CI = unname(theta_hat + z_crit * SE),
    row.names = NULL
  )
  # Mediator-specific indirect effects
  # Exclude the total NIE and the aggregate y1.j1 column
  mediator_names <- setdiff(
    names(ie_hat),
    c("y1.all", "y1.j1")
  )
  missing_boot_cols <- setdiff(
    mediator_names,
    colnames(ie_boot)
  )
  if (length(missing_boot_cols) > 0L) {
    stop(
      "Missing bootstrap indirect-effect columns: ",
      paste(missing_boot_cols, collapse = ", ")
    )
  }
  mediator_results <- lapply(mediator_names, function(med) {
    est <- unname(ie_hat[med])
    b <- as.numeric(ie_boot[, med])
    # Marginal inference uses all finite bootstrap values
    b <- b[is.finite(b)]
    if (length(b) < 2L) {
      stop(
        "Fewer than two valid bootstrap estimates for: ", med
      )
    }
    se <- stats::sd(b)
    data.frame(
      Comparison = NA_character_,
      Effect = paste0("Indirect: ", sub("^y1\\.", "", med)),
      Mediator = sub("^y1\\.", "", med),
      Estimate = est,
      SE = se,
      Lower_CI = est - z_crit * se,
      Upper_CI = est + z_crit * se,
      row.names = NULL
    )
  })
  mediator_results <- if (length(mediator_results)) {
    dplyr::bind_rows(mediator_results)
  } else {
    data.frame(
      Comparison = character(),
      Effect = character(),
      Mediator = character(),
      Estimate = numeric(),
      SE = numeric(),
      Lower_CI = numeric(),
      Upper_CI = numeric()
    )
  }
  # Joint Wald test: NDE = 0 and total NIE = 0
  gamma_hat <- theta_hat[c("NDE", "NIE")]
  gamma_boot <- boot_mat[, c("NDE", "NIE"), drop = FALSE]
  cov_gamma_boot <- cov(gamma_boot)
  cov_gamma_boot <- (
    cov_gamma_boot + t(cov_gamma_boot)
  ) / 2
  V_gamma_hat <- n * cov_gamma_boot
  V_gamma_hat <- (V_gamma_hat + t(V_gamma_hat)) / 2
  dimnames(cov_gamma_boot) <- list(
    c("NDE", "NIE"), c("NDE", "NIE")
  )
  dimnames(V_gamma_hat) <- list(
    c("NDE", "NIE"), c("NDE", "NIE")
  )
  eig_gamma <- eigen(
    cov_gamma_boot,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  tol_gamma <- max(abs(eig_gamma)) * 1e-8
  if (any(!is.finite(eig_gamma)) ||
      max(eig_gamma) <= 0 ||
      min(eig_gamma) <= tol_gamma) {
    stop(
      paste(
        "The NDE-NIE bootstrap covariance matrix is singular",
        "or numerically ill-conditioned."
      )
    )
  }
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
  chi2_crit <- qchisq(conf_level, df = 2)
  joint_region <- list(
    confidence_level = conf_level,
    chi_square_critical_value = chi2_crit,
    center = gamma_hat,
    covariance_estimator_scale = cov_gamma_boot,
    covariance_asymptotic_scale = V_gamma_hat,
    equation = paste0(
      "(gamma_hat - gamma)' inverse(cov_gamma_boot) ",
      "(gamma_hat - gamma) <= ",
      signif(chi2_crit, 6)
    )
  )
  # Joint Wald test for all mediator-specific indirect effects
  mediator_joint_test <- NULL
  mediator_joint_region <- NULL
  cov_mediator_boot <- NULL
  V_mediator_hat <- NULL
  mediator_joint_bootstrap_n <- 0L
  if (length(mediator_names) > 0L) {
    K <- length(mediator_names)
    # Point estimates for all mediator-specific effects
    delta_hat <- unname(ie_hat[mediator_names])
    names(delta_hat) <- mediator_names
    # Preserve paired bootstrap rows across all mediators
    mediator_boot <- ie_boot[
      , mediator_names, drop = FALSE
    ]
    valid_mediators <- complete.cases(mediator_boot) &
      apply(
        mediator_boot,
        1,
        function(z) all(is.finite(z))
      )
    mediator_boot <- mediator_boot[
      valid_mediators, , drop = FALSE
    ]
    mediator_joint_bootstrap_n <- nrow(mediator_boot)
    if (mediator_joint_bootstrap_n < 2L) {
      stop(
        paste(
          "At least two complete bootstrap replicates",
          "are required for mediator joint inference."
        )
      )
    }
    # Full covariance across mediator-specific effects
    cov_mediator_boot <- cov(mediator_boot)
    cov_mediator_boot <- (
      cov_mediator_boot + t(cov_mediator_boot)
    ) / 2
    V_mediator_hat <- n * cov_mediator_boot
    V_mediator_hat <- (
      V_mediator_hat + t(V_mediator_hat)
    ) / 2
    dimnames(cov_mediator_boot) <- list(
      mediator_names, mediator_names
    )
    dimnames(V_mediator_hat) <- list(
      mediator_names, mediator_names
    )
    # Verify full-rank covariance for the ordinary Wald test
    eig_mediator <- eigen(
      cov_mediator_boot,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    tol_mediator <- max(abs(eig_mediator)) * 1e-8
    if (any(!is.finite(eig_mediator)) ||
        max(eig_mediator) <= 0 ||
        min(eig_mediator) <= tol_mediator) {
      stop(
        paste(
          "The mediator-specific bootstrap covariance matrix",
          "is singular or numerically ill-conditioned.",
          "A full-rank joint Wald test cannot be computed."
        )
      )
    }
    # H0: all mediator-specific indirect effects are zero
    mediator_joint_stat <- drop(
      t(delta_hat) %*%
        solve(cov_mediator_boot, delta_hat)
    )
    mediator_joint_p <- pchisq(
      mediator_joint_stat,
      df = K,
      lower.tail = FALSE
    )
    mediator_joint_test <- data.frame(
      Test = paste0(
        "Joint Wald test: all ",
        K,
        " mediator-specific indirect effects = 0"
      ),
      Statistic = mediator_joint_stat,
      DF = K,
      P_value = mediator_joint_p
    )
    # Joint confidence ellipsoid
    mediator_chi2_crit <- qchisq(
      conf_level,
      df = K
    )
    mediator_joint_region <- list(
      confidence_level = conf_level,
      mediator_names = mediator_names,
      center = delta_hat,
      covariance_estimator_scale = cov_mediator_boot,
      covariance_asymptotic_scale = V_mediator_hat,
      chi_square_critical_value = mediator_chi2_crit,
      equation = paste0(
        "(delta_hat - delta)' ",
        "inverse(cov_mediator_boot) ",
        "(delta_hat - delta) <= ",
        signif(mediator_chi2_crit, 6)
      )
    )
  }
  # Joint Wald test: NDE = 0 and every mediator effect = 0
  nde_mediator_joint_test <- NULL
  nde_mediator_joint_region <- NULL
  cov_nde_mediators_boot <- NULL
  V_nde_mediators_hat <- NULL
  nde_mediator_joint_bootstrap_n <- 0L
  if (length(mediator_names) > 0L) {
    # Retain paired replicates for NDE and all mediators
    nde_mediator_boot <- cbind(
      NDE = nde_boot,
      ie_boot[, mediator_names, drop = FALSE]
    )
    valid_nde_mediators <- complete.cases(
      nde_mediator_boot
    ) & apply(
      nde_mediator_boot,
      1,
      function(z) all(is.finite(z))
    )
    nde_mediator_boot <- nde_mediator_boot[
      valid_nde_mediators, , drop = FALSE
    ]
    nde_mediator_joint_bootstrap_n <- nrow(
      nde_mediator_boot
    )
    if (nde_mediator_joint_bootstrap_n < 2L) {
      stop(
        paste(
          "At least two complete bootstrap replicates",
          "are required for joint NDE-mediator inference."
        )
      )
    }
    nde_mediator_names <- c("NDE", mediator_names)
    delta_nde_mediator_hat <- c(
      NDE = nde_hat,
      setNames(
        unname(ie_hat[mediator_names]),
        mediator_names
      )
    )
    cov_nde_mediators_boot <- cov(nde_mediator_boot)
    cov_nde_mediators_boot <- (
      cov_nde_mediators_boot +
        t(cov_nde_mediators_boot)
    ) / 2
    V_nde_mediators_hat <- n * cov_nde_mediators_boot
    V_nde_mediators_hat <- (
      V_nde_mediators_hat +
        t(V_nde_mediators_hat)
    ) / 2
    dimnames(cov_nde_mediators_boot) <- list(
      nde_mediator_names,
      nde_mediator_names
    )
    dimnames(V_nde_mediators_hat) <- list(
      nde_mediator_names,
      nde_mediator_names
    )
    eig_nde_mediators <- eigen(
      cov_nde_mediators_boot,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    tol_nde_mediators <-
      max(abs(eig_nde_mediators)) * 1e-8
    if (any(!is.finite(eig_nde_mediators)) ||
        max(eig_nde_mediators) <= 0 ||
        min(eig_nde_mediators) <= tol_nde_mediators) {
      stop(
        paste(
          "The joint NDE-mediator bootstrap covariance",
          "matrix is singular or numerically ill-conditioned."
        )
      )
    }
    # H0: NDE = 0 and all mediator-specific indirect effects = 0
    nde_mediator_joint_stat <- drop(
      t(delta_nde_mediator_hat) %*%
        solve(
          cov_nde_mediators_boot,
          delta_nde_mediator_hat
        )
    )
    nde_mediator_df <- length(nde_mediator_names)
    nde_mediator_joint_p <- pchisq(
      nde_mediator_joint_stat,
      df = nde_mediator_df,
      lower.tail = FALSE
    )
    nde_mediator_joint_test <- data.frame(
      Test = paste0(
        "Joint Wald test: NDE = 0 and all ",
        length(mediator_names),
        " mediator-specific indirect effects = 0"
      ),
      Statistic = nde_mediator_joint_stat,
      DF = nde_mediator_df,
      P_value = nde_mediator_joint_p
    )
    nde_mediator_chi2_crit <- qchisq(
      conf_level,
      df = nde_mediator_df
    )
    nde_mediator_joint_region <- list(
      confidence_level = conf_level,
      parameter_names = nde_mediator_names,
      center = delta_nde_mediator_hat,
      covariance_estimator_scale = cov_nde_mediators_boot,
      covariance_asymptotic_scale = V_nde_mediators_hat,
      chi_square_critical_value = nde_mediator_chi2_crit,
      equation = paste0(
        "(delta_hat - delta)' ",
        "inverse(cov_nde_mediators_boot) ",
        "(delta_hat - delta) <= ",
        signif(nde_mediator_chi2_crit, 6)
      )
    )
  }
  # Decomposition diagnostic: TE = NDE + NIE
  point_decomposition_error <- unname(
    theta_hat["TE"] -
      theta_hat["NDE"] -
      theta_hat["NIE"]
  )
  bootstrap_decomposition_error <-
    boot_mat[, "TE"] -
    boot_mat[, "NDE"] -
    boot_mat[, "NIE"]
  decomposition_diagnostic <- list(
    point_estimate_error = point_decomposition_error,
    max_absolute_bootstrap_error = max(
      abs(bootstrap_decomposition_error)
    ),
    bootstrap_errors = bootstrap_decomposition_error
  )
  # Estimated rank of TE-NDE-NIE covariance
  eig_full <- eigen(
    cov_boot,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  tolerance <- max(abs(eig_full)) * 1e-8
  estimated_rank <- sum(eig_full > tolerance)
  # Return all results
  list(
    # Overall mediation effects
    results = results,
    theta_hat = theta_hat,
    # Mediator-specific marginal results
    mediator_results = mediator_results,
    # Overall bootstrap covariance
    bootstrap_covariance = cov_boot,
    Sigma_hat = Sigma_hat,
    # Joint NDE and total NIE inference
    gamma_hat = gamma_hat,
    covariance_NDE_NIE = cov_gamma_boot,
    V_gamma_hat = V_gamma_hat,
    joint_wald_test = joint_test,
    joint_confidence_region = joint_region,
    # Joint mediator-specific indirect-effect inference
    mediator_joint_wald_test = mediator_joint_test,
    mediator_joint_confidence_region = mediator_joint_region,
    covariance_mediator_effects = cov_mediator_boot,
    V_mediator_hat = V_mediator_hat,
    mediator_joint_bootstrap_n = mediator_joint_bootstrap_n,
    # Joint NDE and all mediator-specific indirect effects
    nde_mediator_joint_wald_test = nde_mediator_joint_test,
    nde_mediator_joint_confidence_region =
      nde_mediator_joint_region,
    covariance_NDE_mediator_effects =
      cov_nde_mediators_boot,
    V_NDE_mediator_hat = V_nde_mediators_hat,
    nde_mediator_joint_bootstrap_n =
      nde_mediator_joint_bootstrap_n,
    # Diagnostics and metadata
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

