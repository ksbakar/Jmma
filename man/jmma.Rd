\name{jmma}
\alias{jmma}
\title{Joint Multiple Mediation Analysis (JMMA): Competing-Risks RMST using Pseudo-Values and G-Computation}
\description{
Performs joint multiple mediation analysis and provides marginal and joint statistical inference for the total effect (TE), natural direct effect (NDE), and natural indirect effect (NIE).

The function supports binary and categorical exposure variables. For categorical exposures, each specified comparison category is analyzed separately against a reference category.

The outcome may be supplied directly or constructed as pseudo-values of the restricted mean time free of a specified competing event, estimated using the cumulative incidence function. Bootstrap estimates from the fitted mediation models are used to estimate standard errors, confidence intervals, and joint Wald tests for the natural direct and indirect effects.
}
\usage{
jmma(data, exposure, comparison = NULL, reference, covariates, mediators, outcome,
  n_boot = 10, conf_level = 0.95, use_pseudo = TRUE, tau = NULL, event_code = 1,
  time_var, status_var, ...)
}
\arguments{
  \item{data}{
A data frame containing the exposure, covariates, mediators, and outcome variables, or the follow-up time and event status variables when pseudo-values are calculated.
}

  \item{exposure}{
A character string specifying the name of the exposure variable in \code{data}. The exposure may have two or more categories.
}

  \item{comparison}{
Optional character vector specifying the exposure categories to compare against \code{reference}. Each comparison is fitted separately using the comparison category and the reference category. If \code{NULL}, all observed exposure categories other than the reference category are used.
}

  \item{reference}{
A character string specifying the reference category of the exposure. This category must be present in the analysis data.
}

  \item{covariates}{
A character vector containing the names of baseline covariates to adjust for in the mediation analysis. These variables must be present in \code{data} and must not overlap with the specified mediators.
}

  \item{mediators}{
A character vector containing the names of the potential mediators in \code{data}. The current implementation specifies these variables as binary mediators in \code{mma}; therefore, they should be coded as 0/1 variables with 0 as the reference category.
}

  \item{outcome}{
A character string specifying the name of the outcome variable. If \code{use_pseudo = TRUE}, this name is assigned to the generated RMST pseudo-values. Otherwise, the variable must already exist in \code{data}.
}

  \item{n_boot}{
The number of bootstrap resamples used by \code{mma} to estimate mediation effects. Passed to the \code{n2} argument of \code{mma}. The default is 100.
}

  \item{conf_level}{
The confidence level used to construct marginal Wald confidence intervals and the joint Wald confidence region for the NDE and NIE. The default is 0.95.
}

  \item{use_pseudo}{
Logical. If \code{TRUE}, calculate pseudo-values for the restricted mean time free of the event of interest under a competing-risks framework and use them as the outcome. If \code{FALSE}, use the outcome variable supplied in \code{data}. The default is \code{FALSE}.
}

  \item{tau}{
A positive numeric restriction time for calculating the restricted mean time estimand. Required when \code{use_pseudo = TRUE}.
}

  \item{event_code}{
An integer specifying the event of interest in the competing-risks status variable. The default is 1. For the example dataset, 1 represents readmission and 2 represents death before readmission.
}

  \item{time_var}{
A character string specifying the name of the observed follow-up time variable. Used when \code{use_pseudo = TRUE}. The default is \code{"GF_time_cr"}.
}

  \item{status_var}{
A character string specifying the name of the competing-risks event status variable. Used when \code{use_pseudo = TRUE}. The default is \code{"GF_status_cr"}.
}

  \item{...}{
Additional arguments passed to \code{mma::mma()}. See \code{\link[ mma]{mma}} for details.
}
}
\details{
\strong{Data preparation}

The function selects the required variables and removes observations with missing values in those variables. When \code{use_pseudo = TRUE}, the RMST pseudo-values are calculated after this complete-case preparation.

\strong{Exposure comparisons}

For each comparison category, the function retains observations in the reference and comparison groups and fits a separate binary-exposure mediation model. The comparison category is coded as 1 and the reference category as 0. Consequently, each fitted model estimates effects for its specified comparison relative to the reference group.

\strong{Competing-risks RMST pseudo-values}

When \code{use_pseudo = TRUE}, the function estimates the restricted mean time free of the event of interest up to \code{tau} using one minus the cumulative incidence function. It then calculates jackknife pseudo-values from the full-sample estimate and leave-one-out estimates. These pseudo-values are used as the outcome in \code{mma}.

This estimand is not the conventional all-cause survival RMST. It describes the area under one minus the cumulative incidence function for the selected event and therefore concerns the probability of remaining free of that particular event in the presence of competing events.

The pseudo-values are calculated once before the pairwise mediation models are fitted. The standard bootstrap procedure within \code{mma} does not automatically recalculate the competing-risks cumulative incidence estimator and its pseudo-values within every bootstrap resample.

\strong{Mediation inference}

For each fitted model, the function extracts the estimated TE, NDE, and NIE and their bootstrap estimates. It calculates bootstrap standard errors and marginal Wald confidence intervals. It also performs a two-degree-of-freedom joint Wald test of the null hypothesis that both the NDE and NIE are zero.

The joint confidence region for the NDE and NIE is represented by an ellipse based on their estimated bootstrap covariance matrix. A diagnostic is returned to assess whether the estimated total effect is approximately equal to the sum of the direct and indirect effects.
}
\value{
A list containing the following components:
\itemize{
  \item \code{data}: The prepared analysis dataset. When \code{use_pseudo = TRUE}, this includes the generated pseudo-value outcome.
  \item \code{inference}: A named list of inference results for each fitted model. Each element contains marginal estimates and confidence intervals, bootstrap covariance estimates, the joint Wald test, the joint confidence region, and diagnostic information.
  \item \code{results}: A data frame combining the TE, NDE, and NIE estimates, standard errors, and marginal confidence intervals across comparisons.
  \item \code{joint_tests}: A data frame containing the joint Wald test results for the NDE and NIE for each comparison.
  \item \code{settings}: A list recording the analysis settings, including the exposure, reference category, comparisons, covariates, mediators, outcome, pseudo-value settings, model options, bootstrap settings, and analysis sample size.
}
}
\references{

}
\seealso{
\code{\link[mma]{mma}},
\code{\link[cmprsk]{cuminc}},
\code{\link{summary.Jmma}},
\code{\link{plot.Jmma}}
}
\examples{
\dontrun{

head(jmmadata)

## Binary exposure: Intervention versus Control

res_binary <- jmma(
  data = jmmadata,
  exposure = "treatment2",
  comparison = "Intervention",
  reference = "Control",
  covariates = c("age", "sex", "comorbidity"),
  mediators = c("inflammation_high", "poor_function"),
  outcome = "RMST_pseudo",
  use_pseudo = TRUE,
  tau = 5,
  event_code = 1,
  time_var = "GF_time_cr",
  status_var = "GF_status_cr",
  n_boot = 10,
  conf_level = 0.95
)
res_binary
summary(res_binary)
plot(res_binary)

## Three-category exposure: Moderate and Intensive versus Control

res_multi <- jmma(
  data = simdata,
  exposure = "treatment",
  comparison = c("Moderate", "Intensive"),
  reference = "Control",
  covariates = c("age", "sex", "comorbidity"),
  mediators = c("inflammation_high", "poor_function"),
  outcome = "RMST_pseudo",
  use_pseudo = TRUE,
  tau = 5,
  event_code = 1,
  time_var = "GF_time_cr",
  status_var = "GF_status_cr",
  n_boot = 10,
  conf_level = 0.95
)

res_binary
summary(res_binary)
plot(res_binary)

## TE, NDE, NIE and marginal confidence intervals
res_multi$results

## Joint tests of NDE = 0 and NIE = 0
res_multi$joint_tests

## Intensive versus Control
res_multi$inference[["Intensive"]]$results

## Moderate versus Control
res_multi$inference[["Moderate"]]$results

## Bootstrap covariance among TE, NDE and NIE
res_multi$inference[["Intensive"]]$bootstrap_covariance

## Mediation decomposition diagnostic
res_multi$inference[["Intensive"]]$decomposition_diagnostic
}
}
\keyword{models}
