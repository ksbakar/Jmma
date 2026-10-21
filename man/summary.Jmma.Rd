\name{summary.Jmma}
\alias{summary.Jmma}
\alias{print.summary.Jmma}
\title{Summarize Joint Mediation Analysis Results}
\description{
Summarizes the results of a joint mediation analysis fitted using \code{\link{jmma}}. The summary includes model settings, marginal effect
estimates, joint Wald tests for the natural direct and indirect effects, and bootstrap and effect-decomposition diagnostics. Optionally, covariance
matrices can be retained for display by the associated print method.
}
\usage{
\method{summary}{Jmma}(object, cov_matrix = FALSE, diagnostics = FALSE, ...)

\method{print}{summary.Jmma}(x, digits = 4, ...)
}
\arguments{
  \item{object}{
    An object inheriting from class \code{"Jmma"}, typically returned by
    \code{\link{jmma}}.
  }
  \item{cov_matrix}{
    Logical. If \code{TRUE}, the print method displays the bootstrap
    covariance matrix of the total, natural direct, and natural indirect
    effect estimates, the bootstrap covariance matrix of the natural
    direct and indirect effects, and their corresponding asymptotic
    covariance matrices. Defaults to \code{FALSE}.
  }
  \item{diagnostics}{
    Logical. If \code{TRUE}, the print method displays the bootstrap
    and decomposition diagnostics. Defaults to \code{FALSE}.
  }
  \item{x}{
    An object of class \code{"summary.Jmma"}, returned by
    \code{summary.Jmma}.
  }
  \item{digits}{
    Integer. The number of significant digits used when printing numerical
    results. Defaults to \code{4}.
  }
  \item{...}{
    Additional arguments passed to or from methods. Currently not used
    by these methods.
  }
}
\details{
The summary object contains the model settings, marginal effect estimates,
joint Wald test results, diagnostics, and covariance matrices for each
exposure comparison.

The marginal effect estimates include the total effect (TE), natural direct
effect (NDE), and natural indirect effect (NIE), along with their associated
standard errors and confidence intervals, as available in the fitted
object.

For each comparison, the joint Wald test evaluates the null hypothesis
\deqn{H_0: \mathrm{NDE} = 0 \quad \text{and} \quad \mathrm{NIE} = 0}
against the alternative hypothesis
\deqn{H_1: \mathrm{NDE} \ne 0 \quad \text{or} \quad \mathrm{NIE} \ne 0.}
A statistically significant result provides evidence that at least one
of the two effects differs from zero; it does not necessarily imply that
both effects are nonzero.

The diagnostics report the number of bootstrap replicates, analysis sample
size, estimated covariance matrix rank, the point-estimate discrepancy
\eqn{\mathrm{TE} - \mathrm{NDE} - \mathrm{NIE}}, and the maximum absolute
bootstrap decomposition error.

Covariance matrices are stored in the summary object regardless of the
value of \code{cov_matrix}. Setting \code{cov_matrix = TRUE} controls
whether the associated print method displays them; it does not control
whether they are computed or retained.
}
\value{
\code{summary.Jmma} returns an object of class \code{"summary.Jmma"}, a
list containing the following components:
\describe{
  \item{\code{settings}}{Model settings, including the exposure, reference
  category, outcome, analysis sample size, and confidence level.}
  \item{\code{results}}{Marginal effect estimates and their associated
  inferential quantities.}
  \item{\code{joint_tests}}{Joint Wald test results for the natural direct
  and natural indirect effects for each exposure comparison.}
  \item{\code{diagnostics}}{A data frame containing bootstrap and effect-
  decomposition diagnostics for each comparison.}
  \item{\code{bootstrap_covariance}}{A named list of bootstrap covariance
  matrices for the TE, NDE, and NIE estimates, one per comparison.}
  \item{\code{covariance_NDE_NIE}}{A named list of bootstrap covariance
  matrices for the NDE and NIE estimates, one per comparison.}
  \item{\code{Sigma_hat}}{A named list of asymptotic covariance matrices
  for the TE, NDE, and NIE estimates, one per comparison.}
  \item{\code{V_gamma_hat}}{A named list of asymptotic covariance matrices
  for the NDE and NIE estimates, one per comparison.}
  \item{\code{inference}}{The comparison-specific inference objects from
  the original \code{"Jmma"} object.}
  \item{\code{cov_matrix}}{The logical value specifying whether covariance
  matrices should be displayed by the print method.}
}
The \code{print.summary.Jmma} method prints the summary and invisibly
returns \code{x}.
}
\seealso{
\code{\link{jmma}},
\code{\link{plot.Jmma}}
}
\examples{
\dontrun{

## Summarize the results without displaying covariance matrices.
s <- summary(fit)
print(s)

## Display the covariance matrices as well.
print(summary(fit, cov_matrix = TRUE))

## Suppress covariance matrix output explicitly.
print(summary(fit, cov_matrix = FALSE), digits = 3)

}
}
