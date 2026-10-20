\name{plot.Jmma}
\alias{plot.Jmma}
\title{Plot Joint Mediation Analysis Results}
\description{
Creates graphical summaries of the results from a joint mediation analysis
fitted using \code{\link{jmma}}. Two plot types are supported: a forest
plot of total, natural direct, and natural indirect effect estimates with
confidence intervals, and joint confidence ellipse plots for the natural
direct and natural indirect effects.
}
\usage{
\method{plot}{Jmma}(x,
     type = c("forest", "joint"),
     comparison = NULL,
     conf_level = NULL,
     digits = 3,
     ...)
}
\arguments{
  \item{x}{
    An object inheriting from class \code{"Jmma"}, typically returned by
    \code{\link{jmma}}.
  }
  \item{type}{
    Character string specifying the type of plot. The default,
    \code{"forest"}, produces a forest plot of the marginal effect
    estimates and their confidence intervals. The alternative,
    \code{"joint"}, produces joint confidence ellipse plots for the natural
    direct and natural indirect effects.
  }
  \item{comparison}{
    Optional character vector specifying which exposure comparisons to
    display. Comparison names must match those in the fitted
    \code{"Jmma"} object.

    For \code{type = "forest"}, the plot includes the total effect (TE),
    natural direct effect (NDE), and natural indirect effect (NIE) for
    each selected comparison. If \code{NULL}, all available comparisons
    are included.

    For \code{type = "joint"}, each selected comparison is displayed as a
    separate confidence ellipse. If \code{NULL}, all available comparisons
    are plotted in separate panels. A single comparison produces one plot.
  }
  \item{conf_level}{
    Numeric scalar specifying the confidence level for the plotted
    confidence intervals or joint confidence regions. Must be strictly
    between 0 and 1. If \code{NULL}, the confidence level specified in the
    fitted object is used.
  }
  \item{digits}{
    Integer specifying the number of digits for numerical display.
    Defaults to \code{3}. This argument is currently not used by the
    plotting implementation.
  }
  \item{...}{
    Additional graphical parameters passed to the underlying plotting
    functions, such as \code{cex}, \code{cex.main}, or \code{cex.lab}.
  }
}
\details{
When \code{type = "forest"}, the method displays point estimates and
confidence intervals for the TE, NDE, and NIE for each selected exposure
comparison. A vertical reference line at zero is included to facilitate
interpretation of the estimated effects.

When \code{type = "joint"}, the method displays a two-dimensional Wald-type
joint confidence ellipse for the NDE and NIE for each selected comparison.
The horizontal axis represents the NDE, and the vertical axis represents
the NIE. The estimated effects are marked by a point, and reference lines
at zero are drawn on both axes.

The joint confidence ellipse is constructed using the estimated NDE and NIE
and their bootstrap covariance matrix on the estimator scale. Its
boundary is determined using the chi-squared quantile with two degrees of
freedom corresponding to the specified confidence level.

If the joint confidence ellipse excludes the origin, the joint null
hypothesis
\deqn{H_0: \mathrm{NDE} = 0 \quad \text{and} \quad \mathrm{NIE} = 0}
is rejected at the corresponding significance level under the Wald
procedure. This indicates evidence that at least one of the two effects
differs from zero; it does not establish that both effects are nonzero.

When multiple comparisons are selected for a joint plot, the method
arranges the plots in a panel layout. Each panel compares the specified
exposure category with the reference category used in the fitted model.
The plots do not directly test differences in mediation effects between
exposure categories.

The method checks that the estimates and covariance matrices required
for the joint confidence ellipses are valid and that the covariance
matrices are positive definite. An error is raised if these requirements
are not met.

The user's graphical parameters are restored after plotting.
}
\value{
The method returns \code{x} invisibly, leaving the plot on the active
graphics device.
}
\seealso{
\code{\link{jmma}},
\code{\link{summary.Jmma}}
}
\examples{
\dontrun{

## Forest plot of all marginal effect estimates.
plot(fit, type = "forest")

## Forest plot for a selected comparison.
plot(fit, type = "forest", comparison = "Moderate")

## Joint NDE-NIE confidence ellipses for all comparisons.
plot(fit, type = "joint")

## Joint confidence ellipse for a single comparison.
plot(fit, type = "joint", comparison = "Intensive")

## Joint confidence ellipses at the 90% confidence level.
plot(fit, type = "joint", conf_level = 0.90)

## Joint confidence ellipses for selected comparisons.
plot(fit, type = "joint",  comparison = c("Moderate", "Intensive"))

}
}
