\name{jmmadata}
\docType{data}
\alias{jmmadata}
\title{Simulated competing-risks data for mediation analysis}
\description{
This simulated dataset contains baseline characteristics, binary mediators, treatment group assignments, and competing-risks event-time outcomes for 1200 individuals. It was generated to demonstrate mediation analysis with competing risks and restricted mean time estimands.

The treatment variable contains three categories: Control, Moderate, and Intensive. A binary treatment variable, \code{treatment2}, can also be derived by combining the Moderate and Intensive groups into an Intervention category and comparing it with Control.

Two binary mediators are included: high inflammation and poor functional status. The competing events are readmission and death before readmission, with administrative or loss-to-follow-up censoring.

The dataset is simulated using prespecified probabilistic models and is intended for methodological illustration rather than clinical inference.
}
\usage{jmmadata}
\format{
A data frame with 1200 observations and 10 variables:
\itemize{
  \item \code{id}: Unique participant identifier.
  \item \code{treatment}: Three-level factor indicating treatment group: Control, Moderate, or Intensive.
  \item \code{age}: Age in years, simulated from a truncated normal distribution between 40 and 90.
  \item \code{sex}: Factor indicating sex, with levels Female and Male.
  \item \code{comorbidity}: Simulated count of comorbidities, generated from a Poisson distribution.
  \item \code{inflammation_high}: Binary indicator of high inflammation (0 = no, 1 = yes).
  \item \code{poor_function}: Binary indicator of poor functional status (0 = no, 1 = yes).
  \item \code{GF_time_cr}: Observed follow-up time until the first competing event or censoring.
  \item \code{GF_status_cr}: Event status (0 = censored, 1 = readmission, 2 = death before readmission).
  \item \code{treatment2}: Binary treatment grouping derived from \code{treatment} (Control versus Intervention, where Intervention combines Moderate and Intensive).
}
}
\details{
The dataset is generated using simulated cause-specific event hazards.

The readmission hazard depends on treatment, inflammation, functional status, age, and comorbidity. The death hazard depends on treatment, age, comorbidity, and functional status. Event times are generated from exponential distributions, and censoring times are generated independently from a uniform distribution between 3 and 8 time units.

The variable \code{GF_status_cr} represents the observed first-event type:
\itemize{
  \item 0: Censored before either event.
  \item 1: Readmission occurs before death and before censoring.
  \item 2: Death occurs before readmission and before censoring.
}

The binary variable \code{treatment2} is created by combining the Moderate and Intensive treatment groups. It is useful for comparing the combined Intervention group with the Control group in binary-exposure mediation analyses.
}
\seealso{
\code{\link{cuminc}},
\code{\link{mma}}
}
\examples{
\dontrun{
data(jmmadata)

head(jmmadata)

# Examine treatment group sizes
table(jmmadata$treatment)

# Examine the binary treatment grouping
table(jmmadata$treatment2)

# Examine competing event counts
table(jmmadata$GF_status_cr)

# Examine event status by treatment group
table(jmmadata$treatment, jmmadata$GF_status_cr)

# Summarize observed follow-up time
summary(jmmadata$GF_time_cr)

# Examine binary mediator prevalence
prop.table(table(jmmadata$inflammation_high))
prop.table(table(jmmadata$poor_function))
}
}
\keyword{datasets}
