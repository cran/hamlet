\name{mem.plotran}
\alias{mem.plotran}
\title{
Plot random effects histograms for a fitted mixed-effects model
}
\description{
This plot creates histogram plots for the columns extracted from random effects from a model fit. This is useful for model diagnostics, such as observing deviations from normality in the random effects.
}
\usage{
mem.plotran(fit, breaks = 100)
}
\arguments{
  \item{fit}{
A fitted mixed-effects model generated either through the lme4 or the nlme package.
}
  \item{breaks}{
Number of breaks in the histograms (passed to the 'hist'-function)
}
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{mem.getcomp}}, \code{\link{mem.plotresid}}
}
\examples{
data(vcaplong)

exdat <- vcaplong[vcaplong[,"Group"] \%in\% c("Vehicle", "ARN"),]

library(lme4)
f1 <- lmer(log2PSA ~ 1 + DrugWeek + DrugWeek:ARN + (1 + DrugWeek|ID), data = exdat)

ranef(f1) # Histograms are plotted for these columns
mem.plotran(f1)
}
\keyword{ regression }
