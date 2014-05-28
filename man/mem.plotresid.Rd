\name{mem.plotresid}
\alias{mem.plotresid}
\title{
Plot residuals of a mixed-effects model along with trend lines
}
\description{
This function plots stylized residuals of a mixed-effects model. It is possible to obtain fitted values versus errors (XbZu vs e), or original values versus errors (y vs e) in order to obtain different views to the errors in connection to the observations.
}
\usage{
mem.plotresid(fit, linear = T, type = "XbZu", main, xlab, ylab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
A fitted mixed-effects model generated either through the lme4 or the nlme package.
}
  \item{linear}{
Should linear trend lines be drawn to the residual plot
}
  \item{type}{
Type of residual plot; should fitted values (value "XbZu") or original observations (value "y") be plotted against "e" errors
}
  \item{main}{
Main title
}
  \item{xlab}{
x-axis label
}
  \item{ylab}{
y-axis label
}
}
\details{
Notice that the typical residual plot is fitted values (type="XbZu") versus errors ("e").
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{mem.getcomp}}, \code{\link{mem.plotran}}
}
\examples{
data(vcaplong)

exdat <- vcaplong[vcaplong[,"Group"] \%in\% c("Vehicle", "ARN"),]

library(lme4)
f0 <- lmer(log2PSA ~ 1 + DrugWeek + (1 + DrugWeek|ID), data = exdat)
f1 <- lmer(log2PSA ~ 1 + DrugWeek + DrugWeek:ARN + (1 + DrugWeek|ID), data = exdat)
f2 <- lmer(log2PSA ~ 1 + DrugWeek + DrugWeek:ARN + (1|ID) + (0 + DrugWeek|ID), data = exdat)
f3 <- lmer(log2PSA ~ 1 + DrugWeek + DrugWeek:ARN + (1|Submatch) + (0 + DrugWeek|ID), data = exdat)

par(mfrow=c(2,2))
mem.plotresid(f0)
mem.plotresid(f1)
mem.plotresid(f2)
mem.plotresid(f3)

}
\keyword{ regression }
