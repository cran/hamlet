\name{mem.getcomp}
\alias{mem.getcomp}
\title{
Extract per-observation components for fixed and random effects of a mixed-effects model
}
\description{
Assuming a mixed-effects model of form y_fit = Xb + Zu + e, where X is the model matrix for fixed effects, Z is the model matrix for random effects, and b and u are the fixed and random effects respectively, this function returns these components per each fitted value y. These may be useful for model inference and/or diagnostic purposes.
}
\usage{
mem.getcomp(fit)
}
\arguments{
  \item{fit}{
A fitted mixed-effects model generated either through the lme4 or the nlme package.
}
}
\details{
Notice that per-observation model error is e = Xb + Zu - y_observation. Similarly, the components Xb and Zu are extracted.
}
\value{
The function returns per-observation model fit components for a mixed-effects model. The return fields are
\item{Xb}{Fixed effects component Xb}
\item{Zu}{Random effects component Zu}
\item{XbZu}{Full model fit by summing the above two Xb+Zu}
\item{e}{Model error e}
\item{y}{Original observations y}
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{mem.plotran}}
\code{\link{mem.plotresid}}
}
\examples{
data(vcaplong)

exdat <- vcaplong[vcaplong[,"Group"] \%in\% c("Vehicle", "ARN"),]

library(lme4)
f1 <- lmer(log2PSA ~ 1 + DrugWeek + DrugWeek:ARN + (1 + DrugWeek|ID), data = exdat)

mem.getcomp(f1)
}
\keyword{ regression }
