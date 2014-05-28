\name{vcaplong}
\alias{vcaplong}
\docType{data}
\title{
Long-format data of the Castration-resistant Prostate Cancer experiment using the VCaP cell line.
}
\description{
The long-format of the VCaP experiment PSA-measurements may be used to model longitudinal measurements during interventions (Vehicle, ARN, or MDV). Body weights and PSA were measured weekly during the experiment. PSA concentrations were log2-transformed to make data better normally distributed.
}
\usage{data(vcaplong)}
\format{
  A data frame with 225 observations on the following 10 variables.
  \describe{
    \item{\code{log2PSA}}{Log2-transformed PSA (prostate-specific antigen) measurements with unit (log2 ug/l)}
    \item{\code{BW}}{Body weights (g)}
    \item{\code{Submatch}}{A grouping factor for indicating which measurements belong to individuals that were part of the same submatch prior to interventions}
    \item{\code{ID}}{A character vector indicating unique animal IDs}
    \item{\code{Week}}{Week of the experiment, notice that this is not the same as the week of drug administration (see below)}
    \item{\code{DrugWeek}}{Week since beginning administration of the drugs}
    \item{\code{Group}}{Grouping factor for intervention groups of the observations}
    \item{\code{Vehicle}}{Binary indicator for which observations belonged to the group 'Vehicle'}
    \item{\code{ARN}}{Binary indicator for which observations belonged to the group 'ARN-509'}
    \item{\code{MDV}}{Binary indicator for which observations belonged to the group 'MDV3100'}
  }
}
\details{
Notice that the long-format is suitable for modeling longitudinal measurements. The grouping factors ID or Submatch could be used to group observations belonging to a single individual or matched individuals.
}
\source{
Knuuttila M, Yatkin E, Kallio J, Savolainen S, Laajala TD, et al. (2014) Castration induces upregulation of intratumoral androgen biosynthesis and androgen receptor expression in orthotopic VCaP human prostate cancer xenograft model.  Am J Pathol. Accepted
}
\references{
Laajala TD, Kaur A, Knuuttila M, Westermarck J, et al. (2014) Optimal Matching, Randomization and Statistical Modeling of Hierarchical Baseline Variables in Preclinical Cancer Studies. Submitted
}
\examples{
data(vcaplong)

str(vcaplong)
head(vcaplong)

library(lattice)
xyplot(log2PSA ~ DrugWeek | Group, data = vcaplong, type="l", group=ID, layout=c(3,1))
xyplot(BW ~ DrugWeek | Group, data = vcaplong, type="l", group=ID, layout=c(3,1))
}
\keyword{datasets}
