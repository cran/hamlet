\name{orxwide}
\alias{orxwide}
\docType{data}
\title{
Wide-format baseline data for the ORX study
}
\description{
This data frame contains the wide-format data of the ORX study for baseline characteristics of the 
individuals participating in the study. Some fields (Volume, PSA, High, BodyWeight, PSAChange) were 
used to construct the distance matrix in the original matching-based random allocation of individuals
at baseline, while other variables (Group, Submatch) contain these results.
}
\usage{data("orxwide")}
\format{
  A data frame with 109 observations on the following 8 variables.
  \describe{
    \item{\code{ID}}{A unique character indicator for the different individual(s)}
    \item{\code{Group}}{After identifying suitable submatches, the data were distributed to blinded intervention groups. These groups were later then annotated to actual treatments or non-intervention control groups.}
    \item{\code{Submatch}}{Submatches identified at baseline using the methodology presented in this package}
    \item{\code{Volume}}{Tumor volume at baseline in cubic millimeters}
    \item{\code{PSA}}{Raw baseline PSA measurement values in unit (ug/l)}
    \item{\code{High}}{The highest dimension in the tumor in millimeters, giving insight into the shape of the tumor}
    \item{\code{BodyWeight}}{Body weight at baseline in unit (g)}
    \item{\code{PSAChange}}{A fold-change like change in PSA from the prior measurement defined as: (PSA_current - PSA_last)/(PSA_last)}
  }
}
\note{
Notice that while normally the submatches would be distributed equally to the experiment 
groups, here rarely a single submatch may hold multiple instances from a single group. 
This is due to practical constraints in the experiment, that animals had to be manually
moved in order to fulfill groups and to reflect the amount of drug compounds available.
Additionally, the original experiment was performed on 6 intervention groups, while here
only 3 are further presented after the baseline ('ORX+Tx', 'ORX' and 'Intact').
}
\details{
Originally, 3-fold weighting of the baseline 'Volume' and 'PSA' was used in comparison to 
'High', 'BodyWeight' and 'PSAChange' when computing the distance matrix. Furthermore, some 
individuals were annotated prior to matching for exclusion based on outlierish behaviour.
The exclusion criteria were applied before any interventions were given or the matching was 
performed. The excluded tumors had either non-existant PSA, non-detectable tumor volume,
or extremely large tumors (volume above 700 mm^3).
}
\source{
Huhtaniemi R, et al. Manuscript in preparation.
}
\examples{
data(orxwide)
# Construct an example distance matrix based on conventional 
# Euclidean distance and the baseline characteristics
d.orx <- dist(orxwide[,c("Volume", "PSA", "High", "BodyWeight", "PSAChange")])
# Plot a hierarchical clustering of the individuals
plot(hclust(d=d.orx))
# This 'd.orx' may then be further processed by downstream experiment 
# design functions such as match.ga, match.bb, etc.
}
\keyword{datasets}
