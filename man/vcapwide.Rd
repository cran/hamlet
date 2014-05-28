\name{vcapwide}
\alias{vcapwide}
\docType{data}
\title{
Wide-format data of the Castration-resistant Prostate Cancer experiment using the VCaP cell line.
}
\description{
VCaP cancer cells were injected orthotopically into the prostate of mice and PSA (prostate-specific antigen) was followed. The animals were castrated on two subsequent weeks, after which the castration-resistant tumors were allowed to emerge. Since PSA reached pre-castration levels, the animals were non-bipartite matched and allocated to separate intervention arms (at week 10). 3 different interventions are presented here, with 'Vehicle' as a comparison point and MDV3100 and ARN-509 tested for reducing PSA and its correlated tumor size.
}
\usage{data(vcapwide)}
\format{
  A data frame with 45 observations on the following 34 variables.
  \describe{
    \item{\code{CastrationDate}}{A numeric vector indicating week when the animal was castrated, resulting in steep decrease in PSA and subsequent castration-resistant tumors to emerge.}
    \item{\code{CageAtAllocation}}{A factorial vector indicating cage labels for each animal at the intervention allocation.}
    \item{\code{Group}}{A character vector indicating which intervention group the animal was allocated to in the actual experiment (3 alternatives).}
    \item{\code{TreatmentInitiationWeek}}{A character vector indicating at which week the intervention was started.}
    \item{\code{Submatch}}{A character vector indicating which submatch the individual was part of the original non-bipartite matching task.}
    \item{\code{ID}}{A unique character vector indicating the animals.}
    \item{\code{PSAWeek2}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek3}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek4}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek5}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek6}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek7}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek8}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek9}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek10}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek11}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek12}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek13}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{PSAWeek14}}{Numeric vector(s) indicating PSA concentration (ug/l) per each week (2 to 14) of the experiment.}
    \item{\code{BWWeek0}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek1}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek2}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek3}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek4}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek5}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek6}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek7}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek8}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek9}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek10}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek11}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek12}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek13}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
    \item{\code{BWWeek14}}{Numeric vector indicating body weight (g) of the animals per each week (0 to 14) of the experiment.}
  }
}
\details{
The wide-format here presented the longitudinal measurements for PSA and Body Weight per each column. For modeling the PSA growth longitudinally e.g. using mixed-effects models, see the \code{\link{vcaplong}} dataset where the data has been readily transposed into the long-format.
}
\source{
Knuuttila M, Yatkin E, Kallio J, Savolainen S, Laajala TD, et al. (2014) Castration induces upregulation of intratumoral androgen biosynthesis and androgen receptor expression in orthotopic VCaP human prostate cancer xenograft model.  Am J Pathol. Accepted
}
\references{
Laajala TD, Kaur A, Knuuttila M, Westermarck J, et al. (2014) Optimal Matching, Randomization and Statistical Modeling of Hierarchical Baseline Variables in Preclinical Cancer Studies. Submitted
}
\seealso{
\code{\link{vcaplong}}
}
\examples{
data(vcapwide)

str(vcapwide)
head(vcapwide)

mixplot(vcapwide[,c("PSAWeek10", "PSAWeek14", "BWWeek10", "Group")], pch=16)
anv <- aov(PSA ~ Group, data.frame(PSA = vcapwide[,"PSAWeek14"], Group = vcapwide[,"Group"]))
summary(anv)
TukeyHSD(anv)
summary(aov(BW ~ Group, data.frame(BW = vcapwide[,"BWWeek14"], Group = vcapwide[,"Group"])))
}
\keyword{datasets}
