\name{mem.powersimu}
\alias{mem.powersimu}
\title{
Power simulations for the fixed effects of a mixed-effects model through 
structured bootstrapping of the data and re-fitting of the model
}
\description{

Bootstrap sampling is used to investigate the statistical significance of the fixed 
effects terms specified for a readily fitted mixed-effects model as a function of the 
number of individuals participating in the study. User either specifies a
suitable sampling unit, or it is automatically identified based on the random effects
formulation of a readily fitted mixed-effects model. Per each count of individuals in 
vector N, a fixed number of bootstrapped datasets are generated and re-fitted using
the model formulation on the pre-fitted model. Power is then computed as the fraction
of effects identified as statistically significant out of all the bootstrapped datasets.

}
\usage{
mem.powersimu(fit, N = 4:20, boot = 100, level = NULL, strata = NULL, 
	default = FALSE, seed = NULL, plot = TRUE, plot.loess = FALSE, 
	legendpos = "bottomright", return.data = FALSE, verb = 1, ...)
}
\arguments{
  \item{fit}{
A fitted mixed-effects model. Should be either a model produced by the 
lme4-package, or then a modified lme4-fit such as provided by lmerTest 
or similar package that builds on lme4.
}
  \item{N}{
A vector of desired amounts of individuals to be tested, i.e. sample sizes N.
Notice that the N may be either a total N if no strata is spesified, or then an 
N value per each substrata if strata is not NULL. See below the parameter 'strata'.
}
  \item{boot}{
Number of bootstrapped datasets to generate per each N value. The total number of 
generated data frames in the end will be N times boot.
}
  \item{level}{
An unambiguous indicator available in the model data frame that indicates each separate
individual unit in the experiment. For example, this may correspond to a single patient
indicator column ID, where each patient has a unique ID instance. If this parameter is 
given as NULL, then this function automatically attempts to identify the best possible
level of individual indicators based on the random effects specified for the model.
}
  \item{strata}{
If any sampling strata should be balanced, it should be indicated here. For example, if
one is studying the possible effects of an intervention, it is typical to have an equal 
number of individual both in the control and in the intervention arms also in the sampled
datasets. It should be then given as an column name available in the original model data
frame. Each strata will be sampled in equal amounts.
}
  \item{default}{
What is the default statistical significance if a model could not be re-fitted to the 
sampled datasets, which may occur for example due to convergence or redundance issues.
This defaults to FALSE, which means that a coefficient is expected to be statistically 
insignificant if the corresponding model re-fitting fails in lme4.
}
  \item{seed}{
For reproducibility, one may wish to set a numeric seed to produce the exact same results.
}
  \item{plot}{
If set to TRUE, the function will plot a power curve. Each fixed effects coefficient is a 
different curve, with color coding and a legend annotated to separate which one is which.
}
  \item{plot.loess}{
If plot==TRUE, this plot.loess==TRUE adds an additional loess-smoothed approximated curve 
to the existing curves. This is useful if running the simulations with a low number of 
bootstrapped samples, as it may help approximate where the curve reaches critical points, i.e.
power = 0.8.
}
  \item{legendpos}{
Position for the legend in plot==TRUE, defaults to "bottomright". Any legal position similar 
to provided the function 'legend' is allowed.
}
  \item{return.data}{
Should one obtain the bootstrapped data instead of bootstrapping and then re-fitting. This will 
skip the model re-fitting schema and instead return a list of lists with the bootstrapped data 
instead. The outer list corresponds to the values of 'N', while the inner loop corresponds to the
different 'boot' runs of bootstrap. This may be useful to inspecting that the schema is sampling
correct sampling units for example, or if bootstrapping is to be used for something else than 
re-fitting the lme4-models.
}
  \item{verb}{
Numeric value indicating the level of verbosity; 0=silent, 1=normal, 2=debugging.
}
  \item{\dots}{
Additional parameters provided for the function.
}
}
\details{
This function will by default utilizes the lmerTest-package's Satterthwaite approximation for determining the 
p-values for the fixed effects. If this fails, it resorts to the conventional approximation |t|>2 for significance,
which is not accurate, but may provide a reasonable approximation for the power levels.
}
\value{
If return.data==FALSE, this function will return a matrix, where the rows correspond to the different
N values and the columns correspond to the fixed effects. The values [0,1] are the fraction of bootstrapped
datasets where the corresponding fixed effects was detected as statistically significant.
}
\author{
Teemu D. Laajala
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mem.getcomp}}
}
\examples{
# Use the VCaP ARN data as an example
data(vcaplong)
arn <- vcaplong[vcaplong[,"Group"] == "Vehicle" | vcaplong[,"Group"] == "ARN",]

# lme4 is required for mixed-effects models
library(lme4)
# Fit an example fixed effects model
fit <- lmer(PSA ~ 1 + DrugWeek + ARN:DrugWeek + (1|ID) + (0 + DrugWeek|ID), data = arn)

# For reproducibility, set a seed
set.seed(123)
# Run a brief power analysis with only a few selected N values and a limited number of bootstrapping
# Balance strata over the ARN and non-ARN (=Vehicle) so that both contain equal count of individuals
power <- mem.powersimu(fit, N=c(4, 7, 10), boot=20, strata="ARN", plot=TRUE)
# Power curves are plotted, along with returning the power matrix at:
power

# Notice that each column corresponds to a specified fixed effects at the formula 
# "1 + DrugWeek + ARN:DrugWeek"

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ power }
\keyword{ design }
\keyword{ mem }
