\name{hamlet-package}
\alias{hamlet-package}
\alias{hamlet}
\docType{package}
\title{
Hierarchical Optimal Matching and Machine Learning Toolbox
}
\description{
This package provides functions and algorithms for solving optimal matching tasks in the context of preclinical cancer studies. Further, various help and plotting functions are provided for unsupervised and supervised machine learning as well as longitudinal modeling of tumor growth response patterns.
}
\details{
\tabular{ll}{
Package: \tab hamlet\cr
Type: \tab Package\cr
Version: \tab 0.9.5\cr
Date: \tab 2016-03-01\cr
License: \tab GPL (>= 2)\cr
}
The package 'hamlet' offers functions for optimal matching, randomization, and mixed-effects modeling in preclinical cancer studies. The functions are divided to 'match'-prefix indicating optimal matching intended functions, 'mem' indicating mixed-effects modeling, 'mix' for mixed type (numerical and categorical) data analysis, and rest that are various plotting and helper functions for various tasks.
}
\author{
Teemu Daniel Laajala

Maintainer: Teemu Daniel Laajala <teelaa@utu.fi>
}
\references{
Laajala TD, Jumppanen M, Huhtaniemi R, Fey V, Kaur A, et al. (2016) Optimized design and analysis of preclinical intervention studies in vivo. Sci Rep. 2016 Aug 2;6:30723. doi: 10.1038/srep30723.

Knuuttila M, Yatkin E, Kallio J, Savolainen S, Laajala TD, et al. (2014) Castration induces upregulation of intratumoral androgen biosynthesis and androgen receptor expression in orthotopic VCaP human prostate cancer xenograft model.  Am J Pathol. 2014 Aug;184(8):2163-73. doi: 10.1016/j.ajpath.2014.04.010.
}
\keyword{ package }
\examples{

##
## Exploring the VCaP dataset provided alongside the 'hamlet' package
##

data(vcapwide)
data(vcaplong)

# VCaP Castration-resistant prostate cancer (CRPC) PSA-measurements (and body weight) in wide-format
mixplot(vcapwide[,c("PSAWeek10", "PSAWeek14", "BWWeek10", "Group")], pch=16)
anv <- aov(PSA ~ Group, data.frame(PSA = vcapwide[,"PSAWeek14"], Group = vcapwide[,"Group"]))
summary(anv)
TukeyHSD(anv)
summary(aov(BW ~ Group, data.frame(BW = vcapwide[,"BWWeek14"], Group = vcapwide[,"Group"])))

# VCaP Castration-resistant prostate cancer (CRPC) PSA-measurements (and body weight) in long-format
library(lattice)
xyplot(log2PSA ~ DrugWeek | Group, data = vcaplong, type="l", group=ID, layout=c(3,1))
xyplot(BW ~ DrugWeek | Group, data = vcaplong, type="l", group=ID, layout=c(3,1))

##
## Example multigroup (g=3) nbp-matching using the branch and bound algorithm, 
## and subsequent random allocation of submatches to 3 arms
##

# Construct an Euclidean distance example distance matrix using 15 observations from the VCaP study
d <- as.matrix(dist(vcapwide[1:15,c("PSAWeek10", "BWWeek10")]))
# Matching using the b&b algorithm to submatches of size 3 
# (which will result in 3 intervention groups)
bb3 <- match.bb(d, g=3)
str(bb3)

solvec <- bb3$solution 
# matching vector, where each element indicates to which submatch each observation belongs to

# Perform an example random allocation of the above submatches, 
# these will be randomly allocated to 3 arms based on the submatches
set.seed(1)
groups <- match.allocate(solvec)

# Illustrate randomization, no baseline differences in these three artificial groups
by(vcapwide[1:15,c("PSAWeek10", "BWWeek10")], INDICES=groups, FUN=function(x) x)

summary(aov(PSAWeek10 ~ groups, data = data.frame(PSAWeek10 = vcapwide[1:15,"PSAWeek10"], groups)))
summary(aov(BWWeek10 ~ groups, data = data.frame(BWWeek10 = vcapwide[1:15,"BWWeek10"], groups)))

##
## Example mixed-effects modeling of the longitudinal PSA profiles using 
## the actual experimental groups
##

exdat <- vcaplong[vcaplong[,"Group"] \%in\% c("Vehicle", "ARN"),]

library(lme4)
# Model fitting using lme4-package
f1 <- lmer(log2PSA ~ 1 + DrugWeek + DrugWeek:ARN + (1 + DrugWeek|ID), data = exdat)

mem.getcomp(f1)

library(lmerTest)
# Model term testing using the lmerTest-package
summary(f1) 

}
