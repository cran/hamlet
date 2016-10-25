\name{orxlong}
\alias{orxlong}
\docType{data}
\title{
Long-format longitudinal data for the ORX study
}
\description{
Long-format measurements of PSA over the intervention period in the ORX study. Notice that this data.frame
is in suitable format for mixed-effects modeling, where each row should correspond to a single
longitudinal measurement. These measurements are annotated using the individual indicator fields 'ID',
time fields 'Day', 'TrDay', 'Date', and the response values are contained in raw format in 'PSA' or
after log2-transformation in 'log2PSA'. Additional fields are provided for group testing and matched 
inference in 'Group', 'Submatch', and the binary indicators 'ORX+Tx', 'ORX', and 'Intact'.
}
\usage{data("orxlong")}
\format{
  A data frame with 392 observations on the following 11 variables.
  \describe{
    \item{\code{ID}}{A unique character indicator for the different individual(s)}
    \item{\code{PSA}}{Raw longitudinal PSA measurement values in unit (ug/l)}
    \item{\code{log2PSA}}{Log2-transformed longitudinal PSA measurement values in unit (log2 ug/l)}
    \item{\code{Day}}{Day since the first PSA measurement. Notice that there is a single time point prior to interventions.}
    \item{\code{TrDay}}{Day since the interventions began, 0 annotating the point at which surgery was performed or drug compounds were first given.}
    \item{\code{Date}}{A date format when the actual measurement was performed}
    \item{\code{Group}}{The actual intervention groups, after blinded groups were assigned to 'ORX+Tx', 'ORX', or 'Intact'}
    \item{\code{Submatch}}{The submatches that were assigned based on the baseline variables.}
    \item{\samp{ORXTx}}{A binary indicator field indicating which measurements belong to the group 'ORX+Tx'}
    \item{\code{ORX}}{A binary indicator field indicating which measurements belong to the group 'ORX'}
    \item{\code{Intact}}{A binary indicator field indicating which measurements belong to the group 'Intact'}
  }
}
\note{
Group-testing should be performed so that 'ORX+Tx' is tested against 'ORX', in order to infer possible effects 
occurring due to 'Tx' on top of 'ORX'. 'ORX' should be compared to 'Intact', in order to infer if the 'ORX' surgical
procedure has beneficial effects in comparison to intact animals. For statistical modeling of the intervention effects,
one should use observations with the positive 'TrDay'-values, as this indicates the beginning of the interventions.
}
\details{
For mixed-effects modeling, the fields 'ID', 'PSA' (or 'log2PSA'), 'TrDay', and group-specific indicators should be included.
}
\source{
Laajala TD, Jumppanen M, Huhtaniemi R, Fey V, Kaur A, et al. (2016) Optimized design and analysis of preclinical intervention studies in vivo. Sci Rep. 2016 Aug 2;6:30723. doi: 10.1038/srep30723.
}
\examples{
data(orxlong)
# Construct data frames that can be used for testing pairwise group contrasts
orxintact <- orxlong[orxlong[,"Intact"]==1 | orxlong[,"ORX"]==1, 
	c("PSA", "ID", "ORX", "TrDay", "Submatch")]
orxtx <- orxlong[orxlong[,"ORXTx"]==1 | orxlong[,"ORX"]==1, 
	c("PSA", "ID", "ORXTx", "TrDay", "Submatch")]
# Include only observations occurring post-surgery
orxintact <- orxintact[orxintact[,"TrDay"]>=0,]
orxtx <- orxtx[orxtx[,"TrDay"]>=0,]

# Example fits
library(lme4)
library(lmerTest)
# Conventional model fits
fit1a <- lmer(PSA ~ 1 + TrDay + ORXTx:TrDay + (1|ID) + (0 + TrDay|ID), data = orxtx)
fit1b <- lmer(PSA ~ 1 + TrDay + ORXTx:TrDay + (1 + TrDay|ID), data = orxtx)
fit2a <- lmer(PSA ~ 1 + TrDay + ORX:TrDay + (1|ID) + (0 + TrDay|ID), data = orxintact)
fit2b <- lmer(PSA ~ 1 + TrDay + ORX:TrDay + (1 + TrDay|ID), data = orxintact)

# Collate to matched inference for pairwise observations over the submatches
matched.orx <- do.call("rbind", by(orxintact, INDICES=orxintact[,"Submatch"], FUN=function(z){
	z[,"MatchedPSA"] <- z[,"PSA"] - z[z[,"ORX"]==0,"PSA"]
	z <- z[z[,"ORX"]==1,]
	z
}))
# Few examples of matched fits with different model formulations
fit.matched.1 <- lmer(MatchedPSA ~ 0 + TrDay + (1|ID) + (0 + TrDay|ID), data = matched.orx)
fit.matched.2 <- lmer(MatchedPSA ~ 1 + TrDay + (1|ID) + (0 + TrDay|ID), data = matched.orx)
fit.matched.3 <- lmer(MatchedPSA ~ 1 + TrDay + (1 + TrDay|ID), data = matched.orx)
summary(fit.matched.1)
summary(fit.matched.2)
summary(fit.matched.3)
# We notice that the intercept term is highly insignificant 
# if included in the matched model, as expected by baseline balance.
# In contrast, the matched intervention growth coefficient is highly
# statistically significant in each of the models.
}
\keyword{datasets}
