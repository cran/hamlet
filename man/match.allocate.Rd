\name{match.allocate}
\alias{match.allocate}
\title{
Allocation of matched units to intervention arms
}
\description{
This function allocates units belonging to a single submatch to separate intervention arms. This ensures that the resulting intervention groups are homogeneous in respect to the variables that were used to construct the distance/dissimilarity matrix for the non-bipartite matching. The number of resulting intervention groups is equal to the 'g' (i.e. submatch size) used in the multigroup non-bipartite matching.
}
\usage{
match.allocate(xmat)
}
\arguments{
  \item{xmat}{
A binary matching matrix or a matching vector given by match.bb-function.
}
}
\value{
A vector where each element indicates to which group the observation was randomized to. The group names are "Group_A", "Group_B", "Group_C", ... until 'g' letters, where 'g' was the size of submatches.
}
\references{
Laajala TD, Kaur A, Knuuttila M, Westermarck J, et al. (2014) Optimal Matching, Randomization and Statistical Modeling of Hierarchical Baseline Variables in Preclinical Cancer Studies. Submitted
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{match.bb}}
\code{\link{match.mat2vec}}
\code{\link{match.vec2mat}}
\code{\link{match.dummy}}
}
\examples{
data(vcapwide)

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
}
\keyword{ design }
