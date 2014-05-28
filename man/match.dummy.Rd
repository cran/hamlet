\name{match.dummy}
\alias{match.dummy}
\title{
Create dummy individuals or sinks to a data matrix or a distance/dissimilarity matrix
}
\description{
Dummy observations are allowed in order to make the number of observations dividable by the number of elements in each submatch, i.e. for pairwise matching the number of observations should be paired, for triangular matching the number of observations should be dividable by 3, etc. This can be done either by adding column averaged individuals to the original data frame (parameter 'dat'), or by adding zero distance sinks to the distance/dissimilarity matrix (parameter 'd'). The latter approach favors dummies being matched to real extreme observations, while the former favors dummies being matched to close-to-mean real observations.
}
\usage{
match.dummy(dat, d, g = 2)
}
\arguments{
  \item{dat}{
A data.frame of the original observations, to which column averaged new dummy observations are added
}
  \item{d}{
N times N distance/dissimilarity matrix, to which zero distance sinks are added
}
  \item{g}{
The desired number of elements per each submatch, i.e. the size of the clusters. The number of added dummies is the smallest number of additions that fulfills (N+dummy)%%g == 0
}
}
\value{
Depending on if the dat or the d parameter was provided, the function either:
dat: adds new averaged individuals according to column means and then returns the data matrix
d: adds zero distance sinks to the distance/dissimilarity matrix and returns the new distance/dissimilarity matrix
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\note{
Adding zero distance sinks to the distance matrix or averaged individuals to the original data frame produce different results and affect the optimal matching task differently.
}
\seealso{
\code{\link{match.allocate}}
\code{\link{match.mat2vec}}
\code{\link{match.vec2mat}}
\code{\link{match.bb}}
}
\examples{
data(vcapwide)

exdat <- vcapwide[1:10,c("PSAWeek10", "BWWeek10")]
dim(exdat)
avgdummies <- match.dummy(dat=exdat, g=3)
dim(avgdummies)
# Construct an Euclidean distance matrix after adding two dummy individuals 
# (averaged individuals to the original data matrix)
bb3 <- match.bb(as.matrix(dist(avgdummies)), g=3)
str(bb3)

# Construct an Euclidean distance matrix after adding two dummy distances (zero distance sinks)
exd <- as.matrix(dist(vcapwide[1:10,c("PSAWeek10", "BWWeek10")]))
dim(exd)
d <- match.dummy(d=exd, g=3)
dim(d)
# 10 is not dividable by 3, 2 sinks are added to make d 12x12
bb3 <- match.bb(d, g=3)
str(bb3)

# Notice that sinks produce a lot smaller target function costs than averaged individuals
}
\keyword{ manip }
