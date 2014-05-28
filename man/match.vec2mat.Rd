\name{match.vec2mat}
\alias{match.vec2mat}
\title{
Transform a matching vector to a binary matching matrix
}
\description{
This function allows transforming a matching vector to a binary matching matrix. A matching vector is of length N where each element indicates the submatch to which the observation belongs to. Notice that this is not the same as the group allocation vector that is provided by the match.allocate-function. The binary matching matrix is of size N x N where 0 indicates that the observations have been part of a different submatch, and 1 indicates that the observations have been part of the same submatch. Diagonal is always 0 although an observation is always in the same submatch with its self.
}
\usage{
match.vec2mat(x)
}
\arguments{
  \item{x}{
A matching vector 'x'
}
}
\value{
N times N binary matching matrix, where 0 indicates that the observations have been part of a different submatch, and 1 indicates that the observations have been part of the same submatch.}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{match.allocate}}
\code{\link{match.mat2vec}}
\code{\link{match.bb}}
\code{\link{match.dummy}}
}
\examples{
data(vcapwide)

# Construct an Euclidean distance example distance matrix using 15 observations from the VCaP study
d <- as.matrix(dist(vcapwide[1:15,c("PSAWeek10", "BWWeek10")]))

bb3 <- match.bb(d, g=3)
str(bb3)

solvec <- bb3$solution 
# matching vector, where each element indicates to which submatch each observation belongs to

solvec
mat <- match.vec2mat(solvec)
mat
which(mat[1,] == 1) 
# E.g. the first, third and thirteenth observation are part of the same submatch
which(solvec == solvec[1]) 
# Similarly
}
\keyword{ manip }
