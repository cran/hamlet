\name{match.mat2vec}
\alias{match.mat2vec}
\title{
Transform a binary matching matrix to a matching vector
}
\description{
This function transforms a binary matching matrix to a matching vector. A matching vector is of length N where each element indicates the submatch to which the observation belongs to. Notice that this is not the same as the group allocation vector that is provided by the match.allocate-function. The binary matching matrix is of size N x N where 0 indicates that the observations have been part of a different submatch, and 1 indicates that the observations have been part of the same submatch. Diagonal is always 0 although an observation is always in the same submatch with its self.
}
\usage{
match.mat2vec(xmat)
}
\arguments{
  \item{xmat}{
A binary matching matrix 'xmat'
}
}
\value{
A matching vector where each element indicates submatch the observation belongs to
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\note{
Notice that the particular index numbers produced by match.mat2vec may be different to that of the branch and bound solution vector, but that the submatches shared by observations are common.
}
\seealso{
\code{\link{match.allocate}}
\code{\link{match.bb}}
\code{\link{match.vec2mat}}
\code{\link{match.dummy}}
}
\examples{
data(vcapwide)

# Construct an Euclidean distance example distance matrix using 15 observations from the VCaP study
d <- as.matrix(dist(vcapwide[1:15,c("PSAWeek10", "BWWeek10")]))

bb3 <- match.bb(d, g=3)
str(bb3)

mat <- bb3$matrix 
# matching vector, where each element indicates to which submatch each observation belongs to

mat
solvec <- match.mat2vec(mat)
which(mat[1,] == 1) 
# E.g. the first, third and thirteenth observation are part of the same submatch
which(solvec == solvec[1]) 
# Similarly
}
\keyword{ manip }
