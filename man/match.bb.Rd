\name{match.bb}
\alias{match.bb}
\title{
Branch and Bound algorithm implementation for performing multigroup non-bipartite matching
}
\description{
This function performs multigroup non-bipartite matching of observations based on a provided distance/dissimilarity matrix 'd'. The number of elements in each submatch is defined by the parameter 'g'.
}
\usage{
match.bb(d, g = 2, presort = "complete", progress = 1e+05, 
bestknown = Inf, maxbranches = Inf, verb = 0)
}
\arguments{
  \item{d}{
A distance matrix with NxN elements
}
  \item{g}{
Number of elements per each submatch, i.e. how many observations are always matched together
}
  \item{presort}{
If hierarchical clustering should be used for an initial guess, hclust method-options are valid options ("complete", "single", "ward", "average")
}
  \item{progress}{
How many branching operations are done before outputting information to the user
}
  \item{bestknown}{
If a best known solution already exists, this can be used to bound branches of the tree before initiation. The default Inf value causes whole search tree to be potential solution space.
}
  \item{maxbranches}{
Maximum number of branching operations before returning current best solution, by default no cutoff is defined.
}
  \item{verb}{
Level of verbosity
}
}
\details{
See further details in the reference Laajala et al.
}
\value{
The function returns a list of objects, where elements are
\item{branches }{Number of branching operations during the branch and bound algorithm}
\item{bounds }{Number of bounding operations during the branch and bound algorithm}
\item{ends }{Number of end leaf nodes visited during the branch and bound algorithm}
\item{matrix }{The resulting binary matching matrix where rows and columns sum to g}
\item{solution }{The resulting matching vector where each element indicates the submatch where the observation was placed}
\item{cost }{Final cost value of the target function in the minimization task}
}
\references{
Laajala TD, Kaur A, Knuuttila M, Westermarck J, et al. (2014) Optimal Matching, Randomization and Statistical Modeling of Hierarchical Baseline Variables in Preclinical Cancer Studies. Submitted
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\note{
Notice that the solution submatch vector in $solution is not the same as the intervention group allocation. Submatches should be randomly allocated to intervention arms using the match.allocate-function.

The package 'nbpMatching' provides a FORTRAN implementation for computation of paired non-bipartite matching case (g=2).

Computation may be heavy if the number of observations is high, or the number of within-submatch pairwise distances to consider is high (increases quadratically as a function of 'g').
}
\seealso{
\code{\link{match.allocate}}
\code{\link{match.mat2vec}}
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
# binary matching matrix
solvec <- bb3$solution 
# matching vector, where each element indicates to which submatch each observation belongs to

mixplot(data.frame(vcapwide[1:15,c("PSAWeek10", "BWWeek10")], 
 submatch=as.factor(paste("Submatch_",solvec, sep=""))), pch=16, col=rainbow(5))
}
\keyword{ cluster }
\keyword{ design }
