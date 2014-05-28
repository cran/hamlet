\name{extendsymrange}
\alias{extendsymrange}
\title{
Extend range of variable limits while retaining a point of symmetricity
}
\description{
This function serves as an alternative to the R function 'extendrange', when user wishes to conserve a point of symmetricity for the range. For example, this might be desired when the plot should be symmetric around the origin x=0, but that the sides need to extend beyond the actual range of values.
}
\usage{
extendsymrange(x, r = range(x, na.rm = T), f = 0.05, sym = 0)
}
\arguments{
  \item{x}{
Vector of values to compute the range for
}
  \item{r}{
The range of values
}
  \item{f}{
The factor by which the range is extended beyond the extremes
}
  \item{sym}{
The defined point of symmetricity
}
}
\value{
A vector of 2 values for the lower and higher limit of the symmetric extended range
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{extendrange}}
}
\examples{
set.seed(1)
ex <- rnorm(10)+2

hist(ex, xlim=extendsymrange(ex, sym=0), breaks=100)
}
\keyword{ dplot }
