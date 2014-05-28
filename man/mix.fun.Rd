\name{mix.fun}
\alias{mix.fun}
\title{
Apply function to numerical columns of a mixed data.frame while ignoring non-numeric fields
}
\description{
This function is intended for applying functions to numeric fields of a mixed type data.frame. Namely, the function ignores fields that are e.g. factors, and returns FUN function applied to only the numeric fields.
}
\usage{
mix.fun(x, FUN = scale, ...)
}
\arguments{
  \item{x}{
Data.frame x with mixed type fields
}
  \item{FUN}{
Function to apply, for example 'scale', 'cov', or 'cor'
}
  \item{\dots}{
Additional parameters passed on to FUN
}
}
\value{
Return values of FUN when applied to numeric columns of 'x'
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{apply}}
}
\examples{
data(vcapwide)

mix.fun(vcapwide[,c("Group", "PSAWeek4", "PSAWeek10", "PSAWeek14")], FUN=scale) 
# Column 'Group' is ignored
mix.fun(vcapwide[,c("Group", "PSAWeek4", "PSAWeek10", "PSAWeek14")], FUN=cov, use="na.or.complete") 
# ... is used to pass the 'use' parameter to the 'cov'-function
}
\keyword{ manip }
