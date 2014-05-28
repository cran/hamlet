\name{smartjitter}
\alias{smartjitter}
\title{
Smart jittering function for deterministic shifting of overlapping values
}
\description{
This function takes in a vector of measurements and computes overlapping bins of observations, and applies a jittering function within each overlapping bin.
}
\usage{
smartjitter(x, q = seq(from = 0, to = 1, length.out = 10), type = 1, 
amount = 0.1, jitterfuncs = list(function(n) {
    (1:n)/(1/amount)
}, function(n) {
    (((-1)^c(0:(n - 1))) * (0:(n - 1)))/(1/amount)
}), jits = jitterfuncs[[type]])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The values that should be jittered. Notice that these are used to determine which are overlapping, and should not be though of as x-axis positions (see example).
}
  \item{q}{
Probability quantiles where the ends of the bins should be placed
}
  \item{type}{
Type of jittering, by default it is used to choose which element (1 or 2) of the list of jittering functions is chosen as the final jittering function. Customized functions may be provided to the jitterfuncs-parameter.
}
  \item{amount}{
Amount of jittering (here deterministic shifting) for the jittering function
}
  \item{jitterfuncs}{
List of possible jittering functions for n overlapping values. The jittering function at list position 'type' is chosen
}
  \item{jits}{
Final jittering function from the jitterfuncs-list
}
}
\details{
The smart jittering is applied to the x-parameter values, and returns a vector of shifting amounts per each observation. Notice that in the typical case, parameter 'x' are the desired response values e.g. among the y-axis, and the returned value of smartjitter are the amounts of jittering done on the x-axis of a plot.
}
\value{
The function returns a vector of values with same length as x. The values in this vector indicate what should be the shifting per each observation, if the observations should be jittered along an another axis.
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\examples{
data(vcapwide)

plot.new()
plot.window(xlim=extendrange(c(0,1)), ylim=extendrange(vcapwide[,"PSAWeek4"]))
y1 <- vcapwide[vcapwide[,"CastrationDate"]=="100413","PSAWeek4"]
y2 <- vcapwide[vcapwide[,"CastrationDate"]=="170413","PSAWeek4"]
points(x=0+smartjitter(y1, type=2, amount=0.02), y=y1)
points(x=1+smartjitter(y2, type=2, amount=0.02), y=y2)
axis(1, at=c(0,1), labels=c("10.04.13", "17.04.13"))
axis(2); box()
title(ylab="PSA at week 4", xlab="Castration batches")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
