\name{mixplot}
\alias{mixplot}
\title{
Scatterplot for mixed type data
}
\description{
This function plots a scatterplot similar to the default plot-function, with the difference that factor/character fields in input data.frame are handled as categorical variables. These categorical variables are color-coded and handled separately in marginal distributions.
}
\usage{
mixplot(x, 
	main = NA, 
	match, 
	func = function(x, y, par) 
		{ segments(x0 = x[1], y0 = x[2], x1 = y[1], y1 = y[2], col = par)}, 
	legend = T, 
	col = palette(), na.lines = T, 
	origin = F, 
	marginal = F, 
	lhei, 
	lwid, 
	verb = 0,
...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data.frame or a matrix of observations. Typically x should be a data.frame, where columns are of different types, e.g. some of 'numeric' and some of 'factor' class.
}
  \item{main}{
Main title plotted on top of the figure
}
  \item{match}{
A matching matrix (e.g. produced by hamlet::match.vec2mat) or a matching vector (e.g. produced by hamlet::match.mat2vec) that indicates with different values if certain observations should be connected.
}
  \item{func}{
The function to apply to each pair of observations 'x' and 'y'. By default, it is a segment line in 2 dimensions (each individual bivariate panel). Segment line color is indicated by the matching vector or individual element in the matching matrix. Thus 0-values indicate no line, while other values are used to annotate submatches. 'par' is the index of the submatch, and by default indicate the colors.
}
  \item{legend}{
Should an automated legend be generated
}
  \item{col}{
Colors per observation
}
  \item{na.lines}{
Should lines be drawn to represent one of the variables if the other one is missing in a 2-dim scatterplot
}
  \item{origin}{
Should the origin x=0, y=0 be separately indicated using lines
}
  \item{marginal}{
Should marginal distributions be drawn in sides of each scatterplot
}
  \item{lhei}{
Heights for bins in the layout
}
  \item{lwid}{
Widths for bins in the layout
}
  \item{verb}{
Level of verbosity: -1<= (no verbosity), 0/FALSE (warnings) or >=1/TRUE (additional information)
}
  \item{\dots}{
Additional parameters given to the plot-function
}
}
\value{
An invisible return of the measurements and plot layout structure (matrix, heights, and widths)
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\examples{
data(vcapwide)

mixplot(vcapwide[,c("Group", "PSAWeek4", "PSAWeek10", "PSAWeek14")], marginal=TRUE, pch=16, 
 main="PSA at weeks 4, 10 and 14 per intervention group")
}
\keyword{ hplot }
