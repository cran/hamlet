\name{hmap}
\alias{hmap}
\title{
Plot-region based heatmap
}
\description{
This function plots heatmap figure based on the normal plot-region. This is useful if the image-based function 'heatmap' is not suitable, i.e. when multiple heatmaps should be placed in a single device.
}
\usage{
hmap(x, add = F, xlim = c(0, 1), ylim = c(0, 1), col = heat.colors(10), 
border = matrix(NA, nrow = nrow(x), ncol = ncol(x)), 
lty = matrix("solid", nrow = nrow(x), ncol = ncol(x)), 
lwd = matrix(1, nrow = nrow(x), ncol = ncol(x)), ...)
}
\arguments{
  \item{x}{
Matrix to be plotted
}
  \item{add}{
Should the figure be added to the plotting region of an already existing figure
}
  \item{xlim}{
The x limits in which the heatmap is placed horizontally in the plotting region
}
  \item{ylim}{
The y limits in which the heatmap is placed vertically in the plotting region
}
  \item{col}{
Color palette for the heatmap colors
}
  \item{border}{
A matrix of border color definitions (rectangles in the heatmap)
}
  \item{lty}{
A matrix of line type definitions (rectangles in the heatmap)
}
  \item{lwd}{
A matrix of line width definitions (rectangles in the heatmap)
}
  \item{\dots}{
Additional parameters provided to the rectangle plotting function
}
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\note{
The function does not yet implement the hclust-functionality of the normal R heatmap-function.
}
\seealso{
\code{\link{heatmap}}
}
\examples{
set.seed(1)
ex <- matrix(rnorm(100), nrow=10, ncol=10)
highlight <- matrix("black", nrow=10, ncol=10)
highlight[ex>1] <- "red"
thick <- matrix(1, nrow=10, ncol=10)
thick[ex>1] <- 3

plot.new()
plot.window(xlim=c(0,3), ylim=c(0,1))
hmap(ex, add=TRUE, xlim=c(0,1))
hmap(ex, add=TRUE, xlim=c(2,3), border=highlight, lwd=thick)
}
\keyword{ hplot }
\keyword{ aplot }
