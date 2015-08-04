\name{hmap.key}
\alias{hmap.key}
\title{
Add a color key to a plot-region based heatmap built with 'hmap'
}
\description{
A continuous color scale key for a heatmap. By default the key is constructed according to the 'h'-object which is invisibly returned by the original 'hmap'-call. Some customization may be supplied to position the legend or to customize ticks and style.
}
\usage{
hmap.key(h, x0 = h$leftlim[1], x1 = h$leftlim[2], y0 =
                 h$toplim[1], y1 = h$toplim[2], xlim = range(h$valseq),
                 ratio = 0.5, tick = 0.1, at = seq(from =
                 min(h$valseq), to = max(h$valseq), length.out = 5),
                 bty = "c", cex = 0.5, pos = 3)
}
\arguments{
  \item{h}{
The list of heatmap parameters returned invisibly by the original 'hmap'-call.
}
  \item{x0}{
Coordinates for the color key; left border
}
  \item{x1}{
Coordinates for the color key; right border
}
  \item{y0}{
Coordinates for the color key; lower border
}
  \item{y1}{
Coordinates for the color key; upper border
}
  \item{xlim}{
Value range for the x-axis within the key itself, by default extracted from the h-object
}
  \item{ratio}{
Ratio between y-axis coordinates to separate the key box to upper color key box and lower tick and values
}
  \item{tick}{
The vertical length in value ticks
}
  \item{at}{
The values in color key at which to plot ticks and the values at ticks
}
  \item{bty}{
Type of box to plot around the color key
}
  \item{cex}{
The zooming factor for plotting the text and other objects affected by the 'cex' parameter in 'par'
}
  \item{pos}{
The text alignment and position parameter given to the 'text' function in the key
}
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{heatmap}}
\code{\link{hmap}}
\code{\link{hmap.annotate}}
}
\examples{
# Generate some data
set.seed(1)
r1 <- replicate(30, rnorm(20))
lab <- sample(letters[1:2], 20, replace=TRUE)
r1[lab==lab[1],] <- r1[lab==lab[1],] + 2
r2a <- replicate(10, rnorm(10))
r2b <- replicate(10, rnorm(10))

# Set up a new plot region, notice we have a 2-fold wider x-axis
plot.new()
plot.window(xlim=c(0,2), ylim=c(0,1))

# Plot an example plot along with the color key and annotations for our 'lab' vector
h1 <- hmap(r1, add = TRUE)
hmap.key(h1, x1=0.18)
hmap.annotate(h1, rw = lab, rw.wid=c(0.82,0.90))

# Plot the rest to show how the coordinates are adjusted to place the heatmap(s) differently
h2a <- hmap(r2a, add = TRUE, xlim=c(1.2, 1.8), leftlim=c(1.0, 1.2), 
	rightlim=c(1.8,2.0), ylim=c(0.6, 1.0), bottomlim=c(0.5,0.6), Colv=NA)
h2b <- hmap(r2b, add = TRUE, xlim=c(1.2, 1.8), leftlim=c(1.0, 1.2), 
	rightlim=c(1.8,2.0), ylim=c(0.1, 0.5), bottomlim=c(0.0,0.1), Colv=NA)

# Show the normal plot region axes 
axis(1, at=c(0.5,1.5), c("A", "B"))

}
\keyword{ hplot }
\keyword{ aplot }
