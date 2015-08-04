\name{hmap.annotate}
\alias{hmap.annotate}
\title{
Add a row and column annotations to a plot-region based heatmap built with 'hmap'
}
\description{
Annotation of rows or columns in a 'hmap'-plot. By default, rectangles aligned with either rows or columns are plotted to the right-side or lower-side of the heatmap respectively. User-specified customizations may be given to change these annotations in positioning or type.
}
\usage{
hmap.annotate(h, rw, rw.n = length(unique(rw)), rw.col = rainbow(rw.n,
                 start = 0.05, end = 0.5), rw.wid, rw.hei, rw.pch,
                 rw.x = rep(min(h$rightlim), times =
                 length(h$rowtext$xseq)), rw.y = h$rowtext$yseq, rw.shift
                 = c(0.02, 0), cl, cl.n
                 = length(unique(cl)), cl.col = rainbow(cl.n, start =
                 0.55, end = 1), cl.wid, cl.hei, cl.pch, cl.x =
                 h$coltext$xseq, cl.y = rep(max(h$bottomlim), times =
                 length(h$coltext$yseq)), cl.shift = c(0, -0.02), ...)
}
\arguments{
  \item{h}{
The list of heatmap parameters returned invisibly by the original 'hmap'-call.
}
  \item{rw}{
Annotation vector for rows 'r', each unique instance is given a different color (or pch) and plotted right-side of the corresponding heatmap rows
}
  \item{rw.n}{
Number of unique colors (or pch) to give each annotated row
}
  \item{rw.col}{
A vector for color values for unique instances in 'r' for annotating rows
}
  \item{rw.wid}{
The widths for annotation boxes for each row 'r'
}
  \item{rw.hei}{
The heights for annotation boxes for each row 'r'
}
  \item{rw.pch}{
Alternatively, instead of widths and heights user may specify a symbol 'pch' to use for annotating each row
}
  \item{rw.x}{
The x-coordinate locations for the row annotations, by default right side of heatmap itself
}
  \item{rw.y}{
The y-coordinate locations for the row annotations, by default same vertical locations as for the heatmap rows
}
  \item{rw.shift}{
Row annotation shift: a vector of 2 values, where first indicates the amount of x-axis shift desired and the second indicates the amount of y-axis shift
}
  \item{cl}{
Annotation vector for columns 'r', each unique instance is given a different color (or pch) and plotted lower-side of the corresponding heatmap columns
}
  \item{cl.n}{
Number of unique colors (or pch) to give each annotated column
}
  \item{cl.col}{
A vector for color values for unique instances in 'c' for annotating columns
}
  \item{cl.wid}{
The widths for annotation boxes for each column 'c'
}
  \item{cl.hei}{
The heights for annotation boxes for each column 'c'
}
  \item{cl.pch}{
Alternatively, instead of widths and heights user may specify a symbol 'pch' to use for annotating each column
}
  \item{cl.x}{
The x-coordinate locations for the column annotations, by default same horizontal locations as for the heatmap columns
}
  \item{cl.y}{
The y-coordinate locations for the column annotations, by default lower side of heatmap itself
}
  \item{cl.shift}{
Column annotation shift: a vector of 2 values, where first indicates the amount of x-axis shift desired and the second indicates the amount of y-axis shift
}
  \item{...}{
Additional parameters supplied either to 'rect' or 'points' function if user desired rectangles or 'pch'-based points respectively
}
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{heatmap}}
\code{\link{hmap.key}}
\code{\link{hmap}}
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
