\name{hmap}
\alias{hmap}
\title{
Plot-region based heatmap
}
\description{
This function plots heatmap figure based on the normal plot-region. This is useful if the image-based function 'heatmap' is not suitable, i.e. when multiple heatmaps should be placed in a single device.
}
\usage{
hmap(x, add = F, 
	xlim = c(0.2, 0.8), 
	ylim = c(0.2, 0.8),
	col = heat.colors(10), 
	border = matrix(NA, nrow = nrow(x), ncol = ncol(x)), 
	lty = matrix("solid", nrow = nrow(x), ncol = ncol(x)), 
	lwd = matrix(1, nrow = nrow(x), ncol = ncol(x)), 
	hclustfun = hclust, 
	distfun = dist, 
	reorderfun = function(d, w) reorder(d, w), 
	textfun = function(xseq, yseq, labels, type = "row", ...) 
		{ if (type == "col") par(srt = 90); 
		text(x = xseq, y = yseq, labels = labels, ...);  
		if (type == "col") par(srt = 0)}, 
	symm = F, 
	Rowv = NULL, 
	Colv = if (symm) Rowv else NULL, 
	leftlim = c(0, 0.2), toplim = c(0.8, 1), 
	rightlim = c(0.8, 1), bottomlim = c(0, 0.2), 
	type = "rect", 
	scale = c("none", "row", "column"), 
	na.rm = T, 
	nbins = length(col), 
	valseq = 
		seq(from = min(x, na.rm = na.rm), 
		to = max(x, na.rm = na.rm), length.out = nbins), 
	namerows = T, 
	namecols = T, 
	...)
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
  \item{hclustfun}{
The hierarchical clustering function similar to 'stats::heatmap' implementation. Should yield a valid 'hclust' object for a given distance/dissimilarity matrix.
}
  \item{distfun}{
The distance/dissimilarity function similar to 'stats::heatmap' implementation. Should yield a valid 'dist' object for a given data matrix.
}
  \item{reorderfun}{
The function to use to reorder branches of the clustering (notice that same-level branches in a hierarchical clustering may be permutated without violating the solution). The default approach from 'stats::heatmap' is utilized here.
}
  \item{textfun}{
A text function that is used to plot the names of the rows and columns, if desired. The default implementation shows how user could tailor the columns and rows differently, by turning the column labels around 90-degrees. The parameter 'type' is used to distinguish between rows and columns. 
}
  \item{symm}{
Should the given data matrix be treated as symmetric (has to be a square matrix if so), by default 'FALSE'.  
}
  \item{Rowv}{
The row clustering parameter. If 'NA' the row hierarchical clustering is completely omitted. Alternatively, if a numeric vector of ranks, the ordering of branches is tried to be permutated according to the desired order. This can also be a pre-computed dendrogram-object.
}
  \item{Colv}{
The column clustering parameter. If 'NA' the column hierarchical clustering is completely omitted. Alternatively, if a numeric vector of ranks, the ordering of branches is tried to be permutated according to the desired order. This can also be a pre-computed dendrogram-object.
}
  \item{leftlim}{
The horizontal limits of the row hierarchical clustering. The horizontal limits of the heatmap are {a=leftlim[1], b=leftlim[2], c=xlim[1], d=xlim[2]} where the 'a' is where the row dendrogram begins, 'b' is where the row dendrogram ends, 'c' is where the heatmap itself begins, and 'd' is where the heatmap itself ends.
The vertical limits are computed according to 'ylim' to align correctly with the heatmap rectangles.
}
  \item{toplim}{
The vertical limits of the row hierarchical clustering. The horizontal limits of the heatmap are {a=ylim[1], b=ylim[2], c=toplim[1], d=toplim[2]} where the 'a' is where the heatmap begins, 'b' is where the heatmap ends, 'c' is where the column dendrogram begins, and 'd' is where the column dendrogram ends.
The horizontal limits are computed according to 'xlim' to align correctly with the heatmap rectangles.
}
  \item{rightlim}{
The horizontal limits for the row texts.
}
  \item{bottomlim}{
The vertical limits for the column texts.  
}
  \item{type}{
Type of clustering visualization; while default "rect" provides a rectangular-angled tree, the alternate option "tri" provides triangular-angled tree.
}
  \item{scale}{
Should data be scaled according to 'row' or 'column' (or 'col') similarly to 'stats::heatmap'.
}
  \item{na.rm}{
Should missing values be removed, by default TRUE.  
}
  \item{nbins}{
Number of discrete bins the data is divided into using 'seq(from=min(x), to=max(x), length.out=nbins)'.  
}
  \item{valseq}{
The sequence of values by which the data is binned, typically corresponding to the above 'nbins' parameter.
}
  \item{namerows}{
If a single boolean value TRUE, then the default 'rownames(x)' are plotted to the right of the rows with the 'textfun'. If it is a vector of length 'nrow(x)', then this vector is used for plotting the row names instead.
}
  \item{namecols}{
If a single boolean value TRUE, then the default 'colnames(x)' are plotted below the columns with the 'textfun'. If it is a vector of length 'nrow(x)', then this vector is used for plotting the column names instead.
}
  \item{\dots}{
Additional parameters provided to the rectangle plotting function
}
	\if{html}{\figure{hmapex.png}{options: width="35\%" alt="Figure: mai.png"}}
	\if{latex}{\figure{hmapex.pdf}{options: width=7cm}}
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\seealso{
\code{\link{heatmap}}
\code{\link{hmap.key}}
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

\dontrun{
	# Heatmap used as base for the help documentation figure
	set.seed(1)
	hmap(matrix(rnorm(100), nrow=10), xlim=c(0.2,0.8), ylim=c(0.2,0.8), 
		leftlim=c(0.0,0.2), rightlim=c(0.8,1.0), 
		bottomlim=c(0.0,0.2), toplim=c(0.8,1.0))
	axis(1); axis(2); title(xlab="x", ylab="y")
}
}
\keyword{ hplot }
\keyword{ aplot }
