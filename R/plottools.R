# Scatterplot function for mixed type (numerical and categorical) variables with marginal distributions in each scatterplot
mixplot <- function(
	# Data frame or a matrix
	x,
	# Main title
	main = NA,
	# If legend should be constructed and plotted
	legend = T,
	# Colors for the categories or single observations if categories are not present
	col = palette(),
	# Should lines be drawn to represent one of the variables if the other one is missing in a 2-dim scatterplot
	na.lines = T,
	# Plot origin {0,0} using vertical and horizontal ablines
	origin = F,
	# Should marginal distributions be plotted, value of 0/FALSE/"no"/"none", 1/TRUE/"rug", 2/"hist"
	marginal = F,
	# Layout heights
	lhei,
	# Layout widths
	lwid,
	# Level of verbosity: -1<= (no verbosity), 0/FALSE (warnings) or >=1/TRUE (additional information)
	verb = 0,
	# Additional parameters
	...
){
	# Old par-settings (taken from ?layout)
	def.par <- par(no.readonly = TRUE) # save default
	
	# Make sure that some of the parameters are of correct class or cast to one
	verb <- as.numeric(verb)
	marginal <- as.character(marginal)
	if(is.vector(x)) x <- as.matrix(x)

	# Which fields to consider numeric or categorical (treated differently)
	nums <- c("numeric", "integer", "logical")
	ctgr <- c("factor", "character")
	
	# Extract types for the columns, cast x to data.frame if necessary
	if(class(x)=="data.frame"){
		colclass <- lapply(x, FUN=class)
	}else if(class(x)=="matrix"){
		colclass <- class(x[1,1])
	}else{
		x <- as.data.frame(x)
		colclass <- lapply(x, FUN=class)
	}

	# Number of matrix elements in the plotting region ncol_numeric x ncol_numeric
	wnum <- colclass %in% nums
	nnum <- sum(wnum)
	# Number of different categorical variables
	wcat <- colclass %in% ctgr
	ncat <- sum(wcat)

	# If no categorical variables available, introduce artificial category "Observation" which affects all
	if(ncat==0){
		cats <- rep("Observation", times=nrow(x))
		wcats <- 1:nrow(x)
	}else{
	# Else map observations to categories based on ctgr-type columns
		cats <- apply(x[,wcat,drop=F], MARGIN=1, FUN=function(x) paste(x, collapse=", "))
		wcats <- match(cats, unique(cats))
	}
	# Colors for the categories per each observation
	col <- rep(col, length.out=length(wcats))[wcats]

	# If user has not defined alternative lhei (layout heights) and lwid (layout widths) construct them
	if(missing(lhei)){
		lhei <- c()
		if(!is.na(main)) lhei <- c(lhei, 0.1)
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			lhei <- c(lhei, rep(c(0.05,0.3), times=ifelse(nnum>2, nnum, 1)))
		}else{
			lhei <- c(lhei, rep(0.3, times=ifelse(nnum>2, nnum, 1)))
		}
		if(legend) lhei <- c(lhei, 0.1)
	}
	if(missing(lwid)){
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			lwid <- rep(c(0.3, 0.05), times=ifelse(nnum>2, nnum, 1))
		}else{
			lwid <- rep(0.3, times=ifelse(nnum>2, nnum, 1))
		}
	}

	# Division of device plot region to subregions
	subsize <- ifelse(!as.character(marginal) %in% c("0", "FALSE", "no", "none"), 4, 1)
	if(subsize==4){
		temp <- rbind(rep(letters[1:2], times=ifelse(nnum>2, nnum, 1)), rep(letters[3:4], times=ifelse(nnum>2, nnum, 1)))[rep(1:2, times=ifelse(nnum>2, nnum, 1)),]
		temp[temp=="a"] <- seq(from=1+0, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+0, by=4)
		temp[temp=="b"] <- seq(from=1+1, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+1, by=4)
		temp[temp=="c"] <- seq(from=1+2, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+2, by=4)
		temp[temp=="d"] <- seq(from=1+3, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+3, by=4)
		class(temp) <- "numeric"
		lmat <- rbind( 
				rep(1, times=length(lwid)),
				t(temp)+1,
				rep(1+max(temp)+1, times=length(lwid))
			)
	}else{
		lmat <- rbind( 
				rep(1, times=length(lwid)),
				matrix(2:(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1)+1), ncol=nnum, byrow=T),
				rep(1+ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1)+1, times=length(lwid))
			)
	}
	# Remove header row from layout matrix if it is not desired
	if(is.na(main)){
		lmat <- lmat - 1
		lmat <- lmat[-1,]
	}
	# Remove legend row from layout matrix if it is not desired
	if(!legend){
		lmat <- lmat[-nrow(lmat),]
	}
	if(as.numeric(verb)>=1){
		print("lmat"); print(lmat)
		print("lwid"); print(lwid)
		print("lhei"); print(lhei)
	}
	# Set up the layout matrix to the plot device
	l <- layout(lmat, widths=lwid, heights=lhei)

	# If we should plot the title
	if(!is.na(main)){
		par(mar=c(0,0,0,0)); plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))
		text(0,0, main)
	}
	
	# Sub plotting function for each 1x1 or 2x2 subplot depending on marginal type
	subplox <- function(x, marginal, labels = F, ...){
		# Should variable labels be shown (and thus larger margins)
		if(labels){
			mar <- c(4,4,1,1)
		}else{
			mar <- c(2,2,1,1)
		}
		# Top marginal distribution
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			martemp <- mar; martemp[c(1,3)] <- 0; par(mar=martemp)			
			# Marginals
			plot.new(); plot.window(xlim=extendrange(x[,1]), ylim=c(-1,1))
			# Rug
			if(any(as.character(marginal) %in% c("1", "TRUE", "rug"))){
				for(i in 1:nrow(x)) abline(v=x[i,1], col=col[i], lwd=2)
			}
			# Histogram
			if(any(as.character(marginal) %in% c("2", "hist"))){
				uniqs <- unique(wcats)
				hall <- hist(x[,1], plot=F, breaks=seq(from=min(x[,1], na.rm=T), to=max(x[,1], na.rm=T), length.out=50))
				for(i in 1:length(uniqs)){
					h <- hist(x[wcats==uniqs[i],1], plot=F, breaks=seq(from=min(x[,1], na.rm=T), to=max(x[,1], na.rm=T), length.out=50))
					for(j in 1:length(h$density)){
						d <- h$counts[j]/max(hall$counts)
						rect(xleft = h$breaks[j], xright = h$breaks[j+1], ybottom = -1, ytop = unlist(ifelse(d>0, d, -1)), col=unique(col)[i], border=NA)
					}
				}
			}
		}
		# Actual scatterplot
		par(mar=mar)
		plot(x, xlim=extendrange(x[,1]), ylim=extendrange(x[,2]), ...)
		# Origin lines if desired
		if(origin){
			abline(h=0, col="grey")
			abline(v=0, col="grey")
		}
		# Top-right empty box and right marginal distribution
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			# Empty box top-right
			par(mar=c(0,0,0,0))
			plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))
			martemp <- mar; martemp[c(2,4)] <- 0; par(mar=martemp)			
			# Right marginal
			plot.new(); plot.window(xlim=c(-1,1), ylim=extendrange(x[,2]))
			# Rug
			if(any(as.character(marginal) %in% c("1", "TRUE", "rug"))){
				for(i in 1:nrow(x)){ abline(h=x[i,2], col=col[i], lwd=2)}
			}
			# Histogram
			if(any(as.character(marginal) %in% c("2", "hist"))){
				uniqs <- unique(wcats)
				hall <- hist(x[,2], plot=F, breaks=seq(from=min(x[,2], na.rm=T), to=max(x[,2], na.rm=T), length.out=50))
				for(i in 1:length(uniqs)){
					h <- hist(x[wcats==uniqs[i],2], plot=F, breaks=seq(from=min(x[,2], na.rm=T), to=max(x[,2], na.rm=T), length.out=50))
					for(j in 1:length(h$density)){
						d <- h$counts[j]/max(hall$counts)
						rect(ytop = h$breaks[j], ybottom = h$breaks[j+1], xleft = -1, xright = unlist(ifelse(d>0, d, -1)), col=unique(col)[i], border=NA)
					}
				}
			}
		}
	}
	
	# Plot actual scatterplots
	# If only one numeric column handle this exceptionally
	if(nnum==1){
		subplox(x=x[,which(wnum)[1]], marginal=marginal, labels = T, col=col, ...)
	# Else multiple scatterplots
	}else if(nnum==2){
		subplox(x=x[,which(wnum)[1:2]], marginal=marginal, labels = T, col=col, ...)
	}else{
		for(row in 1:nnum){
			for(column in 1:nnum){
				if(verb>=1) print(paste("row", row, "column", column))
				if(row==column){
					par(mar=c(0,0,0,0))
					# Empty top marginal
					if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
						plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
					}
					# Text of the variable
					plot.new()
					plot.window(xlim=c(-1,1), ylim=c(-1,1))
					text(0,0,colnames(x)[which(wnum)[row]])
					# Empty top right box and right marginal
					if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
						plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
						plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
					}
				}else{
					subplox(x=x[,rev(which(wnum)[c(row,column)])], marginal=marginal, col=col, ...)
				}
			}
		}
	}
	
	# Possibly build the legend and plot it	
	if(!as.numeric(legend)==0){
		par(mar=c(0,0,0,0))
		plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
		if(verb>=1){
			print(col)
			print(wcats)
			print(cats)
		}
		legend("bottom", col=unique(col)[unique(wcats)], legend=unique(cats)[unique(wcats)], horiz=T, bty="n", ...)
	}
	
	# Restore old par settings
	par(def.par)
	
	# Return invisibly some information about the plot
	invisible(list(x=x, lmat=lmat, lwid=lwid, lhei=lhei))
}

# Plot-region based heatmap
hmap <- function(
	# Input data matrix to plot
	x,

	# Plotting region settings
	#
	# Whether we want to add to an existing region or create a new one
	add = F,
	# x and y axis limits in the heatmap
	xlim=c(0,1),
	ylim=c(0,1),
	
	# Colors to use
	col = heat.colors(10),
	
	#
	# Border Settings for the heatmap bins
	# Should be a matrix of equal dimensions to x
	#
	# Color for borders in the heatmap bins
	border = matrix(NA, nrow=nrow(x), ncol=ncol(x)),
	# Line type for borders in the heatmap bins
	lty = matrix("solid", nrow=nrow(x), ncol=ncol(x)),
	# Line width for borders in the heatmap bins
	lwd = matrix(1, nrow=nrow(x), ncol=ncol(x)),
	
	# Additional parameters
	...
){
	# Number of rows and columns
	nr = nrow(x)
	nc = ncol(x)
	
	# Recycle border color, line type and line width settings if dimensions are not equal to x
	if(!identical(dim(border),dim(x))){
		print("Border not identical")
		border = matrix(border, nrow=nrow(x), ncol=ncol(x))
	}
	if(!identical(dim(lty),dim(x))){
		print("lty not identical")
		lty = matrix(lty, nrow=nrow(x), ncol=ncol(x))
	}
	if(!identical(dim(lwd),dim(x))){
		print("lwd not identical")
		lwd = matrix(lwd, nrow=nrow(x), ncol=ncol(x))
	}
	
	# Value range
	nbins = length(col)
	valseq = seq(from=min(x, na.rm=T), to=max(x, na.rm=T), length.out=nbins)
	# Finding color bins for the values according to the interval number in palette
	intervalmat = matrix(findInterval(x=x, vec=valseq), nrow=nr, ncol=nc)
	
	xseq = seq(from=xlim[1], to=xlim[2], length.out=nc+1)
	# Reversing the y-coordinates, otherwise the heatmap will be upside-down
	yseq = rev(seq(from=ylim[1], to=ylim[2], length.out=nr+1))
	
	# If requested, create a new plotting region
	if(!add){
		plot.new()
		plot.window(xlim=xlim, ylim=ylim)
	}
	
	# Plotting the rectangles
	for(i in 1:(nr-1)){
		for(j in 2:nc){
			rect(xleft=xseq[i], ybottom=yseq[j], xright=xseq[i+1], ytop=yseq[j-1], col=col[intervalmat[i,j]], border=border[i,j], lty=lty[i,j], lwd=lwd[i,j])
		}
	}
}

# Extend range -function (based on 'extendrange' from grDevices) with forced symmetry around a point
extendsymrange <- function(
	x, # extendrange x
	r = range(x, na.rm=T), # extendrange r
	f = 0.05, # extendrange f
	sym = 0 # Axis of symmetry
){
	ex <- extendrange(x = x, r = r, f = f)
	ran <- max(abs(ex - sym))
	c(sym - ran, sym + ran)
}

# Give in vector (or matrix/data.frame) of y-values, and obtain smart jittering for separating the measurements on the x-axis (obtain x-axis values)
smartjitter <- function(
	# Original data vector (or matrix or data.frame)
	x,
	# Quantile splits
	q = seq(from=0, to=1, length.out=10),
	# Type of jittering;
	# type == 1: consecutive jitters in same bin are values 0,0.1,0.2,0.3, ...
	# type != 1: consecutive jitters in the same bin are alternating -1 coefficiented values 0,-0.1,0.2,-0.3,0.4, ...
	type = 1,
	# Amount of jittering per overlapping values
	amount = 0.1,
	# Jittering function for values that reside in the same split bin
	jitterfuncs = list( 
		# type == 1 option function
		function(n){
			(1:n)/(1/amount)
		},
		# type != 1 option function
		function(n){
			(((-1)^c(0:(n-1)))*(0:(n-1)))/(1/amount)
		}
	),
	jits = jitterfuncs[[type]]
	)
{
	w <- as.numeric(cut(x, breaks=c(-Inf,unique(quantile(x, probs=q, na.rm=T)),Inf)))
	jit <- unlist(
		apply(
			do.call("rbind", 
				lapply(unique(w), FUN=function(z) 
					{ 
						res <- rep(NA, times=length(w)); 
						res[which(w==z)] <- jits(length(which(w==z))); 
						res 
					} 
				)), 
		MARGIN=2, FUN=function(y) 
			ifelse(all(is.na(y)), 0 , y[!is.na(y)]))
	)
	jit
}
