###
# Teemu Daniel Laajala
# February 2016
# Confounder simulations - Supplementary Figure 9 draft
###

.simurun <- function(...){

# Number of simulations per each parameter combination
nsim <- 1000
# Combination grid of parameters
grid <- expand.grid(q=c(1,3,10,20), p=c(1,3,10,20), s=c(0,0.4,0.7), n=c(5,10,15), difference=c(0,1,2))
# For reproducibility
set.seed(1)
result <- list()
for(i in 1:nrow(grid)){
	cat("\n\n")
	print(paste(i, "of", nrow(grid)))
	print(Sys.time())
	cat("\n\n")
	result[[length(result)+1]] <- lapply(1:nsim, FUN=function(...){
		# Extract parameters
		q <- grid[i,"q"]
		p <- grid[i,"p"]
		s <- grid[i,"s"]
		n <- grid[i,"n"]
		dif <- grid[i,"difference"]
		# Construct Sigma
		covmat <- matrix(0, nrow=q+p, ncol=q+p)
		covmat <- cbind(rbind(covmat, c(s*rep(1, times=q), rep(0, times=p))), c(s*rep(1, times=q), rep(0, times=p), 1))
		diag(covmat) <- 1
		
		#library(MASS)
		#library(Matrix)
		# nearPD provides nearest positive definite covariance matrix with the desired properties
		basedat <- MASS::mvrnorm(n=n*2, mu=rep(0, times=q+p+1), Sigma=cov2cor(Matrix::nearPD(covmat)$mat))

		y <- basedat[,ncol(basedat)]
		x <- basedat[,-ncol(basedat)]
		d <- as.matrix(dist(x))
		library(hamlet)
		sink("bbtemp.txt")
		# Branch and Bound matching with constrained maximal branching
		m <- match.bb(d, g=2, maxbranches=100)
		sink()
		# Allocation based on the matching
		g_matched <- match.allocate(m$solution)
		# Allocation through conventional random permutation
		g_nomatch <- LETTERS[1:2][sample(c(rep(1, times=n), rep(2, times=n)))]
		# Conduct 3 different testing types
		y_matched <- y
		y_nomatch <- y
		y_matched[g_matched=="Group_A"] <- y_matched[g_matched=="Group_A"] + dif
		y_nomatch[g_nomatch=="A"] <- y_nomatch[g_nomatch=="A"] + dif
		matches <- by(g_matched, INDICES=m$solution, FUN=function(z){ 
			if(z[1]=="Group_A"){ c(1,2) }
			else{ c(2,1) }
		})
		y_paired <- do.call("rbind", by(1:length(m$solution), INDICES=m$solution, FUN=function(z){
			if(g_matched[z[1]]=="Group_A"){
				c(y_matched[z[1]], y_matched[z[2]])
			}else{
				c(y_matched[z[2]], y_matched[z[1]])
			}
		}))

		# Matched design, post-intervention pairing
		p1 <- t.test(y_paired[,1], y_paired[,2], paired=T)$p.value
		# Matched design, no post-intervention pairing
		p2 <- t.test(y_matched ~ g_matched)$p.value
		# Non-matched design, no post-intervention pairing
		p3 <- t.test(y_nomatch ~ g_nomatch)$p.value
		# Return the different p-values
		c(p1, p2, p3)
	})
}


### FINISHED RUNS - PLOT FIGURE DRAFT ###


# Average detection of effects
restab <- do.call("rbind", lapply(result, FUN=function(z){
	apply(do.call("rbind", z), MARGIN=2, FUN=function(q){
		sum(q<0.05)/length(q)
	})
}))
colnames(restab) <- c("MatchedPaired", "MatchedNonpaired", "ConventionalNonpaired")
# Bind results
tab <- cbind(grid, restab)

# Generate a Cairo PDF figure draft
#library(Cairo)
Cairo::CairoPDF("Confounder_Simulations.pdf", width=8, height=8)
l <- cbind(c(82,92:104), rbind(
	83,
	cbind( 0+matrix(1:9, byrow=T,ncol=3), 86,  9+matrix(1:9, byrow=T,ncol=3), 87, 18+matrix(1:9, byrow=T,ncol=3)),
	84,
	cbind(27+matrix(1:9, byrow=T,ncol=3), 88, 36+matrix(1:9, byrow=T,ncol=3), 89, 45+matrix(1:9, byrow=T,ncol=3)),
	85,
	cbind(54+matrix(1:9, byrow=T,ncol=3), 90, 63+matrix(1:9, byrow=T,ncol=3), 91, 72+matrix(1:9, byrow=T,ncol=3))
))
layout(l)
library(hamlet)
#par(mar=c(0,0,0,0), mfrow=c(9,9), cex.main=0.7)
par(mar=c(0,0,0,0), cex.main=0.7)
# par 1 = difference in groups
# par 2 = sample size n
# par 3 = predictiviness
# q = number of predictive baseline variables
# p = number of non-predictive baseline variables

subtablist <- list()

# Combine the palette of two weighted components
palet <- c(
	colorRampPalette(c("orange","red","black"), bias=4)(51),
	colorRampPalette(c("black", "blue", "cyan"), bias=0.4)(50)
)

methodvec <- c("MatchedPaired", "MatchedNonpaired", "ConventionalNonpaired")
for(par1 in unique(tab[,"difference"])){
	for(method in methodvec){
		for(par3 in unique(tab[,"s"])){
			for(par2 in unique(tab[,"n"])){
				subtab <- tab[tab[,"difference"]==par1 & tab[,"n"]==par2 & tab[,"s"]==par3,]
				subtablist[[length(subtablist)+1]] <- c(diff = par1, n = par2, s = par3, method = method)
				#print(subtab)
				submat <- matrix(subtab[,method], nrow=4)
				submat[1,1] <- subtab[subtab[,"q"]==1 & subtab[,"p"]==1,method]
				submat[2,1] <- subtab[subtab[,"q"]==3 & subtab[,"p"]==1,method]
				submat[3,1] <- subtab[subtab[,"q"]==10 & subtab[,"p"]==1,method]
				submat[4,1] <- subtab[subtab[,"q"]==20 & subtab[,"p"]==1,method]
				submat[1,2] <- subtab[subtab[,"q"]==1 & subtab[,"p"]==3,method]
				submat[2,2] <- subtab[subtab[,"q"]==3 & subtab[,"p"]==3,method]
				submat[3,2] <- subtab[subtab[,"q"]==10 & subtab[,"p"]==3,method]
				submat[4,2] <- subtab[subtab[,"q"]==20 & subtab[,"p"]==3,method]
				submat[1,3] <- subtab[subtab[,"q"]==1 & subtab[,"p"]==10,method]
				submat[2,3] <- subtab[subtab[,"q"]==3 & subtab[,"p"]==10,method]
				submat[3,3] <- subtab[subtab[,"q"]==10 & subtab[,"p"]==10,method]
				submat[4,3] <- subtab[subtab[,"q"]==20 & subtab[,"p"]==10,method]
				submat[1,4] <- subtab[subtab[,"q"]==1 & subtab[,"p"]==20,method]
				submat[2,4] <- subtab[subtab[,"q"]==3 & subtab[,"p"]==20,method]
				submat[3,4] <- subtab[subtab[,"q"]==10 & subtab[,"p"]==20,method]
				submat[4,4] <- subtab[subtab[,"q"]==20 & subtab[,"p"]==20,method]
				#print(submat)
				rownames(submat) <- paste("q=",c(1,3,10,20),sep="")
				colnames(submat) <- paste("p=",c(1,3,10,20),sep="")
				h <- hmap(submat, Colv=NA, Rowv=NA, valseq=seq(from=0.00, to=1.00, by=.01), col=palet, namerows=F, namecols=F, bottomlim=c(0,0.05), ylim=c(0,0.85), rightlim=c(0,0.05))
				title(main=paste("\n",par1,par2,par3,which(method==methodvec), collapse=","))
			}
		}
	}
}

plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
hmap.key(h, x0=0.1, x1=1, y0=0.1, y1=0.9, at=c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95))
text(0.5,0.05, "% significant", cex=0.7)

plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.2, "mu1 - mu2 = 0 (no difference)", cex=1)
	text(0.15,0.8, "Matched randomization\nPaired testing", cex=1.2)
	text(0.5,0.8, "Matched randomization\nNon-paired testing", cex=1.2)
	text(0.85,0.8, "Conventional randomization\nNon-paired testing", cex=1.2)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.2, "mu1 - mu2 = 1 (mild difference)", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.2, "mu1 - mu2 = 2 (strong difference)", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text5", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text6", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text7", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text8", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text9", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.5, "s=0\nNo predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.5, "s=0.4\nSlight predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
	text(0.5,0.5, "s=0.7\nStrong predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text13", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.5, "s=0\nNo predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.5, "s=0.4\nSlight predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
	text(0.5,0.5, "s=0.7\nStrong predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); #text(0.5,0.5, "text17", cex=1)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.5, "s=0\nNo predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); 
	text(0.5,0.5, "s=0.4\nSlight predictivity", cex=0.7)
plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
	text(0.5,0.5, "s=0.7\nStrong predictivity", cex=0.7)

subtablist <- do.call("rbind", subtablist)

# Close DPF
dev.off()

}