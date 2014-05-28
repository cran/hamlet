# Plot residuals of a mixed-effects model along with linear trend fit
mem.plotresid <- function(fit, linear=T, type="XbZu", main, xlab, ylab){
	res <- mem.getcomp(fit)
	plot(res[,type], res[,"e"], pch=16, xlab="", ylab="")
	
	if(linear){
		residlm <- lm(e ~ y, data.frame(y = res[,type], e = res[,"e"]))
		abline(residlm)	
		legend("bottomright", lwd=1, col="black", legend=c("Linear trend"), bty="n")	
	}
	if(!missing(main)) title(main=main)
	if(!missing(xlab)){ title(xlab=xlab)}
	else{ title(xlab=type)}
	if(!missing(ylab)){ title(ylab=ylab)}
	else{ title(ylab="e")}
}

# Plot histogram distributions of random effects
mem.plotran <- function(fit, breaks=100){
	count <- 0
	if(class(fit) %in% c("lmerMod", "mer", "merModLmerTest")){
		lapply(ranef(fit), FUN=function(x) count <<- count + ncol(x))
		par(mfrow=c(ceiling(sqrt(count)),ceiling(sqrt(count))))
		lapply(ranef(fit), FUN=function(x) apply(x, MARGIN=2, FUN=function(y) hist(y, breaks=breaks)))
	}else if(class(fit) %in% c("lme", "merModLmerTest")){
		count <- count + ncol(ranef(fit))
		par(mfrow=c(ceiling(sqrt(count)),ceiling(sqrt(count))))
		apply(ranef(fit), MARGIN=2, FUN=function(y) hist(y, breaks=breaks))
	}else{
		stop(paste("Invalid class of fit:", class(fit)))
	}
	
}

# Get per-observation components Xb and Zu from a mixed-effects model
mem.getcomp <- function(fit){
	if(class(fit) %in% c("lmerMod", "glmerMod", "mer", "merModLmerTest")){
		# Xb 
		fix <- getME(fit,'X') %*% fixef(fit)
		# Zu
		ran <- t(as.matrix(getME(fit,'Zt'))) %*% unlist(ranef(fit))
		# Xb + Zu
		fixran <- getME(fit,'y') - resid(fit)
		# y
		y = getME(fit, 'y')
	}else if(class(fit)=="lme"){
		# Xb
		fix <- fit$fitted[,1]
		# Xb + Zu
		fixran <- fit$fitted[,2]
		# Zu = Xb + Zu - Xb
		ran <- fixran - fix
		# y
		y <- fit$data[,as.character(fit$terms[[2]])]
	}else if(class(fit)=="lm"){
		fix = fixran = fitted(fit)
		ran = 0
		y <- fit$model[,as.character(fit$terms[[2]])]
	}else{
		stop(paste("Invalid class of fit:", class(fit)))
	}
	result <- cbind(fix = fix, ran = ran, fixran = fixran, e = resid(fit), y = y)
	colnames(result) <- c("Xb", "Zu", "XbZu", "e", "y")
	result
}
