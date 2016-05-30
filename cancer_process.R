source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")

	

load_expression_libraries <- function()
{
	library(limma)
	library(gcrma)
	library(MASS)
	library(gamlss)
	library(stats4)
	require(R.utils)
}

expression_preprocess <- function(dir, outputfile=NULL)
{
	setwd(dir)
	ab <- ReadAffy()
	eset <- rma(ab)
	e_intensity <- exprs(eset)
	e_samples <- length(e_intensity[1,])
	e_spots <- length(e_intensity[,1])
	transformed_exp <- apply(e_intensity, 2, function(e_intensity){s <- sum(e_intensity); v<- e_intensity/s; v})
	if(!missing(outputfile)){
		write.matrix(transformed_exp, file=outputfile)
	}
	transformed_exp
}

expression_preprocess <- function(rmas, outputfile=NULL)
{
	e_intensity <- exprs(rmas)
	e_samples <- length(e_intensity[1,])
	e_spots <- length(e_intensity[,1])
	transformed_exp <- apply(e_intensity, 2, function(e_intensity){s <- sum(e_intensity); v<- e_intensity/s; v})
	if(!missing(outputfile)){
		write.matrix(transformed_exp, file=outputfile)
	}
	transformed_exp
}


null_model <- function(x_cancer, x_normal, x_embryonic, shape1, shape2){
	x_total <- c(x_cancer, x_normal, x_embryonic)
	-sum(dbeta(x_total, shape1=shape1, shape2=shape2, log=TRUE))
}

difference_in_regulation <- function(x_cancer, x_normal, x_embryonic, shape1_cancer, shape1_normal, shape1_embryonic, shape2_cancer, shape2_normal, shape2_embryonic){
	-sum(dbeta(x_cancer, shape1=shape1_cancer, shape2=shape2_cancer))-(-sum(dbeta(x_normal, shape1=shape1_normal, shape2=shape2_normal))-(-sum(dbeta(x_embryonic, shape1=shape1_embryonic, shape2=shape2_embryonic))
}

Different_RatesLRT <- function(x_cancer, x_normal){
	x_total <- c(x_cancer, x_normal)
	starting_point_total_shape1 <- mean(x_total)*((mean(x_total)*(1-mean(x_total))/(var(x_total))-1)
	starting_point_total_shape2 <- (1-mean(x_total))*((mean(x_total)*(1-mean(x_total))/(var(x_total))-1)
	starting_point_cancer_shape1 <- mean(x_cancer)*((mean(x_cancer)*(1-mean(x_cancer))/(var(x_cancer))-1)
	starting_point_cancer_shape2 <- (1-mean(x_cancer))*((mean(x_cancer)*(1-mean(x_cancer))/(var(x_cancer))-1)
	starting_point_normal_shape1 <- mean(x_normal)*((mean(x_normal)*(1-mean(x_normal))/(var(normal))-1)
	starting_point_normal_shape2 <- (1-mean(x_normal))*((mean(x_normal)*(1-mean(x_normal))/(var(normal))-1)
	null_estimation <- mle(minuslog=null_model, start=list(shape1=starting_point_total_shape1, shape2=starting_point_total_shape2), fixed=list(x_cancer=x_cancer, x_normal=x_normal))
	change_estimation <- mle(minuslog_change_in_regulation, start=list(shape1_cancer=starting_point_cancer_shape1, shape2_cancer=starting_point_cancer_shape2, shape1_normal=starting_point_normal_shape1, shape2_normal=starting_point_normal_shape2), fixed=list(x_cancer=x_cancer, x_normal=x_normal))
	D = (-2*logLik(null_estimation))-(-2*logLik(change_estimation))
	print(D)
	pchisq(D[1], 2)
}


testNormEXPWithKS <- function(count, trans_exp){
	pvalue <- double(count)
	for(i in 1:count){
		print(i)
		x <- trans_exp[i,]
		fit <- try(normexp.fit(x))
		if(any(pexGAUS(x, fit$par[1], exp(fit$par[2]), exp(fit$par[3]))<0)){pvalue[i] <- NA}
		else{
			pvalue[i] <- try(ks.test(x, "pexGAUS", fit$par[1], exp(fit$par[2]), exp(fit$par[3]))$p.value)
		}
	}
	return(pvalue)
}

testNormalWithKS <- function(count, trans_exp){
	pvalue <- double(count)
	for(i in 1:count){
		print(i)
		x <- trans_exp[i,]
		fit <- try(fitdistr(x, "normal"))
		if(inherits(fit, "try-error") == TRUE){
			pvalue[i] <- NA
		}
		else if(!all(pgamma(x, fit$estimate[1], fit$estimate[2])>0)){pvalue[i] <- NA}
		else{
			pvalue[i] <- try(ks.test(x, "pnorm", fit$estimate[1], fit$estimate[2])$p.value)
		}
	}
	return(pvalue)
}

normal_ps <- testNormalWithKS(e_spots, transformed_exp)

testLogNormalWithKS <- function(count, trans_exp){
	pvalue <- double(count)
	for(i in 1:count){
		print(i)
		x <- trans_exp[i,]
		fit <- try(fitdistr(x, "lognormal"))
		if(inherits(fit, "try-error") == TRUE){
			pvalue[i] <- NA
		}
		else if(!all(plnorm(x, fit$estimate[1], fit$estimate[2])>0)){pvalue[i] <- NA}
		else{
			pvalue[i] <- try(ks.test(x, "plnorm", fit$estimate[1], fit$estimate[2])$p.value)
		}
	}
	return(pvalue)
}

lognormal_ps <- testLogNormalWithKS(e_spots, transformed_exp)
summary(lognormal_ps)
hist(lognormal_ps)
length(lognormal_ps[lognormal_ps<0.05])/length(lognormal_ps)

testExponentialWithKS <- function(count, trans_exp){
	pvalue <- double(count)
	for(i in 1:count){
		print(i)
		x <- trans_exp[i,]
		fit <- try(fitdistr(x, "exponential"))
		if(inherits(fit, "try-error") == TRUE){
			pvalue[i] <- NA
		}
		else if(!all(pexp(x, fit$estimate[1])>0)){pvalue[i] <- NA}
		else{
			pvalue[i] <- try(ks.test(x, "pexp", fit$estimate[1])$p.value)
		}
	}
	return(pvalue)
}

exponential_ps <- testExponentialWithKS(e_spots, transformed_exp)
summary(exponential_ps)
hist(exponential_ps)
length(exponential_ps[exponential_ps<0.05])/length(exponential_ps)

testGammaWithKS <- function(count, trans_exp){
	pvalue <- double(count)
	for(i in 1:count){
		print(i)
		x <- trans_exp[i,]
		fit <- try(fitdistr(x, "gamma"))
		if(inherits(fit, "try-error") == TRUE){
			pvalue[i] <- NA
		}
		else if(!all(pgamma(x, fit$estimate[1], fit$estimate[2])>0)){pvalue[i] <- NA}
		else{
			pvalue[i] <- try(ks.test(x, "pgamma", fit$estimate[1], fit$estimate[2])$p.value)
		}
	}
	return(pvalue)
}

testGammaWithKS <- function(count, trans_exp){
	pvalue <- double(count)
	for(i in 1:count){
		print(i)
		x <- trans_exp[i,]
		fit <- try(fitdistr(x, "gamma"))
		print(fit)
	}
	return(pvalue)
}

gamma_ps <- testGammaWithKS(e_spots, transformed_exp)

evalWithTimeout(testGammaWithKS(e_spots, transformed_exp) timeout=6, onTimeout="warning")

evalWithTimeout(testWithKS(e_spots, transformed_exp, timeout=6, onTimeout="warning")

null_model <- function(x_cancer, x_normal, shape1, shape2){
	x_total <- c(x_cancer, x_normal)
	-sum(dbeta(x_total, shape1=shape1, shape1=shape2))
}

change_in_regulation <- function(x_cancer, x_normal, shape1_cancer, shape1_normal, shape2_cancer, shape2_normal){
	-sum(dbeta(x_cancer, shape1=shape1_cancer, shape2=shape2_cancer))-(-sum(dbeta(x_normal, shape1=shape1_normal, shape2=shape2_normal))
}

Different_RatesLRT <- function(x_cancer, x_normal){
	x_total <- c(x_cancer, x_normal)
	starting_point_total_shape1 <- mean(x_total)*((mean(x_total)*(1-mean(x_total))/(var(x_total))-1)
	starting_point_total_shape2 <- (1-mean(x_total))*((mean(x_total)*(1-mean(x_total))/(var(x_total))-1)
	starting_point_cancer_shape1 <- mean(x_cancer)*((mean(x_cancer)*(1-mean(x_cancer))/(var(x_cancer))-1)
	starting_point_cancer_shape2 <- (1-mean(x_cancer))*((mean(x_cancer)*(1-mean(x_cancer))/(var(x_cancer))-1)
	starting_point_normal_shape1 <- mean(x_normal)*((mean(x_normal)*(1-mean(x_normal))/(var(normal))-1)
	starting_point_normal_shape2 <- (1-mean(x_normal))*((mean(x_normal)*(1-mean(x_normal))/(var(normal))-1)
	null_estimation <- mle(minuslog=null_model, start=list(shape1=starting_point_total_shape1, shape2=starting_point_total_shape2), fixed=list(x_cancer=x_cancer, x_normal=x_normal))
	change_estimation <- mle(minuslog_change_in_regulation, start=list(shape1_cancer=starting_point_cancer_shape1, shape2_cancer=starting_point_cancer_shape2, shape1_normal=starting_point_normal_shape1, shape2_normal=starting_point_normal_shape2), fixed=list(x_cancer=x_cancer, x_normal=x_normal))
	D = (-2*logLik(null_estimation))-(-2*logLik(change_estimation))
	print(D)
	pchisq(D[1], 2)
}

library(limma)
library(gcrma)
setwd("~/Desktop/cancer_R/GSE23402_RAW/")
ab <- ReadAffy()
eset <- rma(ab)
e_intensity <- exprs(eset)
e_samples <- length(e_intensity[1,])
e_spots <- length(e_intensity[,1])
transformed_exp <- apply(e_intensity, 2, function(e_intensity){ s <- sum(e_intensity); v <- e_intensity/s; v})
library(MASS)
library(gamlss)
require(R.utils)

testWithAD <- function(count, trans_exp, distribution){
	pvalue <- double(count)
	if(distribution == "exponential"){
		for(i in 1:count){
			x <- trans_exp[i,]
			fit <- try(fitdistr(x, "exponential"))
			if(inherits(fit, "try-error") == TRUE){
				pvalue[i] <- NA
			} else if(!all(pexp(x, fit$estimate[1])>0)){pvalue[i] <- NA}
			else{
				pvalue[i] <- try(ad.test(x, pexp, fit$estimate[1])$p.value)
			}
		}
	}
	else if(distribution == "normal"){
		for(i in 1:count){
			x <- trans_exp[i,]
			fit <- try(fitdistr(x, "normal"))
			if(inherits(fit, "try-error") == TRUE){
				pvalue[i] <- NA
			} else if(!all(pnorm(x, fit$estimate[1], fit$estimate[2])>0)){pvalue[i] <- NA}
			else{
				pvalue[i] <- try(ad.test(x, pnorm, fit$estimate[1], fit$estimate[2])$p.value)
			}
		}
	}
	else if(distribution == "lognormal"){
		for(i in 1:count){
			x <- trans_exp[i,]
			fit <- try(fitdistr(x, "lognormal"))
			if(inherits(fit, "try-error") == TRUE){
				pvalue[i] <- NA
			} else if(!all(plnorm(x, fit$estimate[1], fit$estimate[2])>0)){pvalue[i] <- NA}
			else{
				pvalue[i] <- try(ad.test(x, plnorm, fit$estimate[1], fit$estimate[2])$p.value)
			}
		}
	}
	else if(distribution == "gamma"){
		for(i in 1:count){
			x <- trans_exp[i,]
			fit <- try(fitdistr(x, "gamma"))
			if(inherits(fit, "try-error") == TRUE){
				pvalue[i] <- NA
			} else if(!all(pgamma(x, fit$estimate[1], fit$estimate[2])>0)){pvalue[i] <- NA}
			else{
				pvalue[i] <- try(ad.test(x, pgamma, fit$estimate[1], fit$estimate[2])$p.value)
			}
		}
	}
	else if(distribution == "beta"){
		for(i in 1:count){
			x <- trans_exp[i,]
			alfa <- mean(x)*((mean(x)*(1-mean(x))/var(x))-1)
			b <- (1-mean(x))*((mean(x)*(1-mean(x))/var(x))-1)
			fit <- try(fitdistr(x, "beta", start=list(shape1=alfa, shape2=b)))
			if(inherits(fit, "try-error") == TRUE){
				pvalue[i] <- NA
			} else if(!all(pbeta(x, shape1=fit$estimate[1], shape2=fit$estimate[2])>0)){pvalue[i] <- NA}
			else{
				pvalue[i] <- try(ad.test(x, pbeta, shape1=fit$estimate[1], shape2=fit$estimate[2])$p.value)
			}
		}
	}
	return(pvalue)
}