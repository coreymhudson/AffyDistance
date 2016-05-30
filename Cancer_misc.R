library(MASS)
library(stat4)
library(truncgof)

normal = "normal.txt"

read_in_expression <- function(filename){
	mat <- as.matrix(read.table(filename, as.is=T))
	mat <- mat[-1,]
	class(mat) <- "numeric"
	return(mat)
}

normal_matrix = read_in_expression(normal)
vec_beta <- numeric(length(normal_matrix[,1]))
for(i in 1:length(vec_beta)){
	v <- normal_matrix[i,]
	k <- ks.test(v, "pbeta", list(shape1=mean(v)*((mean(v)*(1-mean(v)))/(var(v))-1), shape2=(1-mean(v))*((mean(v)*(1-mean(v)))/(var(v))-1)))
	vec_beta[i] <- k$p.value
}
write.table(filename="Normal_beta_stats.txt", summary(vec_beta))
postscript(file="Beta_hist_normal.eps")
hist(vec_beta, main="KS_pvalues", xlab="pvalues")
dev.off()
stat = numeric(3)
stat[1] = length(vec_beta[vec_beta<0.01])/length(vec_beta)
stat[2] = length(vec_beta[vec_beta<0.05])/length(vec_beta)
stat[3] = length(vec_beta[vec_beta<0.1])/length(vec_beta)
write.table(filename="Normal_beta_error.txt", stat)

ML_null <- function(x_cancer, x_normal, x_embryonic){
	x_total <- c(x_cancer, x_normal, x_embryonic)
	v <- x_total
	total_shape1 <- mean(v)*((mean(v)*(1-mean(v)))/(var(v))-1)
	total_shape2 <- (1-mean(v))*((mean(v)*(1-mean(v)))/(var(v))-1)
	null_estimation <- mle(minuslog=null_model, start=list(shape1=total_shape1, shape2=total_shape2), fixed=list(x_cancer=x_cancer, x_normal=x_normal, x_embryonic=x_embryonic))
	null_estimation
}

ML_null(tumor_matrix[1,], normal_matrix[1,], embryonic_matrix[1,])

full_model <- function(x_cancer, x_normal, x_embryonic, shape1_cancer, shape1_normal, shape1_embryonic, shape2_cancer, shape2_normal, shape2_embryonic){
	(-sum(dbeta(x_cancer, shape1=shape1_cancer, shape2=shape2_cancer, log=T)))+(-sum(dbeta(x_normal, shape1=shape1_normal, shape2=shape2_normal, log=T)))+(-sum(dbeta(x_embryonic, shape1=shape1_embryonic, shape2=shape2_embryonic, log=T)))
}

ML_full <- function(x_cancer, x_normal, x_embryonic){
	a <- x_cancer
	b <- x_normal
	c <- x_embryonic
	a_shape1 <- mean(a)*((mean(a)*(1-mean(a)))/(var(a))-1)
	a_shape2 <- (1-mean(a))*((mean(a)*(1-mean(a)))/(var(a))-1)
	b_shape1 <- mean(b)*((mean(b)*(1-mean(b)))/(var(b))-1)
	b_shape2 <- (1-mean(b))*((mean(b)*(1-mean(b)))/(var(b))-1)
	c_shape1 <- mean(c)*((mean(c)*(1-mean(c)))/(var(c))-1)
	c_shape2 <- (1-mean(c))*((mean(c)*(1-mean(c)))/(var(c))-1)
	full_estimation <- mle(minuslog=full_model, start=list(shape1_cancer=a_shape1, shape2_cancer=a_shape2, shape1_normal=b_shape1, shape2_normal=b_shape2, shape1_embryonic=c_shape1, shape2_embryonic=c_shape2), fixed=list(x_cancer=x_cancer, x_normal=x_normal, x_embryonic=x_embryonic))
	full_estimation
}

LRT_full_null <- function(x_cancer, x_normal, x_embryonic){
	a <- x_cancer
	b <- x_normal
	c <- x_embryonic
	NL <- ML_null(a, b, c)
	NLL <- logLik(NL)
	FL <- ML_full(a, b, c)
	FLL <- logLik(FL)
	D = (-2*NLL)-(-2*FLL)
	return(pchisq(D[1], 2))
}

vec <- double(length(tumor_matrix[,1]))
for(i in 1:length(tumor_matrix[,1]))
{
	vec[i] <- try(LRT_full_null(tumor_matrix[i,], normal_matrix[i,], embryonic_matrix[i,]))
	cat(i, "\t", vec[i], "\n")
}

shape_best <- function(x){
	v <- x
	z <- (2)
	z[1] <- mean(v)*((mean(v)*(1-mean(v)))/(var(v))-1)
	z[2] <- (1-mean(v))*((mean(v)*(1-mean(v)))/(var(v))-1)
	return(z)
}

pair_model <- function(x, y, shape1_x, shape2_x, shape1_y, shape2_y){
	(-sum(dbeta(x, shape1=shape1_x, shape2=shape2_x, log=T)))+(-sum(dbeta(y, shape1=shape1_y, shape2=shape2_y, log=T)))
}


ML_pairs <- function(x, y){
	a <- x
	b <- y
	a_shape1 <- mean(a)*((mean(a)*(1-mean(a)))/(var(a))-1)
	a_shape2 <- (1-mean(a))*((mean(a)*(1-mean(a)))/(var(a))-1)
	b_shape1 <- mean(b)*((mean(b)*(1-mean(b)))/(var(b))-1)
	b_shape2 <- (1-mean(b))*((mean(b)*(1-mean(b)))/(var(b))-1)
	pair_estimation <- mle(minuslog=pair_model, start=list(shape1_x=a_shape1, shape2_x=a_shape2, shape1_y=b_shape1, shape2_y=b_shape1), fixed=list(x=x, y=y))
	pair_estimation
}

MedianVector_tumor <- sapply(tumor_matrix, median)

Check: Values are in the right rows.

Assumptions:
	Between normal and tumor a limited number of genes change in expression.
	
vec <- numeric(length(normal_matrix[,1]))
for(i in 1:length(vec)){
	print(i)
	k <- ks.test(normal_matrix[i,], "pnorm", list(mean=mean(normal_matrix[i,]), sd=sd(normal_matrix[i,])))
	vec[i] <- k$p.value
}

vec_exp <- numeric(length(normal_matrix[,1]))
for(i in 1:length(vec_exp)){
	print(i)
	k <- ks.test(normal_matrix[i,], "pexp", list(rate=fitdistr(normal_matrix[i,], "exponential")$estimate[1]))
	vec_exp[i] <- k$p.value
}

vec_beta <- numeric(length(normal_matrix[,1]))
for(i in 1:length(vec_beta)){
	print(i)
	v <- normal_matrix[i,]
	k <- ks.test(v, "pbeta", list(shape1=mean(v)*((mean(v)*(1-mean(v)))/(var(v))-1), shape2=(1-mean(v))*((mean(v)*(1-mean(v)))/(var(v))-1)))
	vec_beta[i] <- k$p.value
	

optimize_move <- function(x, data1, data2){
	-log(ks.test(data1+x, data2)$p.value)
}

move_distribution <- function(x, data1, data2){
	f1 <- fitdistr(data1, "norm")
	f2 <- fitdistr(data2, "norm")
	mu_x <- f1$estimate[1]-f2$estimate[1]
	sd_x <- sqrt((f1$estimate[2]^2)+(f2$estimate[2]^2))
	p <- suppressWarnings(optim(as.numeric(mu_x), optimize_move, method="Nelder-Mead", hessian=T, data1=data1, data2=data2))
}



optimize_move <- function(x){
	-log(ks.test(normal_data[1,]+x, cancer_data[1,])$p.value)
}
