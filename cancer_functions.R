sumAIC <- function(z){
	sAIC = (2)
	maxex <- paste("sum(AIC(")
	for(i in z){
		r <- sample(0:1, 1)
		if(r==1){
			ex1 <- paste("lm(",i, "~")
			start = 0
	 		for (j in z){
	 			if(i != j){
	 				r <- sample(0:1, 1)
 					if(r == 1){
 						if(start != 0){
				 			ex1 <- paste(ex1, "+", j)
				 		}
				 		else{
				 			ex1 <- paste(ex1, j)
				 			start = start+1
				 		}
	 				}
 				}
 			}
  			ex1 <- paste(ex1, ")")
 			maxex <- paste(maxex, ex1, ",")
 		}
	}
	sumAIC = 0
	maxex <- substr(maxex, 1, nchar(maxex)-1)
	maxex <- paste(maxex,"))")
	try(sumAIC <- eval(parse(text=maxex)), silent=TRUE)
	sAIC[1] = sumAIC
	sAIC[2] = maxex
	sAIC
}

iteratesumAIC <- function(x, z){
	t <- double(x)
	ms <- character(x)
	for(i in 0:x){
		v <- sumAIC(z)
		like <- as.double(v[1])
		model <- v[2]
		t[i] <- like
		ms[i] <- model
	}
	print(min(t))
	print(ms[which.min(t)])
	t
}

vec <- rep(1, length(z))
mat <- ones(length(z))
for(i in 1:length(z)){
	for(j in 1:length(z)){
	if(i == j){mat[i,j]=0}
	}
}
slice <- z[mat[i,]==1]
itmodel <- function(p){
	ex1 <- paste(p[1])
	for(i in 2:length(p)){
		ex1 <- paste(ex1, "+" , p[i])
	}
}