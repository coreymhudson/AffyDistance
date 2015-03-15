.packageName <- 'AffyDistance'

getDistance <-
function(condition1, condition2){
	v <- cbind(condition1, condition2)
	nCondition1 <- ncol(condition1)
	nCondition2 <- ncol(condition2)
	distances <- apply(v, 1, function(v){as.numeric(ks.test(v[1:nCondition1], v[nCondition1+1:nCondition2])$statistic)})
	return(distances)
}

