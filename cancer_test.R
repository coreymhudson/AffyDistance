loadAffylibraries <- function(){
	library(limma)
	library(gcrma)
	library(MASS)
	library(gamlss)
	library(stats4)
}

load_libraries()

transformAffyData <- function(affy.list){
	ab <- ReadAffy(filenames=affy.list)
	geneNames <- geneNames(ab[,1])
	eset <- rma(ab)
	remove(ab)
	eIntensity <- exprs(eset)
	remove(eset)
	transformedAffy <- apply(eIntensity, 2, function(eIntensity){ s <- sum(eIntensity); v <- eIntensity/s; v})
	remove(eIntensity)
	return(transformedAffy)
}

transformedEmbryo <- transformAffyData(embryo.list)
transformedNormal <- transformAffyData(normal.list)
transformedCancer <- transformAffyData(cancer.list)

getDistance <- function(condition1, condition2){
	v <- cbind(condition1, condition2)
	nCondition1 <- ncol(condition1)
	nCondition2 <- ncol(condition2)
	distances <- apply(v, 1, function(v){as.numeric(ks.test(v[1:nCondition1], v[nCondition1+1:nCondition2])$statistic)})
	return(distances)
}

cancer_embryo_distance = getDistance(cancer, embryo)
cancer_normal_distance = getDistance(cancer, normal)
embryo_normal_distance = getDistance(normal, embryo)

for(i in 1:length(cancer_embryo_distance)){
	if(which.min(c(cancer_embryo_distance[i],cancer_normal_distance[i],embryo_normal_distance[i]))==1){
		test <- as.numeric(wilcox.test(c(transformedEmbryo[i,], transformedCancer[i,]), transformedNormal[i,])$p.value)			
		if(test < 0.05){
			print(geneNames[i])
		}
	}
}


sampleRow <- function(condition1, condition2, condition3){
	vec <- c(condition1, condition2, condition3)
	index <- seq(from=1, to=length(vec))
	a <- sample(index, length(condition1))
	b <- sample(setdiff(index, a), length(condition2))
	c <- setdiff(index, c(a, b))
	aPrime <- sapply(a, function(a){vec[a]})
	bPrime <- sapply(b, function(b){vec[b]})
	cPrime <- sapply(c, function(c){vec[c]})
	d <- list(group1=aPrime, group2=bPrime, group3=cPrime)
	return(d)
}

significantSpot <- function(conditionTest, conditionTestAlt, conditionNull){
	b = FALSE
	dist1 <- as.numeric(ks.test(conditionTest, conditionTestAlt)$statistic)
	dist2 <- as.numeric(ks.test(conditionTest, conditionNull)$statistic)
	dist3 <- as.numeric(ks.test(conditionTestAlt, conditionNull)$statistic)
	if(which.min(c(dist1, dist2, dist3))==1){
		test <- as.numeric(wilcox.test(c(conditionTest, conditionTestAlt), conditionNull)$p.value)
		if(test<0.05){
			b=TRUE
		}
	}
	return(b)
}

significantSpots <- logical(nrow(transformedEmbryo))
for(i in 1:nrow(transformedEmbryo)){
	significantSpots[i] <- significantSpot(transformedEmbryo[i,], transformedCancer[i,], transformedNormal[i,])
}
length(significantSpots[significantSpots==TRUE])

bootstrapSize <- function(conditionTest, conditionTestAlt, conditionNull, size){
	bootstrapCounts <- numeric(size)
	for(i in 1:size){
		significantSpots <- logical(nrow(conditionTest))
		for(j in 1:nrow(conditionTest)){
			sampledrow <- sampleRow(conditionTest[j,], conditionTestAlt[j,], conditionNull[j,])
			significantSpots[j] <- significantSpot(sampledrow$group1, sampledrow$group2, sampledrow$group3)
		}
		print(i)
		bootstrapCounts[i] <- length(significantSpots[significantSpots==TRUE])
	}
	return(bootstrapCounts)
}
	