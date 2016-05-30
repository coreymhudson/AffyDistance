.packageName <- 'AffyDistance'

transformAffyData <-
function(affy.list){
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

