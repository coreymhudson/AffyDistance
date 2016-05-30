source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
cat("Reading Table\n")
data <- read.table("~/Desktop/cancer_R/stem_cell_HiPS.txt", header=T, stringsAsFactors=FALSE)
data.spotid <- data[,1]
expression_vals <- data[sapply(data, function(data) !all(is.character(data)))]
cat("Transforming Data\n")
transformed_exp <- apply(expression_vals, 2, function(expression_vals){ s <- sum(expression_vals); v <- expression_vals/s; v})
transvalue = double(nrow(transformed_exp))
cat("Fitting Values\n")
for(i in 1:nrow(transformed_exp)){
	v <- normexp.fit(transformed_exp[i,])
	if (v$convergence != 0){
		cat(data.spotid[i], "\n")
		transvalue[i] <- NA
	}
	else{
		cat(i,"\n")
		transvalue[i] <- exp(v$par[3])
	}
}
mat <- cbind(data.spotid, transvalue)
mat.spotid <- mat[,1]
mat.lambda <- mat[,2]
cat("Writing table\n")
write.table(mat, file="stem_cell.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


> sig <- sqrt(0.5 * var(SBP1-SBP2))
> sig
[1] 9.148137
> bw0 <- bw.dnrd(SBP2, sig=sig, error="normal")
> temp <- bw0
> ibw <- rep(0, 20)
> for (i in 1:20){
+ temp <- bw.dboot2(SBP2, sig=sig, h0=temp, error="normal", B=1000)
+ ibw[i] <- temp
+ }
> ibw1 <- mean(ibw)
> SBP2.dec <- DeconPdf(SBP2, sig=sig, error="normal", bw=ibw1, fft=TRUE)
> plot(SBP2.dec, lwd=3, main="")
> lines(density(SBP2, adjust=1.6), lty=3, lwd=3, col="blue")
> library(MASS)
> fitdistr(SBP2, "exponential")

gsm <- readGEO
GSEmat <- matrix(data=NA, nrow=length(probesets),  ncol=length(gsmplatforms))
colnames(GSEmat) <- names(GSMList(gse))
for(i in 1:gsmplatforms){
	coli <- as.double(Table(GSMList(gse)[[]])$VAL)
	GSEmat[,i] <- coli
}