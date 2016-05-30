library(limma)
library(marray)
v <- list.files(pattern="^GSM")
for(i in 1:length(v)){
	file = paste(v[i],".jpg",sep="")
	ab <- ReadAffy(filenames=v[i])
	jpeg(filename=file)
	image(ab, transfo=log)
	dev.off()
	print(i)
	remove(ab)
}