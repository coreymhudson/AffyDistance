\name{transformAffyData}
\alias{transformAffyData}
\title{
Transform affy spot data to probabilities
}
\description{
This function takes an affy array and transforms the data to a probability that a given mRNA is found within a random experiment.
}
\usage{
transformAffyData(affy.list)
}
\arguments{
  \item{affy.list}{
a list of filenames pointing to a set of affymetrix experiments.
}
}
\details{
}
\value{
\item{transformedArray }{A matrix of size (col=experiments, rows=spots, values=probabilities)}
}
\references{
}
\author{
Corey M. Hudson <coreymhudson@gmail.com>
}
\note{
}

\seealso{
}
\examples{
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
}
\keyword{ affy }
