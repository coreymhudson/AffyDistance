\name{getDistance}
\alias{getDistance}
\title{
Get the Kolmorogov distance between two paired arrays.
}
\description{
This function allows the user to get the Kolmorogov distance between two equal spots on an array.
}
\usage{
getDistance(condition1, condition2)
}
\arguments{
  \item{condition1}{
A transformed affymetrix array
}
  \item{condition2}{
A transformed affymetrix array
}
}
\details{
This function computes the Kolmorogov distance between two affy sets, transformed into probabilities.
}
\value{
A list of distances
}
\references{
}
\author{
Corey M. Hudson <coreymhudson@gmail.com>
}
\note{
}

\seealso{
\code{\link{transformAffyData}}
}
\examples{
function(condition1, condition2){
	v <- cbind(condition1, condition2)
	nCondition1 <- ncol(condition1)
	nCondition2 <- ncol(condition2)
	distances <- apply(v, 1, function(v){as.numeric(ks.test(v[1:nCondition1], v[nCondition1+1:nCondition2])$statistic)})
	return(distances)
  }
}
\keyword{ affy }
