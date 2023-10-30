\name{UniquePoints}
\alias{UniquePoints}
\title{Unique Points}
\description{
return only the unique points in Datapoints
}
\usage{UniquePoints(Datapoints, Cls, Eps=1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Datapoints}{[1:n,1:d] numeric matrix of Datapoints points of dimension d, the
points are in the  rows}
\item{Cls}{[1:n] numeric vector of classes for each datapoint.}
\item{Eps}{Optional,scalar above zero that defines minimum non-identical
euclidean distance between two points}
}
\details{
Euclidean distance is computed and used within. Setting \code{Eps} to a very small number results in the identification of unique data points. Setting epsilon to a higher number results in the definition of mesh points within an d-dimensional R-ball graph.
}
\value{
List with
\item{Unique}{[1:k,1:d] Datapoints points  without duplicate points}
\item{UniqueInd}{[1:k] index vector such that Unique ==  Datapoints[UniqueInd,],
it has k non-consecutive numbers or labels, each label defines a row number
within Datapoints[1:n,1:d] of a unique data point}
\item{Uniq2DatapointsInd }{ [1:n] index vector. It has k unique index numbers
representing the arbitrary labels. Each labels is mapped uniquely to a point in
\code{Unique}. Logically in a way such that
Datapoints ==  Unique[Uniq2DatapointsInd,] (will not work directly in R this
way)}
\item{NewUniqueInd}{[1:k] index vector stating the index of the newly defined
datastructure Unique.}
\item{NewUniq2DataIdx}{[1:k] index vector such that
Unique[NewUniq2DataIdx,] ==  Datapoints[Uniq2DatapointsInd,], it has n
non-consecutive numbers or labels, each label defines a row number within
Unique[1:k,1:d] of a unique data point}
\item{IsDuplicate}{[1:n,1:n] matrix,for i!=j IsDuplicate[i,j]== 1
if Datapoints[i,] == Datapoints[j,]    IsDuplicate[i,i]==0}
\item{Eps}{Numeric stating the neighborhood radius around unique points.}
}
\author{
Michael Thrun
}
\examples{
Datapoints  = rbind(c(0,0), c(1,1), c(2,2))
Datapoints2 = rbind(Datapoints, Datapoints+0.001)
Datapoints3 = rbind(Datapoints2, c(1,1)-0.001)

Datapoints  = rbind(c(0,0), c(0,0.015), c(0,0.01), c(0,0.015))

V1 = UniquePoints(Datapoints = Datapoints, Eps = 0.01)
V2 = UniquePoints(Datapoints = Datapoints2, Eps = 0.01)
V3 = UniquePoints(Datapoints = Datapoints3, Eps = 0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{points}
\keyword{unique}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{data points}
% \concept{ ~cpt2 }
% Use only one concept per line.
