\name{UniquePoints}
\alias{UniquePoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unique Points
}
\description{
return only the unique points in Datapoints
}
\usage{
UniquePoints(Datapoints,Eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Datapoints}{
[1:n,1:d] matrix of Datapoints points of dimension d, the points are in the  rows
}
  \item{Eps}{
Optional,scalar above zero that defines minimum non-identical euclidean distance between two points
}
}
\details{
Euclidean distance is computed.
}
\value{
List with
\item{Unique }{[1:u,1:d] Datapoints points  without duplicate points}
\item{IsDuplicate }{[1:n,1:n] matrix,for i!=j IsDuplicate[i,j]== 1  if Datapoints[i,] == Datapoints[j,]    IsDuplicate[i,i]==0}
\item{UniqueInd }{[1:u] index vector such that Unique ==  Datapoints[UniqueInd,]}
\item{Uniq2DatapointsInd }{[1:n] index vector such that   Datapoints ==  Unique[Uniq2DatapointsInd,]}
}

\author{
Michael Thrun
}

\examples{
Datapoints2D=rbind(c(1,2),c(1,2),c(1,3),c(3,1))
V=UniquePoints(Datapoints2D)
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
