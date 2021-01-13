\name{DijkstraSSSP}
\alias{DijkstraSSSP}
\title{
Internal function: Dijkstra SSSP
}
\description{
Dijkstra's SSSP (Single source shortest path) algorithm:

gets the shortest path (geodesic distance) from source vertice(point) to all other vertices(points) defined by 
 the edges of the adjasency matrix
}
\usage{
DijkstraSSSP(Adj, Costs, source)
}
\arguments{
  \item{Adj}{
[1:n,1:n]         0/1 adjascency matrix, e.g. from delaunay graph or gabriel graph
}
  \item{Costs}{
[1:n,1:n]       matrix, distances between n points (normally euclidean)
}
  \item{source}{
int, vertice(point) from which to calculate the geodesic distance to all other points
}
}
\details{
Preallocating space for DataStructures accordingly to the maximum possible number of vertices which is fixed set at 
the number 10001.
This is an internal function of \code{\link{ShortestGraphPathsC}}, no errors or mis-usage is caught here.
}
\value{
ShortestPaths[1:n]   vector, shortest paths (geodesic) to all other vertices including the source vertice itself

}
\references{
uses a changed code which is inspired by Shreyans Sheth 28.05.2015, see 
\url{https://ideone.com/qkmt31}
}
\author{
Michael Thrun
}
\note{
runs in O(E*Log(V))
}
\keyword{Single source shortest path}
\keyword{Dijkstra's SSSP}
\keyword{Dijkstra}
\keyword{SSSP}
