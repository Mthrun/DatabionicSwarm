\name{ShortestGraphPathsC}
\alias{ShortestGraphPathsC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shortest GraphPaths = geodesic distances
}
\description{
Dijkstra's SSSP (Single source shortest path) algorithm, from all points to all points
}
\usage{
ShortestGraphPathsC(Adj, Cost)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Adj}{
[1:n,1:n]         0/1 adjascency matrix, e.g. from delaunay graph or gabriel graph
}
  \item{Cost}{
[1:n,1:n]       matrix, distances between n points (normally euclidean)
}
}
\details{
Vertices are the points, edges have the costs defined by weights (normally a distance).
The algorithm runs in runs in O(n*E*Log(V)), see also [Jungnickel, 2013, p. 87]. Further details can be foubd in 
[Jungnickel, 2013, p. 83-87] and [Thrun, 2018, p. 12].
}
\value{
ShortestPaths[1:n,1:n]   
vector, shortest paths (geodesic) to all other vertices including the source vertice itself
from al vertices to all vertices, stored as a matrix
}
\references{
[Dijkstra,1959]	Dijkstra, E. W.: A note on two problems in connexion with graphs, Numerische mathematik, Vol. 1(1), pp. 269-271. 1959.

[Jungnickel, 2013]	Jungnickel, D.: Graphs, networks and algorithms, (4th ed ed. Vol. 5), Berlin, Heidelberg, Germany, Springer, ISBN: 978-3-642-32278-5, 2013.

[Thrun/Ultsch, 2017]	Thrun, M.C., Ultsch, A.: Projection based Clustering, Conf. Int. Federation of Classification Societies (IFCS),DOI:10.13140/RG.2.2.13124.53124, Tokyo, 2017.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

}
\author{
Michael Thrun
}
\seealso{
\code{\link{DijkstraSSSP}}
}
\note{
require C++11 standard (set flag in Compiler, if not set automatically)
}
\keyword{ShortestGraphPaths}% use one of  RShowDoc("KEYWORDS")
\keyword{ShortestPaths}% __ONLY ONE__ keyword per line
