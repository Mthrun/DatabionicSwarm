\name{ProjectedPoints2Grid}
\alias{ProjectedPoints2Grid}

\title{
Transforms ProjectedPoints to a grid
}
\description{
quantized xy cartesianncoordinates of ProjectedPoints
}
\usage{
ProjectedPoints2Grid(ProjectedPoints, Lines, Columns,PlotIt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ProjectedPoints}{[1:n,1:2] matrix of cartesian xy coordinates}
  \item{Lines}{double, length of small side of the rectangular grid}
  \item{Columns}{double, length of big side of the rectangular grid}
  \item{PlotIt}{optional, bool, shows the result if TRUE}
}
\details{
intern function, described in [Thrun, 2018, p.47]
}
\value{
BestMatches[1:n,1:3] columns in order: Key,Lines,Columns
}

\author{
Michael Thrun

}
\seealso{\code{\link{GeneratePswarmVisualization}}
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 
}
\keyword{GeneratePswarmVisualization}% use one of  RShowDoc("KEYWORDS")
\keyword{generlized Umatrix}% __ONLY ONE__ keyword per line
