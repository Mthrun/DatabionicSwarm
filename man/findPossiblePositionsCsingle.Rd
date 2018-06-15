\name{findPossiblePositionsCsingle}
\alias{findPossiblePositionsCsingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intern function, do not use yourself
}
\description{
Finds all possible jumping position regarding a grid anda Radius for DataBots
}
\usage{
findPossiblePositionsCsingle(RadiusPositionsschablone,
 jumplength, alpha, Lines)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RadiusPositionsschablone}{NumericMatrix, see \code{\link{setPolarGrid}}}
   \item{jumplength}{double radius of databots regarding neighborhood, they can jump to}
    \item{alpha}{double, zu streichen}
     \item{Lines}{double, jumpinglength has to smaller than Lines/2 and Lines/2 has to yield to a integer number. }
}

\value{
 \item{OpenPositions}{NumericMatrix, indizes of open positions}

}
\details{

Algorithm is described in [Thrun, 2018, p. 95, Listing 8.1].
}
\author{
Michael Thrun
}

\seealso{
\code{\link{setPolarGrid}}
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}
\keyword{positions}
\keyword{swarms}
