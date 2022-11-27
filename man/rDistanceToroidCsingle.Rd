\name{rDistanceToroidCsingle}
\alias{rDistanceToroidCsingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intern function for \code{\link{Pswarm}}
}
\description{
toroid distance calculation
}
\usage{
rDistanceToroidCsingle(AllDataBotsPosX, AllDataBotsPosY, AllallowedDBPosR0,
Lines, Columns,  Nullpunkt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{AllDataBotsPosX}{NumericVector [1:n,1], positions of on grid}
\item{AllDataBotsPosY}{NumericVector [1:n,1], positions of on grid}
\item{AllallowedDBPosR0}{NumericMatrix}
\item{Lines}{double}
\item{Columns}{double}
\item{Nullpunkt}{NumericVector}
}
\note{do not use yourself}
\value{
numeric matrix of toroid Distances[1:n,1:n]
}

\author{
Michael Thrun
}

\seealso{
\code{\link{Pswarm}}
}

\details{
Part of the algorithm described in [Thrun, 2018, p. 95, Listing 8.1].
}

\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through
Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer,
Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018.}
\keyword{toroid}% use one of  RShowDoc("KEYWORDS")
\keyword{distances}% __ONLY ONE__ keyword per line
