\name{PswarmCurrentRadiusC2botsPositive}
\alias{PswarmCurrentRadiusC2botsPositive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
intern function, do not use yourself
}
\description{
Finds the weak Nash equilibirium for DataBots in one epoch(Radius), requires the setting of constants, grid, and so on in \code{\link{Pswarm}}
}
\usage{
PswarmCurrentRadiusC2botsPositive( AllDataBotsPosOld,
 Radius, DataDists,
 IndPossibleDBPosR,
 RadiusPositionsschablone,  pp,
 Nullpunkt, Lines,  Columns,
 nBots,  limit, steigungsverlaufind,  StressConstAditiv,debug)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{AllDataBotsPosOld}{ComplexVector [1:n,1], DataBots position in the last Nash-Equlibriuum}
\item{Radius}{double, Radius of payoff function, neighborhood, where other DatsBots can be smelled}
\item{DataDists}{NumericMatrix, Inputdistances[1:n,1:n]}
\item{IndPossibleDBPosR}{ComplexVector, see output of \code{\link{findPossiblePositionsCsingle}}}
\item{RadiusPositionsschablone}{NumericMatrix, see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
\item{pp}{NumericVector, number of jumping simultaneously DataBots of one eppoch (per nash-equilibirum), this vector is linearly monotonically decreasing}
\item{Nullpunkt}{NumericVector, equals \code{which(AllallowedDBPosR0==0,arr.ind=T)}, see see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
\item{Lines}{double, small edge length of rectangulare grid}
\item{Columns}{double, big edge length of rectangulare grid}
\item{nBots}{double, intern constant, equals \code{round(pp[Radius]*DBAnzahl)}}
\item{limit}{int, intern constant, equals \code{ceiling(1/pp[Radius])}}
\item{steigungsverlaufind}{int, intern constant}
\item{StressConstAditiv}{double, intern constant, sum of payoff of all databots in random condition before the algorithm starts}
\item{debug}{optional, bool: If TRUE prints status every 100 iterations}
}

\value{
list of
\item{AllDataBotsPos}{ComplexVector, indizes of DataBot Positions after a weak Nash equlibrium is found}
\item{stressverlauf}{NumericVector, intern result, for debugging only}
\item{fokussiertlaufind}{NumericVector, intern result, for debugging only}
}
\details{

Algorithm is described in [Thrun, 2018, p. 95, Listing 8.1].
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}
\author{
Michael Thrun
}
\keyword{swarm}
\keyword{nash}
\keyword{polar}