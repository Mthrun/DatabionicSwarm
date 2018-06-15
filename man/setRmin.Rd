\name{setRmin}
\alias{setRmin}
\title{Intern function: Estimates the minimal radius for the Databot scent}
\arguments{
\item{Lines}{x-value determining the size of the map, i.e. how many open places for DataBots will be available  on the 2-dimensional grid BEWARE: has to be able to be divided by 2}
\item{Columns}{y-value  determining the size of the map, i.e. how many open places for DataBots will be available on the 2-dimensional grid Columns>Lines}
\item{AllallowedDBPosR0}{[1:Lines+1,1:Lines+1]Matrix of radii in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump}
\item{p}{percent of gitterpositions, which should be considered}
}
\value{Rmin Minimum Radius}
\description{
estimates the minimal radius on apolar grid in the automated annealing process of Pswarm, details of how can be read in [Thrun, 2018, p. 97]}
\author{Michael Thrun
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}