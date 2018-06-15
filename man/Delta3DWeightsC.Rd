\name{Delta3DWeightsC}
\alias{Delta3DWeightsC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intern function
}
\description{
Implementation of the main equation for SOM, ESOM or the sESOM algorithms}
\usage{
Delta3DWeightsC(vx,Datasample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vx}{ array of weights [1:Lines,1:Columns,1:Weights]}
    \item{Datasample}{NumericVector of one Datapoint[1:n]}
}
\details{
intern function in case of \code{ComputeInR==FALSE} in \code{\link{GeneratePswarmVisualization}},
see chapter 5.3 of [Thrun, 2018] for generalized Umatrix and especially the \code{sESOM4BMUs} algorithm.
}
\value{
modified array of weights [1:Lines,1:Columns,1:]
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}
\author{
Michael Thrun
}