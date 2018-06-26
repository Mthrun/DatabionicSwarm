\name{DelaunayClassificationError}
\alias{DelaunayClassificationError}
\title{
Delaunay Classification Error (DCE)
}
\description{
DCE searches for the k-nearest neighbors of the first delaunay neighbors weighted by the Euclidean Distances of the Inputspace. DCE evaluates these neighbors in the Output space. A low value indicates a better two-dimensional projection of the high-dimensional Input space.

}
\usage{
DelaunayClassificationError(Data,ProjectedPoints,Cls,LC=NULL)
}
\arguments{
  \item{Data}{
[1:n,1:d]
}
  \item{ProjectedPoints}{
[1:n,1:2]
}
  \item{Cls}{
[1:n,1]
}
  \item{LC}{
Optional, default NULL,  Wenn toroid, muss c(Lines,Columns) angeben werden
}
}
\details{
Delaunay classification error (DCE) makes an unbiased evaluation of distance and densitiybased structure which ma be even non-linear seperable. First, DCE utilizes the information provided by a prior classification to assess projected structures. Second, DCE applies the insights drawn from graph theory. Details are described in [Thrun/Ultsch, 2018]

}
\note{
see also chapter 6 of [Thrun, 2018]
}
\value{
list of 

 DCE                            DelaunayClassificationError 
                                NOTE the rest is just for development purposes
																
 DCEperPoint(1:n)               unnormalized DCE of each point: DCE = mean( DCEperPoint)
 
 nn                             the number of points in a relevant neghborhood: 0.5 * 85percentile(AnzNN)
 
 AnzNN(1:n)                     the number of points with a delaunay graph neighborhood
 
 NNdists(1:n,1:nn)              the distances within the relevant nehborhoot, 0 for inner cluster distances
 
 HD(1:nn)                       HD = HarmonicDecay(nn) i.e weight function for the NNdists: DCEperPoint = HD*NNdists
 

}
\references{
[Thrun/Ultsch, 2018]  Thrun, M. C., & Ultsch, A. : Investigating Quality measurements of projections for the Evaluation of Distance and Density-based Structures of High-Dimensional Data, Proc. European Conference on Data Analysis (ECDA), pp. accepted, Paderborn, Germany, 2018. 
}
\author{
Michael Thrun
}

\examples{
data(Hepta)
\donttest{
InputDistances=as.matrix(dist(Hepta$Data))
projection=Pswarm(InputDistances)
DelaunayClassificationError(Hepta$Data,projection$ProjectedPoints,Hepta$Cls,LC=projection$LC)$DCE
}
\dontshow{
n=nrow(Hepta$Data)
DelaunayClassificationError(Hepta$Data,cbind(runif(n),runif(n)),Hepta$Cls)$DCE
}
}

\keyword{Classification Error}

\keyword{Delaunay Classification Error}


\keyword{DCE}


\keyword{Quality measure}


\keyword{QM}


\keyword{Quality measurement}


\keyword{Evaluation of projection methods}




