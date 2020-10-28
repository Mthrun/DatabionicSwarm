\name{Pswarm}
\alias{pswarmCpp}
\alias{Pswarm}
\alias{pswarm}

\title{
A Swarm of Databots based on polar coordinates (Polar Swarm).
}
\description{
This projetion method is a part of the databionic swarm which uses the nash equlibrium.
Using polar coordinates for agents (here Databots) in two dimensions has many advantages, for further details see [Thrun, 2018] and [Thrun/Ultsch, 2020].
}
\usage{
Pswarm(DataOrDistance,PlotIt=F,Cls=NULL,Silent=T,

Debug=FALSE,LC=c(NULL,NULL),method= "euclidean",...)
}
\arguments{
  \item{DataOrDistance}{
	matrix, DataOrDistance[1:n,1:n] symmetric matrix of dissimilarities, if variable unsymmetric 
	DataOrDistance[1:d,1:n] is assumed as a dataset and the euclidean distances are calculated of d variables and n cases
}

  \item{PlotIt}{
	Optional, bool, default=FALSE, If =TRUE, Plots the projection during the computation prozess after every nash equlibirum
}
  \item{Cls}{
	Optional, numeric vector [1:n], given Classification in numbers, only for plotting if PlotIt=TRUE, irrelevant for computations
}


  \item{Silent}{
		Optional, bool, default=FALSE, If =TRUE results in various console messages
}
  \item{Debug}{
		Optional, Debug, default=FALSE, =TRUE results in various console messages, depricated for CRAN, because cout is not allowed.
}
  \item{LC}{
		Optional, grid size c(Lines, Columns),

		sometimes it is better to call \code{\link{setGridSize}} separately.
}
  \item{method}{Optional, one of 39 distance methods of \code{parDist} of package parallelDist,  if Data matrix is chosen above}
    \item{\dots}{Further arguments passed on to the \code{parDist} function, e.g. user-defined distance functions}
}

\details{
DBS is a flexible and robust clustering framework that consists
 of three independent modules. The first module is the parameter-free
 projection method Pswarm \code{\link{Pswarm}}, which exploits the concepts of self-organization
 and emergence, game theory, swarm intelligence and symmetry considerations. 
 The second module is a parameter-free high-dimensional data visualization technique, 
 which generates projected points on a topographic map with hypsometric colors \code{\link{GeneratePswarmVisualization}}, 
 called the generalized U-matrix. The third module is a clustering method with no 
 sensitive parameters \code{\link{DBSclustering}}. The clustering can be verified by the visualization and vice versa. 
 The term DBS refers to the method as a whole. 
}
\value{
List with
  \item{ProjectedPoints}{
		[1:n,1:2] xy cartesian coordinates of projection
}
  \item{LC}{
		number of Lines and Columns in c(Lines,Columns)
}
  \item{Control}{
		List, only for intern debugging
}
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Thrun/Ultsch, 2020]  Thrun, M. C., and Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Artificial Intelligence, in press, \url{https://doi.org/10.1016/j.artint.2020.103237}, 2020.
}
\note{
LC is now automatically estimated; LC is the size of the grid c(Lines,Columns), number of Lines and Columns, default c(NULL,NULL) and automatic calculation by \code{\link{setGridSize}}
}
\author{
Michael Thrun
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
Cls=Lsun3D$Cls
InputDistances=as.matrix(dist(Data))
#If not called separately setGridSize() is called in Pswarm
LC=setGridSize(InputDistances)
\donttest{res=Pswarm(InputDistances,LC=LC,Cls=Cls,PlotIt=TRUE)}

\dontshow{
data=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data))
res=Pswarm(distance,LC = c(10,12))
}
}
\keyword{DBS}
\keyword{swarm}
\keyword{swarms}
\concept{Databionic swarm}
\concept{Polar Swarm}
\keyword{pswarm}
\keyword{pSwarm}
\keyword{PSwarm}
\concept{swarm intelligence}
\concept{game theory}
\concept{cluster analysis}
\keyword{projection}
\keyword{projection method}
\concept{dimensionality reduction}
\keyword{visualization}
\keyword{DR}