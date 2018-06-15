\name{DeepSwarming}
\alias{DeepSwarming}

\title{
Multiple Swarms based on polar coordinates using the nash equlibrium
}
\description{
a part of the databionic swarm, see thesis
}
\usage{
DeepSwarming(InputDistances,NumberOfSwarms=12,

LogicalProcessors=4,PlotIt=F,Cls=NULL,Silent=T,Debug=FALSE)
}
\arguments{
  \item{InputDistances}{
	matrix, InputDistances[1:n,1:n] symmetric matrix of dissimilarities, if variable unsymmetric 
	InputDistances[1:d,1:n] is assumed as a dataset and the euclidean distances are calculated of d variables and n cases
}
  \item{NumberOfSwarms}{
Number of pararell Pswarms
}

  \item{LogicalProcessors}{
Number of LogicalProcessors. How many games should be simoultansly computated.
}

  \item{PlotIt}{
	Optional, bool, default=FALSE, If =TRUE, Plots the projection during the computation prozess after every nash equlibirum
}
  \item{Cls}{
	Optional, numeric vector, given Classification in numbers, only for plotting if PlotIt=TRUE, irrelevant for computations
}


  \item{Silent}{
		Optional, bool, default=FALSE, If =TRUE results in various console messages
}
  \item{Debug}{
		Optional, Debug, default=FALSE, =TRUE results in various console messages
}
}

\details{
see thesis
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
Thrun, M. C.:A System for Projection Based Clustering through Self-Organization and Swarm Intelligence, (Doctoral dissertation), Philipps-Universität Marburg, Marburg, 2016.
}
\author{
Michael Thrun
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
Cls=Lsun3D$Cls
InputDistances=as.matrix(dist(Data))
#res=DeepSwarming(InputDistances,4)
}
\keyword{DBS}
\keyword{swarm}
\keyword{MultipleSwarms}
\keyword{DeepSwarming}

\keyword{swarms}
\keyword{Databonic swarm}
\keyword{Databionic Swarm}
\keyword{pswarm}
\keyword{pSwarm}
\keyword{PSwarm}
\keyword{swarm intelligence}