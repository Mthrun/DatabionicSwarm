\name{GeneratePswarmVisualization}
\alias{GeneratePswarmVisualization}

\title{
Generates the Umatrix for Pswarm algorithm
}
\description{
Generates the special case of the generalized Umatrix with the help of an unsupervised neural network (here self-organizing map).
From the generalized Umatrix a topographic map with hypsometric tints can be derived.
To see this visualization use \code{plotTopographicMap} of the package GeneralizedUmatrix.
}
\usage{
GeneratePswarmVisualization(Data,ProjectedPoints,LC,PlotIt=FALSE,ComputeInR=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Data}{[1:n,1:d] array of data: n cases in rows, d variables in columns}
\item{ProjectedPoints}{matrix, ProjectedPoints[1:n,1:2]  n by 2 matrix containing coordinates of the Projection: A matrix of the fitted configuration. see output of \code{\link{Pswarm}} for further details}
  \item{LC}{size of the grid c(Lines,Columns), number of Lines and Columns automatic calculated by \code{\link{setGridSize}} in \code{\link{Pswarm}}
  
Sometimes is better to choose a different grid size, e.g. to to reduce computional effort
contrary to SOM, here the grid size defined only the resolution of the visualizations
The real grid size is predfined by Pswarm, but you may choose a factor x*res$LC if you so desire.
Therefore, The resulting grid size is given back in the Output.
}
\item{PlotIt}{Optional, default(FALSE), If TRUE than uses \code{plotTopographicMap} of the package GeneralizedUmatrix}
\item{ ComputeInR}{Optional, =TRUE: Rcode, =FALSE C++ implementation}
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

see
chapter 8 of [Thrun, 2018] for Pswarm

chapter 5.3 of [Thrun, 2018] for generalized Umatrix and especially the \code{sESOM4BMUs} algorithm.
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}
\note{
If you used pswarm with distance matrix instead of a data matrix you can mds transform your distances into data (see the \code{MDS} function of the ProjectionBasedClustering package.). The correct dimension can be found through the Sheppard diagram or kruskals stress.
}
\value{
list of
\item{Bestmatches}{matrix [1:n,1:2], BestMatches of the Umatrix, contrary to ESOM they are always fixed, because predefined by GridPoints.}
\item{Umatrix}{matrix [1:Lines,1:Columns],}
\item{WeightsOfNeurons}{array [1:Lines,1:Columns,1:d], d is the dimension of the weights, the same as in the ESOM algorithm}
\item{GridPoints}{matrix [1:n,1:2],quantized projected points: projected points now lie on a predefined grid.}
\item{LC}{c(Lines,Columns), normally equal to grid size of Pswarm, sometimes it a better or a lower resolution for the visualization is better.
Therefore here the grid size of the neurons is given back.
}

}
\author{
Michael Thrun
}
\note{

For interactive Island Generation of an generalized Umatrix see \code{interactiveGeneralizedUmatrixIsland} function in the package GeneralizedUmatrix.

The main code of both functions \code{GeneralizedUmatrix} and
 \code{GeneratePswarmVisualization} is the same C++ function
 \code{sESOM4BMUs} which is described in [Thrun, 2018].
}
\seealso{
\code{\link{Pswarm}} and
\code{plotTopographicMap} and \code{GeneralizedUmatrix} of the package GeneralizedUmatrix
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
Cls=Lsun3D$Cls
InputDistances=as.matrix(dist(Data))

\donttest{
res=Pswarm(InputDistances)
resUmatrix=GeneratePswarmVisualization(Data,res$ProjectedPoints,res$LC)
library(GeneralizedUmatrix)
plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,Cls)
}

\dontshow{
data2=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data2))
res2=Pswarm(distance,LC = c(10,12))
resUmatrix=GeneratePswarmVisualization(matrix(runif(n = 100),10,10),res2$ProjectedPoints,res2$LC)
}

}
\keyword{Databonic swarm}
\keyword{generalized Umatrix}
\keyword{cluster analysis}
\keyword{ESOM}
\keyword{U-matrix}
\keyword{Umatrix}
\keyword{sESOM}
\keyword{SOM}
\keyword{swarm}
\keyword{self-organizing-map}
\keyword{visualization}
\keyword{DR}