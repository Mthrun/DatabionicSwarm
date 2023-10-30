\name{GeneratePswarmVisualization}
\alias{GeneratePswarmVisualization}

\title{
Generates the Umatrix for Pswarm algorithm
}
\description{
DBS is a flexible and robust clustering framework that consists of three
independent modules. The first module is the parameter-free projection method
Pswarm \code{\link{Pswarm}}, which exploits the concepts of self-organization 
and emergence, game theory, swarm intelligence and symmetry considerations. The
second module is a parameter-free high-dimensional data visualization technique,
which generates projected points on a topographic map with hypsometric colors
\code{\link{GeneratePswarmVisualization}}, called the generalized U-matrix. The
third module is a clustering method with no sensitive parameters
\code{\link{DBSclustering}}. The clustering can be verified by the visualization
and vice versa. The term DBS refers to the method as a whole.

The \code{\link{GeneratePswarmVisualization}} function generates the special
case (please see [Thrun, 2018]) of the generalized Umatrix with the help of an
unsupervised neural network (simplified emergent self-organizing map  published
in [Thrun/Ultsch, 2020]). From the generalized Umatrix a topographic map with
hypsometric tints can be visualized. To see this visualization use
\code{\link[GeneralizedUmatrix]{plotTopographicMap}} of the package
\pkg{GeneralizedUmatrix}.
}
\usage{
GeneratePswarmVisualization(Data,ProjectedPoints,LC,PlotIt=FALSE,
ComputeInR=FALSE,Parallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Data}{[1:n,1:d] array of data: n cases in rows, d variables in columns}
\item{ProjectedPoints}{matrix, ProjectedPoints[1:n,1:2]  n by 2 matrix
containing coordinates of the Projection: A matrix of the fitted configuration.
See output of \code{\link{Pswarm}} for further details}
\item{LC}{size of the grid c(Lines,Columns), number of Lines and Columns
automatic calculated by \code{\link{setGridSize}} in \code{\link{Pswarm}}

Sometimes is better to choose a different grid size, e.g. to to reduce
computional effort contrary to SOM, here the grid size defined only the
resolution of the visualizations. The real grid size is predefined by Pswarm,
but you may choose a factor x*res$LC if you so desire. Therefore, The resulting
grid size is given back in the Output.}
\item{PlotIt}{Optional, default(FALSE), If TRUE than uses
\code{plotTopographicMap} of the package GeneralizedUmatrix is plotted as a
topview in the tiled option, see details for explanation.}
\item{ComputeInR}{Optional, =TRUE: Rcode, =FALSE C++ implementation}

\item{Parallel}{Optional, =TRUE: Parallel C++ implementation, =FALSE C++ implementation}

}
\details{
Tiled: The topographic map is visualized 4 times because the projection is
toroidal. The reason is that there are no border in the visualizations and
clusters (if they exist) are not disrupted by borders of the plot.

If you used \code{\link{Pswarm}} with distance matrix instead of a data matrix
(in the sense that you do not have any data matrix available), you may transform
your distances into data by using \code{MDS} of the
\pkg{ProjectionBasedClustering} package in order to use the
\code{\link{GeneratePswarmVisualization}} function. The correct dimension can be
found through the Sheppard diagram or kruskals stress.
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through
Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer,
Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 

[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A.: Uncovering High-Dimensional
Structures of Projections from Dimensionality Reduction Methods, MethodsX,
Vol. 7, pp. 101093, \doi{10.1016/j.mex.2020.101093}, 2020. 
}
\note{
If you used pswarm with distance matrix instead of a data matrix you can mds
transform your distances into data (see the \code{MDS} function of the
ProjectionBasedClustering package.). The correct dimension can be found through
the Sheppard diagram or kruskals stress.
}
\value{
list of
\item{Bestmatches}{Numeric matrix [1:n,1:2], BestMatches of the Umatrix,
contrary to ESOM they are always fixed, because predefined by GridPoints.}
\item{Umatrix}{Numeric matrix [1:Lines,1:Columns],}
\item{WeightsOfNeurons}{Numeric 3D array [1:Lines,1:Columns,1:d], d is the
dimension of the weights, the same as in the ESOM algorithm}
\item{GridPoints}{Integer matrix [1:n,1:2], quantized projected points:
projected points now lie on a predefined grid.}
\item{LC}{c(Lines,Columns), normally equal to grid size of Pswarm, sometimes it
a better or a lower resolution for the visualization is better. Therefore here
the grid size of the neurons is given back.
}
\item{PlotlyHandle}{If PlotIt=FALSE: NULL, otherwise plotly object for ploting
topview of topographic map.
}
}
\author{
Michael Thrun
}
\note{
The extraction of an island out of the generalized Umatrix can be performed
using the \code{interactiveGeneralizedUmatrixIsland} function in the package
\pkg{ProjectionBasedClustering}.

The main code of both functions \code{GeneralizedUmatrix} and
\code{GeneratePswarmVisualization} is the same C++ function
\code{sESOM4BMUs} which is described in [Thrun/Ultsch, 2020].
}
\seealso{
\code{\link{Pswarm}} and
\code{\link[GeneralizedUmatrix]{plotTopographicMap}} and
\code{\link[GeneralizedUmatrix]{GeneralizedUmatrix}} of the package
\pkg{GeneralizedUmatrix}
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
Cls=Lsun3D$Cls
InputDistances=as.matrix(dist(Data))
\donttest{
projList=Pswarm(InputDistances)
genUmatrixList=GeneratePswarmVisualization(Data,
  projList$ProjectedPoints,projList$LC,
  Parallel=FALSE)#CRAN guidelines do not allow =TRUE for testing
library(GeneralizedUmatrix)
plotTopographicMap(genUmatrixList$Umatrix,genUmatrixList$Bestmatches,Cls)
}
\dontshow{data2=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data2))
projList2=Pswarm(distance,LC = c(10,12))
genUmatrixList=GeneratePswarmVisualization(matrix(runif(n = 100),10,10),
  projList2$ProjectedPoints,
  projList2$LC,
  Parallel=FALSE)#CRAN guidelines do not allow =TRUE for testing
}

}
\concept{Databionic swarm}
\concept{generalized Umatrix}
\concept{cluster analysis}
\keyword{ESOM}
\keyword{U-matrix}
\keyword{Umatrix}
\keyword{sESOM}
\keyword{SOM}
\keyword{swarm}
\keyword{self-organizing-map}
\keyword{visualization}
\keyword{DR}
\keyword{Umatrix}
\keyword{U-matrix}