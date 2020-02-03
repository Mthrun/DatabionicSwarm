\name{DBSclustering}
\alias{DBSclustering}
\title{
Databonic swarm clustering (DBS)
}
\description{
Automated Clustering approach of the Databonic swarm with abstact U distances, which
are the geodesic distances based on high-dimensional distances combined with low dimensional graph paths by using \code{ShortestGraphPathsC}, see [Thrun/Ultsch, 2020].
}
\usage{
DBSclustering(k, DataOrDistance, BestMatches, LC, StructureType = TRUE, PlotIt = FALSE,
                 method = "euclidean",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{number of clusters, how many to you see in the topographic map (3D landscape)?}
  \item{DataOrDistance}{ Either [1:n,1:d] Matrix of Data (n cases, d dimensions) that will be used. One DataPoint per row
  
  or symmetric Distance matrix [1:n,1:n]}
  \item{BestMatches}{[1:n,1:2] Matrix with positions of Bestmatches=ProjectedPoints, one matrix line per data point}
  \item{LC}{grid size c(Lines,Columns)}
  \item{StructureType}{Optional, bool; =TRUE: compact structure of clusters assumed, =FALSE: connected structure of clusters assumed. For the two options vor Clusters, see [Thrun, 2018] or Handl et al. 2006}
  \item{PlotIt}{Optional, bool, Plots Dendrogramm}
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
 sensitive parameters \code{\link{DBSclustering}} (see [Thrun, 2018, p. 104 ff]). The clustering can be verified by the visualization and vice versa. 
 The term DBS refers to the method as a whole. 
}
\value{
Cls                 [1:n] vector with selected classes of the bestmatches. You can use \code{plotTopographicMap(Umatrix,Bestmatches,Cls)} for verification.
}
\references{
[Thrun/Ultsch, 2020]  Thrun, M. C., and Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Artificial Intelligence, in press, \url{https://doi.org/10.1016/j.artint.2020.103237}, 2020.
}
\author{
Michael Thrun
}
\note{
If you used pswarm with distance matrix instead of a data matrix you may transform your distances into data by using \code{MDS} of the ProjectionBasedClustering package. The correct dimension can be found through the Sheppard diagram or kruskals stress.

Often it is better to mark the outliers  manually after the prozess of clustering and sometimes a clustering can be improved through human interaction [Thrun/Ultsch,2017] <DOI:10.13140/RG.2.2.13124.53124>; use in this case the visualization \code{plotTopographicMap} of the package GeneralizedUmatrix. If you would like to mark the outliers interactivly in the visualization use the Umatrix package in \url{https://www.uni-marburg.de/fb12/datenbionik/software-en}
or the function \code{interactiveClustering()} of the ProjectionBasedClustering package on CRAN.

 If you want to verifiy your clustering result externally, you can use \code{Heatmap} or \code{SilhouettePlot} of the CRAN package \code{DataVisualizations}.
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
InputDistances=as.matrix(dist(Data))
\donttest{
projection=Pswarm(InputDistances)
#autmatic Clustering without GeneralizedUmatrix visualization
Cls=DBSclustering(k=3, Data, 

projection$ProjectedPoints, projection$LC,PlotIt=TRUE)
}
\dontrun{
visualization=GeneratePswarmVisualization(Data,

projection$ProjectedPoints,projection$LC)
## Sometimes an automatic Clustering can be improved 
## thorugh an interactive approach, 
## e.g. if Outliers exist (see [Thrun/Ultsch, 2017])
library(ProjectionBasedClustering)
Cls2=ProjectionBasedClustering::interactiveClustering(visualization$Umatrix, 
visualization$Bestmatches, Cls)
}
\dontshow{
data2=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data2))
res2=Pswarm(distance,LC = c(10,12))
Cls=DBSclustering(k=2, data2, 

res2$ProjectedPoints, res2$LC,PlotIt=FALSE)
}

}
\keyword{swarm}
\keyword{cluster analysis}
\keyword{Databonic swarm}
\keyword{clustering}% __ONLY ONE__ keyword per line
\keyword{cluster analysis}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line