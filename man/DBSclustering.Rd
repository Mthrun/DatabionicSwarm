\name{DBSclustering}
\alias{DBSclustering}
\title{
Databonic swarm clustering (DBS)
}
\description{
 DBS is a flexible and robust clustering framework that consists
 of three independent modules. The first module is the parameter-free
 projection method Pswarm \code{\link{Pswarm}}, which exploits the concepts of self-organization
 and emergence, game theory, swarm intelligence and symmetry considerations [Thrun/Ultsch, 2021]. 
 The second module is a parameter-free high-dimensional data visualization technique, 
 which generates projected points on a topographic map with hypsometric colors \code{\link{GeneratePswarmVisualization}}, 
 called the generalized U-matrix. The third module is a clustering method with no 
 sensitive parameters \code{\link{DBSclustering}} (see [Thrun, 2018, p. 104 ff]). The clustering can be verified by the visualization and vice versa. 
 The term DBS refers to the method as a whole.
 
The \code{\link{DBSclustering}} function applies the automated Clustering approach of the Databonic swarm using abstract U distances, which
are the geodesic distances based on high-dimensional distances combined with low dimensional graph paths by using \code{ShortestGraphPathsC}.
}
\usage{
DBSclustering(k, DataOrDistance, BestMatches, LC, StructureType = TRUE, PlotIt = FALSE,
                 method = "euclidean",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{number of clusters, how many to you see in the topographic map (3D landscape)?}
  \item{DataOrDistance}{Either [1:n,1:d] Matrix of Data (n cases, d dimensions) that will be used. One DataPoint per row
  
  or symmetric Distance matrix [1:n,1:n]}
  \item{BestMatches}{[1:n,1:2] Matrix with positions of Bestmatches or ProjectedPoints, one matrix line per data point}
  \item{LC}{grid size c(Lines,Columns), please see details}
  \item{StructureType}{Optional, bool; =TRUE: compact structure of clusters assumed, =FALSE: connected structure of clusters assumed. For the two options vor Clusters, see [Thrun, 2018] or Handl et al. 2006}
  \item{PlotIt}{Optional, bool, Plots Dendrogramm}
  \item{method}{Optional, one of 39 distance methods of \code{parDist} of package parallelDist,  if Data matrix is chosen above}
   \item{\dots}{Further arguments passed on to the \code{parDist} function, e.g. user-defined distance functions}

  }
\details{
The input of the \code{LC} parameter depends on the choice of \code{Bestmatches} input argument. Usually as the name of the argument states, the Bestmatches of the \code{\link{GeneratePswarmVisualization}} function are used which is define in the notation of self-organizing map. In this case please see example one. However, as written above, clustering and visualization can be applied independently of each other. In this case the places of Lines L and Columns C are switched. Hence, one should give  \code{\link{DBSclustering}} the argument \code{LC[2,1]} as shown in example 2. 

Often it is better to mark the outliers  manually after the prozess of clustering and sometimes a clustering can be improved through human interaction [Thrun/Ultsch,2017] <DOI:10.13140/RG.2.2.13124.53124>; use in this case the visualization \code{\link[GeneralizedUmatrix]{plotTopographicMap}} of the package GeneralizedUmatrix. If you would like to mark the outliers interactivly in the visualization use the \pkg{ProjectionBasedClustering} package with the function \code{interactiveClustering()}, or for full interactive clustering \code{IPBC()}. The package is available on CRAN. An example is shown in case of  \code{interactiveClustering()} function in the third example.
}
\value{
[1:n] numerical vector of numbers defining the classification as the main output of this cluster analysis for the n cases of data corresponding to the n bestmatches. It has k unique numbers representing the arbitrary labels of the clustering. You can use \code{plotTopographicMap(Umatrix,Bestmatches,Cls)} for verification.
}
\references{
[Thrun/Ultsch, 2021]  Thrun, M. C., and Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Artificial Intelligence, Vol. 290, pp. 103237, \doi{10.1016/j.artint.2020.103237}, 2021.
}
\author{
Michael Thrun
}
\note{
If you want to verifiy your clustering result externally, you can use \code{Heatmap} or \code{SilhouettePlot} of the package \pkg{DataVisualizations} available on CRAN.
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
InputDistances=as.matrix(dist(Data))
\donttest{projection=Pswarm(InputDistances)}
## Example One
\donttest{genUmatrixList=GeneratePswarmVisualization(Data,

projection$ProjectedPoints,projection$LC)

Cls=DBSclustering(k=3, Data, 

genUmatrixList$ProjectedPoints, genUmatrixList$LC,PlotIt=TRUE)
}

## Example Two
\donttest{#automatic Clustering without GeneralizedUmatrix visualization
Cls=DBSclustering(k=3, Data, 

projection$ProjectedPoints, projection$LC[c(2,1)],PlotIt=TRUE)
}
\dontrun{
## Example Three
## Sometimes an automatic Clustering can be improved 
## thorugh an interactive approach, 
## e.g. if Outliers exist (see [Thrun/Ultsch, 2017])
library(ProjectionBasedClustering)
Cls2=ProjectionBasedClustering::interactiveClustering(genUmatrixList$Umatrix, 
genUmatrixList$Bestmatches, Cls)
}
\dontshow{
data2=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data2))
res2=Pswarm(distance,LC = c(10,12))
Cls=DBSclustering(k=2, data2, 

res2$ProjectedPoints, res2$LC[c(2,1)],PlotIt=FALSE)
}

}
\keyword{DBS}
\keyword{swarm}
\keyword{clustering}% __ONLY ONE__ keyword per line
\keyword{cluster analysis}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line
\concept{Databonic swarm}
\concept{cluster analysis}