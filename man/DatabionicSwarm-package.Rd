\name{DatabionicSwarm-package}
\alias{DatabionicSwarm-package}
\alias{DatabionicSwarm}
\docType{package}
\title{
\packageTitle{DatabionicSwarm}
}
\description{
\packageDescription{DatabionicSwarm}
}
\details{
For a brief introduction to \pkg{DatabionicSwarm} please see the vignette \href{../doc/DatabionicSwarm.html}{Short Intro to the Databionic Swarm (DBS)}.

\tabular{ll}{
Package: \tab Databonic swarm\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2018-06-11\cr
License: \tab CC BY-NC-SA 4.0\cr
}

\packageIndices{DatabionicSwarm}

}
\author{
Michal Thrun

Maintainer: \packageMaintainer{DatabionicSwarm}
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Ultsch/Thrun, 2017]  Ultsch, A., & Thrun, M. C.: Credible Visualizations for Planar Projections, in Cottrell, M. (Ed.), 12th International Workshop on Self-Organizing Maps and Learning Vector Quantization, Clustering and Data Visualization (WSOM), IEEE Xplore, France, 2017.

[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Loetsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.



Successfully applied in

[Thrun, 2018]  Thrun, M. C. : Cluster Analysis of the World Gross-Domestic Product Based on Emergent Self-Organization of a Swarm, in Papiez, M. & Smiech, S. (eds.), Proc. 12th Professor Aleksander Zelias International Conference on Modelling and Forecasting of Socio-Economic Phenomena, pp. 523-532, Foundation of the Cracow University of Economics, Cracow, Poland, 2018. 

[Thrun et al., 2018]   Thrun, M. C., Breuer, L., & Ultsch, A. : Knowledge discovery from low-frequency stream nitrate concentrations: hydrology and biology contributions, Proc. European Conference on Data Analysis (ECDA), pp. accepted, Paderborn, Germany, 2018. 

[Weyer-Menkhoff et al., 2018]  Weyer-Menkhoff, I., Thrun, M. C., & Loetsch, J.: Machine-learned analysis of quantitative sensory testing responses to noxious cold stimulation in healthy subjects, European Journal of Pain, Vol. 22(5), pp. 862-874, DOI: 10.1002/ejp.1173, 2018. 

[Kringel et al., 2018]  Kringel, D., Geisslinger, G., Resch, E., Oertel, B. G., Thrun, M. C., Heinemann, S., & Loetsch, J. : Machine-learned analysis of the association of next-generation sequencing based human TRPV1 and TRPA1 genotypes with the sensitivity to heat stimuli and topically applied capsaicin, Pain, in press, doi 10.1097/j.pain.0000000000001222, 2018 
}


\examples{
data('Lsun3D')
##2d projection, without instant visualization of steps
#DistanceMatrix hast to be defined by the user.
InputDistances=as.matrix(dist(Lsun3D$Data))
\donttest{
projection=Pswarm(InputDistances)
#2d projection, with instant visualization 
#of steps and DataMatrix (Distance is Euclidean per default)

projection=Pswarm(Lsun3D$Data,Cls=Lsun3D$Cls,PlotIt=T)
#
##Computation of Generalized Umatrix
# If Non Euclidean Distances are used, Please Use \code{SammonsMapping}
# from the ProjectionBasedClustering package with the correct OutputDimension
# to generate a new DataMatrix from the distances (see SheppardDiagram
# or KruskalStress)
visualization=GeneratePswarmVisualization(Data = Lsun3D$Data,

projection$ProjectedPoints,projection$LC)
## Visualizuation of GenerelizedUmatrix, 
# Estimation of the Number of Clusters=Number of valleys
library(GeneralizedUmatrix)#install if not installed
GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches)
## Automatic Clustering
# number of Cluster from dendrogram (PlotIt=TRUE) or visualization 
Cls=DBSclustering(k=3, Lsun3D$Data, 

visualization$Bestmatches, visualization$LC,PlotIt=FALSE)
# Verification, often its better to mark Outliers manually

GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)
}
\dontrun{
# To generate the 3D landscape in the shape of an island 
# from the toroidal topograpic map visualization
# you may cut your island interactivly around high mountain ranges
Imx = ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(visualization$Umatrix,
visualization$Bestmatches,Cls)

GeneralizedUmatrix::plotTopographicMap(visualization$Umatrix,
visualization$Bestmatches, Cls=Cls,Imx = Imx)
}
\dontrun{
library(ProjectionBasedClustering)#install if not installed
Cls2=ProjectionBasedClustering::interactiveClustering(visualization$Umatrix, 
visualization$Bestmatches, Cls)
}
\dontshow{
data2=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data2))
res=Pswarm(distance,LC = c(10,12))
Cls=DBSclustering(k=2, data2, 

res$ProjectedPoints, res$LC,PlotIt=FALSE)

resUmatrix=GeneratePswarmVisualization(data2,res$ProjectedPoints,res$LC)
}

}
\note{
For interactive Island Generation of a generalized Umatrix  
see \code{interactiveGeneralizedUmatrixIsland} function in the package \pkg{ProjectionBasedClustering}.

If you want to verifiy your clustering result externally, you can use \code{Heatmap} or \code{SilhouettePlot} of the CRAN package \pkg{DataVisualizations}.
}
\keyword{DBS}
\keyword{DataBionic}
\keyword{swarm}
\keyword{swarm intelligence}
\keyword{swarm intelligence}
\keyword{Databionic}
\keyword{cluster analysis}
\keyword{clustering}
\keyword{visualization}
\keyword{DR}
\keyword{Dimensionality Reduction}
\keyword{projection}
\keyword{Pswarm}
\keyword{self-organization}
\keyword{emergence}
\keyword{nash}
\keyword{equilibrium}
\keyword{game theory}
\keyword{projection}
\keyword{projection method}
\keyword{ package }