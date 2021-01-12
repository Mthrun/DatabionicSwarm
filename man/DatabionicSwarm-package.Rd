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
Version: \tab 1.1.5\cr
Date: \tab 2020-10-28\cr
License: \tab CC BY-NC-SA 4.0\cr
}

\packageIndices{DatabionicSwarm}

}
\author{
Michal Thrun

Maintainer: \packageMaintainer{DatabionicSwarm}
}
\references{
[Thrun/Ultsch, 2021]  Thrun, M. C., and Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Artificial Intelligence, Vol. 290, pp. 103237, \doi{10.1016/j.artint.2020.103237}, 2021.

[Thrun/Ultsch, 2021] Thrun, M. C., & Ultsch, A.: Swarm Intelligence for Self-Organized Clustering (Extended Abstract), in Bessiere, C. (Ed.), 29th International Joint Conference on Artificial Intelligence (IJCAI), Vol. IJCAI-20, pp. 5125--5129, \doi{10.24963/ijcai.2020/720}, Yokohama, Japan, Jan., 2021. 

[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A.: Uncovering High-Dimensional Structures of Projections from Dimensionality Reduction Methods, MethodsX, Vol. 7, pp. 101093, DOI \doi{10.1016/j.mex.2020.101093}, 2020. 

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 

[Ultsch/Thrun, 2017]  Ultsch, A., & Thrun, M. C.: Credible Visualizations for Planar Projections, in Cottrell, M. (Ed.), 12th International Workshop on Self-Organizing Maps and Learning Vector Quantization, Clustering and Data Visualization (WSOM), IEEE Xplore, France, 2017.

[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Loetsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.


Successfully applied in

[Thrun et al., 2018]   Thrun, M. C., Breuer, L., & Ultsch, A. : Knowledge discovery from low-frequency stream nitrate concentrations: hydrology and biology contributions, Proc. European Conference on Data Analysis (ECDA), pp. 46-47, Paderborn, Germany, 2018. 

[Weyer-Menkhoff et al., 2018]  Weyer-Menkhoff, I., Thrun, M. C., & Loetsch, J.: Machine-learned analysis of quantitative sensory testing responses to noxious cold stimulation in healthy subjects, European Journal of Pain, Vol. 22(5), pp. 862-874, DOI \doi{10.1002/ejp.1173}, 2018. 

[Kringel et al., 2018]  Kringel, D., Geisslinger, G., Resch, E., Oertel, B. G., Thrun, M. C., Heinemann, S., & Loetsch, J. : Machine-learned analysis of the association of next-generation sequencing based human TRPV1 and TRPA1 genotypes with the sensitivity to heat stimuli and topically applied capsaicin, Pain, Vol. 159 (7 ), pp. 1366-1381, DOI \doi{10.1097/j.pain.0000000000001222}, 2018 

[Thrun, 2019]  Thrun, M. C.: : Cluster Analysis of Per Capita Gross Domestic Products, Entrepreneurial Business and Economics Review (EBER), Vol. 7(1), pp. 217-231, DOI: \doi{10.15678/EBER.2019.070113}, 2019.

[López-García et al., 2020]  López-García, P., Argote, D. L., & Thrun, M. C.: Projection-based Classification of Chemical Groups and Provenance Analysis of Archaeological Materials, IEEE Access, Vol. 8, pp. 152439-152451, DOI \doi{10.1109/ACCESS.2020.3016244}, 2020. 
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
# If Non Euclidean Distances are used, Please Use \code{MDS}
# from the ProjectionBasedClustering package with the correct OutputDimension
# to generate a new DataMatrix from the distances (see SheppardDiagram
# or KruskalStress)
genUmatrixList=GeneratePswarmgenUmatrixList(Data = Lsun3D$Data,

projection$ProjectedPoints,projection$LC)
## Visualizuation of GenerelizedUmatrix, 
# Estimation of the Number of Clusters=Number of valleys
library(GeneralizedUmatrix)#install if not installed
GeneralizedUmatrix::plotTopographicMap(genUmatrixList$Umatrix,genUmatrixList$Bestmatches)
## Automatic Clustering
# number of Cluster from dendrogram (PlotIt=TRUE) or visualization 
Cls=DBSclustering(k=3, Lsun3D$Data, 

genUmatrixList$Bestmatches, genUmatrixList$LC,PlotIt=FALSE)
# Verification, often its better to mark Outliers manually

GeneralizedUmatrix::plotTopographicMap(genUmatrixList$Umatrix,genUmatrixList$Bestmatches,Cls)
}
\dontrun{
# To generate the 3D landscape in the shape of an island 
# from the toroidal topograpic map visualization
# you may cut your island interactivly around high mountain ranges
Imx = ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(genUmatrixList$Umatrix,
genUmatrixList$Bestmatches,Cls)

GeneralizedUmatrix::plotTopographicMap(genUmatrixList$Umatrix,
genUmatrixList$Bestmatches, Cls=Cls,Imx = Imx)
}
\dontrun{
library(ProjectionBasedClustering)#install if not installed
Cls2=ProjectionBasedClustering::interactiveClustering(genUmatrixList$Umatrix, 
genUmatrixList$Bestmatches, Cls)
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