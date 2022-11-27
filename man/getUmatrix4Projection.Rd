\name{getUmatrix4Projection}
\alias{getUmatrix4Projection}
\title{depricated! see GeneralizedUmatrix() Generalisierte U-Matrix fuer
Projektionsverfahren }
\usage{
getUmatrix4Projection(Data,ProjectedPoints,
PlotIt=TRUE,Cls=NULL,toroid=T,Tiled=F,ComputeInR=F)
}
\arguments{
\item{Data}{[1:n,1:d] Numeric matrix: n cases in rows, d variables in columns}
\item{ProjectedPoints}{[1:n,2]n by 2 matrix containing coordinates of the
Projection: A matrix of the fitted configuration.}
\item{PlotIt}{Optional,bool, defaut=FALSE, if =TRUE: U-Marix of every current
Position of Databots will be shown}
\item{Cls}{Optional, For plotting, see \code{plotUmatrix} in package Umatrix}
\item{toroid}{Optional, Default=FALSE,
==FALSE planar computation
==TRUE: toroid borderless computation, set so only if projection method is also toroidal
}
\item{Tiled}{Optional,For plotting see \code{plotUmatrix} in package Umatrix}
\item{ComputeInR}{Optional,  =T: Rcode, =F Cpp Code}}
\value{List with
\item{Umatrix}{[1:Lines,1:Columns] (see \code{ReadUMX} in package DataIO)}
\item{EsomNeurons}{[Lines,Columns,weights] 3-dimensional numeric array
(wide format), not wts (long format)}
\item{Bestmatches}{[1:n,OutputDimension] GridConverted Projected Points
information converted by convertProjectionProjectedPoints() to predefined Grid
by Lines and Columns}
\item{gplotres}{Ausgabe von ggplot}
\item{unbesetztePositionen}{Umatrix[unbesetztePositionen] = NA}
}
\description{
depricated! see GeneralizedUmatrix()}
\author{Michael Thrun}

\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, ISBN: 978-3-658-20539-3, Heidelberg, 2018.
}
\examples{
data("Lsun3D")
Data=Lsun3D$Data
Cls=Lsun3D$Cls
InputDistances=as.matrix(dist(Data))
res=cmdscale(d=InputDistances, k = 2, eig = TRUE, add = FALSE, x.ret = FALSE)
ProjectedPoints=as.matrix(res$points)
# Stress = KruskalStress(InputDistances, as.matrix(dist(ProjectedPoints)))
#resUmatrix=GeneralizedUmatrix(Data,ProjectedPoints)
#plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,Cls)
}