\name{DelaunayClassificationError}
\alias{DelaunayClassificationError}
\title{
Delaunay Classification Error (DCE)
}
\description{
DCE searches for the k-nearest neighbors of the first delaunay neighbors
weighted by the Euclidean Distances of the Inputspace. DCE evaluates these
neighbors in the Output space. A low value indicates a better two-dimensional
projection of the high-dimensional Input space.
}
\usage{
DelaunayClassificationError(Data,ProjectedPoints,Cls,LC,Gabriel=FALSE,
PlotIt=FALSE,Plotter = "native", Colors = NULL,LineColor= 'grey',
main = "Name of Projection", mainSize = 24,xlab = "X", ylab = "Y", xlim, ylim,
pch,lwd,Margin=list(t=50,r=0,l=0,b=0))
}
\arguments{
\item{Data}{[1:n,1:d] Numeric matrix with n cases and d variables}
\item{ProjectedPoints}{[1:n,1:2] Numeric matrix with 2D points in cartesian coordinates}
\item{Cls}{[1:n] Numeric vector with class labels}
\item{LC}{Optional, Numeric vector of two values determining grid size of the underlying projection}
\item{Gabriel}{Optional, Boolean: TRUE/FALSE => Gabriel/Delauny graph (Default: FALSE => Delaunay)}
\item{PlotIt}{Optional, Boolean: TRUE/FALSE => Plot/Do not plot (Default: FALSE)}
\item{Plotter}{Optional, Character with plot technique (native or plotly)}
\item{Colors}{Optional, Character vector of class colors for points}
\item{LineColor}{Optional, Character of line color used for edges of graph}
\item{main}{Optional, Character plot title}
\item{mainSize}{Optional, Numeric size of plot title}
\item{xlab}{Optional, Character name of x ax}
\item{ylab}{Optional, Character name of y ax}
\item{xlim}{Optional, Numeric vector with two values defining x ax range}
\item{ylim}{Optional, Numeric vector with two values defining y ax range}
\item{pch}{Optional, Numeric of point size (graphic parameter)}
\item{lwd}{Optional, Numeric of linewidth (graphic parameter)}
\item{Margin}{Optional, Margin of plotly plot}
}
\details{
Delaunay classification error (DCE) makes an unbiased evaluation of distance and
densitiybased structure which ma be even non-linear seperable. First, DCE
utilizes the information provided by a prior classification to assess projected
structures. Second, DCE applies the insights drawn from graph theory. Details
are described in [Thrun/Ultsch, 2018]}
\note{see also chapter 6 of [Thrun, 2018]}
\value{
list of 
\item{DCE}{DelaunayClassificationError 
                            NOTE the rest is just for development purposes}
\item{DCEperPoint}{[1:n] unnormalized DCE of each point: DCE =
mean(DCEperPoint)}
\item{nn}{the number of points in a relevant neghborhood:
0.5 * 85percentile(AnzNN)}
\item{AnzNN}{[1:n] the number of points with a delaunay graph neighborhood}
\item{NNdists}{[1:n,1:nn] the distances within the relevant neighborhood, 0 for
inner cluster distances}
\item{HD}{[1:nn] HD = HarmonicDecay(nn) i.e weight function for the NNdists:
DCEperPoint = HD*NNdists}
}
\references{
[Thrun/Ultsch, 2018]  Thrun, M. C., & Ultsch, A. : Investigating Quality
measurements of projections for the Evaluation of Distance and Density-based
Structures of High-Dimensional Data, Proc. European Conference on Data Analysis
(ECDA), pp. accepted, Paderborn, Germany, 2018. 
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

