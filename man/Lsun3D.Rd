\name{Lsun3D}
\alias{Lsun3D}

\docType{data}
\title{
Lsun3D is part of the Fundamental Clustering Problem Suit (FCPS) [Thrun/Ultsch, 2020].
}
\description{
clearly defined clusters, different variances
}
\usage{data("Lsun3D")}

\details{
Size 404, Dimensions 3

Dataset defined discontinuites, where the clusters have different variances.
Three main Clusters, and four Outliers (in Cluster 4). See for a more detailed description in [Thrun, 2018].
}

\references{
[Thrun/Ultsch, 2020]	Thrun, M. C., & Ultsch, A.: Clustering Benchmark Datasets Exploiting the Fundamental Clustering Problems, Data in Brief,Vol. 30(C), pp. 105501, DOI 10.1016/j.dib.2020.105501 , 2020. 
}
\examples{
data(Lsun3D)
str(Lsun3D)
Cls=Lsun3D$Cls
Data=Lsun3D$Data
}
\keyword{datasets}
\keyword{Lsun3D}
\keyword{FCPS}
