\name{RobustNormalization}
\alias{RobustNormalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RobustNormalization
}
\description{
RobustNormalization as described in [Milligan/Cooper, 1988].
}
\usage{
RobustNormalization(Data,Centered=FALSE,Capped=FALSE,

na.rm=TRUE,WithBackTransformation=FALSE,

pmin=0.01,pmax=0.99) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{[1:n,1:d] data matrix of n cases and d features

}
  \item{Centered}{ centered data around zero by median if TRUE
}
  \item{Capped}{TRUE: outliers are capped above 1 or below -1 and set to 1 or -1.
}
  \item{na.rm}{If TRUE, infinite vlaues are disregarded
}
\item{WithBackTransformation}{
If in the case for forecasting with neural networks a backtransformation is required, this parameter can be set to 'TRUE'.
}

  \item{pmin}{defines outliers on the lower end of scale}
  \item{pmax}{defines outliers on the higher end of scale}
}
\details{
Normalizes features either between -1 to 1 (Centered=TRUE) or 0-1 (Centered=TRUE) without changing the distribution of a feature itself. For a more precise description please read [Thrun, 2018, p.17].

"[The] scaling of the inputs determines the effective scaling of the weights in the last layer of a MLP with BP neural netowrk, it can have a large effect on the quality of the final solution. At the outset it is besto to standardize all inputs to have mean zero and standard deviation 1 [(or at least the range under 1)]. This ensures all inputs are treated equally in the regularization prozess, and allows to choose a meaningful range for the random starting weights."[Friedman et al., 2012]
}
\value{
if \code{WithBackTransformation=FALSE}: TransformedData[1:n,1:d] i.e., normalized data matrix of n cases and d features

if \code{WithBackTransformation=TRUE}: List with
\item{TransformedData}{[1:n,1:d]  normalized data matrix of n cases and d features}
\item{MinX}{[1:d] numerical vector used for manual back-transformation of each feature}
\item{MaxX}{[1:d] numerical vector used for manual back-transformation of each feature}
\item{Denom}{[1:d] numerical vector used for manual back-transformation of each feature}
\item{Center}{[1:d] numerical vector used for manual back-transformation of each feature}
}
\references{
[Milligan/Cooper, 1988]  Milligan, G. W., & Cooper, M. C.: A study of standardization of variables in cluster analysis, Journal of Classification, Vol. 5(2), pp. 181-204. 1988.

[Friedman et al., 2012]  Friedman, J., Hastie, T., & Tibshirani, R.: The Elements of Statistical Learning, (Second ed. Vol. 1), Springer series in statistics New York, NY, USA:, ISBN, 2012.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 

}
\author{
Michael Thrun
}
\examples{
Scaled = RobustNormalization(rnorm(1000, 2, 100), Capped = TRUE)
hist(Scaled)

m = cbind(c(1, 2, 3), c(2, 6, 4))
List = RobustNormalization(m, FALSE, FALSE, FALSE, TRUE)
TransformedData = List$TransformedData

mback = RobustNorm_BackTrafo(TransformedData, List$MinX, List$Denom, List$Center)

sum(m - mback)
}

\seealso{
\code{\link{RobustNorm_BackTrafo}}
}

\keyword{RobustNormalization}% use one of  RShowDoc("KEYWORDS")
