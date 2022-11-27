\name{RobustNorm_BackTrafo}
\alias{RobustNorm_BackTrafo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforms the Robust Normalization back
}
\description{
Transforms the Robust Normalization back if \code{Capped=FALSE}
}
\usage{
RobustNorm_BackTrafo(TransformedData,

MinX,Denom,Center=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{TransformedData}{[1:n,1:d]  matrix}
\item{MinX}{scalar}
\item{Denom}{scalar}
\item{Center}{scalar}
}
\details{
For details see \code{\link{RobustNormalization}}
}
\value{[1:n,1:d] Data matrix}
\author{
Michael Thrun
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{RobustNormalization}}
}
\examples{
data(Hepta)
Data = Hepta$Data
TransList = RobustNormalization(Data, Centered = TRUE, WithBackTransformation = TRUE)

HeptaData = RobustNorm_BackTrafo(TransList$TransformedData,
                                 TransList$MinX,
                                 TransList$Denom,
                                 TransList$Center)

sum(HeptaData - Data) #<e-15
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RobustNormalization}% use one of  RShowDoc("KEYWORDS")
\keyword{BackTransformation_RobustNormalization}% __ONLY ONE__ keyword per line
\keyword{RobustNorm_BackTrafo}