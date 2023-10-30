\name{RelativeDifference}
\alias{RelativeDifference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relative Difference
}
\description{
Calculates the difference between positive x and y values 
}
\usage{
RelativeDifference(X, Y, epsilon = 10^-10,

na.rm=FALSE,Silent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{either a value or numerical vector of [1:n]}
\item{Y}{either a value or numerical vector of [1:n]}
\item{epsilon}{Optional, If both x and y are approximatly zero the output is
also zero}
\item{na.rm}{Optional, function does not work with non finite values. If these
cases should be automatically removed, set parameter TRUE}

\item{Silent}{Optional, if TRUE message abput values below epsilon is not given back}
}
\details{
Contrary to other approaches in this cases the range of values lies between
[-2,2]. The approach is only valid for positive values of\code{X} and \code{Y}.
The realtive difference \code{R} is defined with \deqn{R=\frac{Y-X}{0.5*(X+Y)}}

Negative value indicate that \code{X}  is higher than \code{Y} and positive
values that \code{X} is lower than \code{Y}.
}
\value{R}
\references{
Ultsch, A.: Is Log Ratio a Good Value for Measuring Return in Stock Investments?
GfKl 2008, pp, 505-511, 2008.
}
\author{
Michael Thrun
}
\note{
It can be combined with the \code{GabrielClassificationError} if a clear baseline is defined.
}
\seealso{
\code{\link[DRquality]{GabrielClassificationError}}
}
\examples{
x=c(1:5)
y=runif(5,min=1,max=10)
RelativeDifference(x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RelativeDifference}% use one of  RShowDoc("KEYWORDS")
\keyword{difference}% __ONLY ONE__ keyword per line
\keyword{relative}% __ONLY ONE__ keyword per line