\name{plotSwarm}
\alias{plotSwarm}

\title{
Intern function for plotting during the Pswarm annealing process
}
\description{
Intern function, generates a scatter plot of the progess of the Pswarm algorithm after every nash equlibirum.
Every point symbolizes a Databot. If a prior classification is given (\code{Cls}) then the Databots have the colors defined by the class labels. 
}
\usage{
plotSwarm(Points,Cls,xlab,ylab,main)
}

\arguments{
  \item{Points}{ProjectedPoints or DataBot positions in cartesian coordinates}
   \item{Cls}{optional, Classification as a numeric vector, if given}
    \item{xlab}{='X', optional, string}
     \item{ylab}{='Y', optional, string}
     \item{main}{="DataBots", optional, string}
}

\author{
Michael Thrun
}

\seealso{
\code{\link{Pswarm}} with \code{PlotIt}=TRUE
}

\keyword{Pswarm}
\keyword{swarm}
