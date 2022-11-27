\name{setdiffMatrix}
\alias{setdiffMatrix}
\title{setdiffMatrix shortens Matrix2Curt by those rows that are in both
matrices.}
\arguments{
\item{Matrix2Curt}{[n,k] matrix, which will be shortened by x rows}
\item{Matrix2compare}{[m,k] matrix whose rows will be compared to those of
Matrix2Curt x rows in Matrix2compare equal rows of Matrix2Curt (order of rows is
irrelevant). Has the same number of columns as Matrix2Curt.}
}
\value{V$CurtedMatrix[n-x,k] Shortened Matrix2Curt}
\description{
setdiffMatrix shortens Matrix2Curt by those rows that are in both matrices.}
\author{
CL,MT 12/2014}
