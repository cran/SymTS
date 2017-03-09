\name{qCTS}
\alias{qCTS}
\title{Quantile Function of CTS Distribution}
\description{
Evaluates the quantile function for the symmetric classical tempered stable distribution. When alpha=0 this is the symmetric variance gamma distribution.}
\usage{
qCTS(x, alpha, c = 1, ell = 1, mu = 0)
}
\arguments{
  \item{x}{
Vector of quantiles.}
  \item{alpha}{
Number in [0,2)}
  \item{c}{
Parameter c >0}
  \item{ell}{
Parameter ell>0}
  \item{mu}{
Location parameter, any real number}
}
\details{
For details about this distribution see the the describtion of dCTS.
}
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
x = (1:9)/10
qCTS(x,.5)
}