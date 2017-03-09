\name{qPowTS}
\alias{qPowTS}
\title{Quantile Function of PowTS Distribution}
\description{
Evaluates the quantile function for the symmetric power tempered stable distribution.}
\usage{
qPowTS(x, alpha, c = 1, ell = 1, mu = 0)
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
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
x = (1:9)/10
qPowTS(x,.5)
}
