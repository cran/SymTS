\name{rPowTS}
\alias{rPowTS}
\title{Simulation from PowTS Distribution}
\description{Simulates from the symmetric power tempered stable distribution. The simulation is performed by numerically evaluating the quantile function.
}
\usage{
rPowTS(r, alpha, c = 1, ell = 1, mu = 0)
}
\arguments{
  \item{r}{
Number of observations.}
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
For this distribution the Rosinski measure R(dx) = c*(alpha+ell+1)*(alpha+ell)*(1+|x|)^(-2-alpha-ell)(dx). }
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
pPowTS(10,.5)
}