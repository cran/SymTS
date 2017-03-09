\name{dPowTS}
\alias{dPowTS}
\title{PDF of PowTS Distribution}
\description{
Evaluates the pdf for the symmetric power tempered stable distribution.}
\usage{
dPowTS(x, alpha, c = 1, ell = 1, mu = 0)
}
\arguments{
  \item{x}{
Vector of points}
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
The integration is preformed using the QAWF method in the GSL library for C. For this distribution the Rosinski measure R(dx) = c*(alpha+ell+1)*(alpha+ell)*(1+|x|)^(-2-alpha-ell)(dx). }
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.}
\author{Michael Grabchak and Lijuan Cao}
\note{
We do not allow for the case alpha=0 and c<=.5*(1+ell), as, in this case, the pdf is unbounded. This does not affect pPowTS, qPowTS, or rPowTS.}
\examples{
x = (-10:10)/10
dPowTS(x,.5)
}
