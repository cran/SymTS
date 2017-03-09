\name{pPowTS}
\alias{pPowTS}
\title{PDF of PowTS Distribution}
\description{
Evaluates the cdf for the symmetric power tempered stable distribution.}
\usage{
pPowTS(x, alpha, c = 1, ell = 1, mu = 0)
}
\arguments{
  \item{x}{
Vector of probabilities.}
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
\examples{
x = (-10:10)/10
pPowTS(x,.5)
}
