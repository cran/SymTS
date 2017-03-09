\name{dCTS}
\alias{dCTS}
\title{PDF of CTS Distribution}
\description{
Evaluates the pdf for the symmetric classical tempered stable distribution. When alpha=0 this is the symmetric variance gamma distribution. 
}
\usage{
dCTS(x, alpha, c = 1, ell = 1, mu = 0)
}
\arguments{
  \item{x}{
Vector of points.}
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
The integration is preformed using the QAWF method in the GSL library for C. For this distribution the Rosinski measure R(dx) = c*delta_ell(dx) + c*delta_(-ell)(dx), where delta is the delta function. The Levy measure is M(dx) = c*ell^(alpha) *e^(-x/ell)*x^(-1-alpha) dx. The characteristic function is, for alpha not equal 0,1:

f(t) = exp( 2*c*gamma(-alpha)*(1+ell^2 t^2)^(alpha/2)*(cos(alpha*atan(ell*t))-1)) *e^(i*t*mu),

for alpha = 1 it is

f(t) = (1+ell^2 t^2)^c*exp(-2*c*ell*t*atan(ell*t)) *e^(i*t*mu),

and for alpha=0 it is

f(t) = (1+t^2 ell^2)^(-c) *e^(i*t*mu).
}
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.}
\author{Michael Grabchak and Lijuan Cao}
\note{
When alpha=0 and c<=.5, the pdf is unbounded. It is infinite at mu and the method returns Inf in that case. This does not affect pCTS, qCTS, or rCTS.}
\examples{
x = (-10:10)/10
dCTS(x,.5)
}

