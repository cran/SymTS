\name{pSaS}
\alias{pSaS}
\title{CDF of Symmetric Stable Distribution}
\description{
Evaluates the cdf for the symmetric alpha stable distribution. For alpha=1 this is the Cauchy distribution.}
\usage{
pSaS(x, alpha, c = 1, mu = 0)
}
\arguments{
  \item{x}{
Vector of probabilities.}
  \item{alpha}{
Index of stability; Number in (0,2)}
  \item{c}{
Scale parameter, c>0}
  \item{mu}{
Location parameter, any real number}
}
\details{
The integration is preformed using the QAWF method in the GSL library for C. The characteristic function is

f(t) = e^(-c |t|^alpha) *e^(i*t*mu).

}
\references{
G. Samorodnitsky and M. Taqqu (1994). Stable Non-Gaussian Random Processes: Stochastic Models with Infinite Variance. Chapman & Hall, Boca Raton.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
x = (-10:10)/10
pSaS(x,.5)
}