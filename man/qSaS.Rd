\name{qSaS}
\alias{qSaS}
\title{Quantile Function of Symmetric Stable Distribution}
\description{
Evaluates the quantile function for the symmetric alpha stable distribution. For alpha=1 this is the Cauchy distribution.}
\usage{
qSaS(x, alpha, c = 1, mu = 0)
}
\arguments{
  \item{x}{
Vector of points.}
  \item{alpha}{
Index of stability; Number in (0,2)}
  \item{c}{
Scale parameter, c>0}
  \item{mu}{
Location parameter, any real number}
}
\details{
The characteristic function is

f(t) = e^(-c |t|^alpha) *e^(i*t*mu).

}
\references{
G. Samorodnitsky and M. Taqqu (1994). Stable Non-Gaussian Random Processes: Stochastic Models with Infinite Variance. Chapman & Hall, Boca Raton.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
x = (1:9)/10
qSaS(x,.5)
}