\name{rSaS}
\alias{rSaS}
\title{Simulation from Symmetric Stable Distribution}
\description{
Simulates from the symmetric alpha stable distribution. When alpha=1 this is the Cauchy distribution. The simulation is performed using a well-known approah. See for instance Proposition 1.7.1 in Samorodnitsky and Taqqu (1994).}
\usage{
rSaS(r, alpha, c = 1, mu = 0)
}
\arguments{
  \item{r}{
Number of observations.}
  \item{alpha}{
Index of stability; Number in (0,2)}
  \item{c}{
Scale parameter, c>0}
  \item{mu}{
Location parameter, any real number}
}
\details{
The characteristic function is

f(t) = e^(-c |t|^alpha)*e^(i*t*mu).

}
\references{
G. Samorodnitsky and M. Taqqu (1994). Stable Non-Gaussian Random Processes: Stochastic Models with Infinite Variance. Chapman & Hall, Boca Raton.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rSaS(10,.5)
}