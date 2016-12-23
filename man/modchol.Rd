% File src/library/base/man/modchol.Rd
\name{modchol}
\alias{modchol}
\alias{PosDefEst}
\title{Modified cholesky decomposition.}
\description{
  Evaluates the modified cholesky decomposition for a 
  near-positive definite matrix.
}
\usage{
  modchol(A, tol = .Machine$double.eps)
  PosDefEst(A, tol = .Machine$double.eps)
}

\arguments{
  \item{A}{Square symmetric real matrix.}
  \item{tol}{Tolerance. At \code{.Machine$double.eps} by default.}
}

\details{
  \code{modchol} evaluates the modified cholesky decomposition for a 
  near-positive defintie matrix. See reference for procedure. 
  \code{PosDefEst} is a wrapper command for modchol that returns the 
  positive definite matrix approximation to \code{A}.
}

\value{
  The command \code{modchol} returns either a list, described below, 
  or \code{NULL} if the routine attempts to take the square-root of 
  a negative number (hence \code{A} is assumed to not be near-positive 
  definite). If a decomposition is evaluated successfully, then item in 
  the returned list are:
  \item{L}{Lower triangular modified cholesky estimate (subject to 
    pivoting)}
  \item{g}{Lower Gerschgorin bounds}
  \item{P}{Pivot transformatiton}
  \item{E}{Diagonal elements of modifying matrix to make \code{A} 
    positive definite}
  \item{delta}{maximum adjustment applied in routine}
  
  The command \code{PosDefEst} returns the positive definite symmetric 
  square matrix if \code{modchol} is successful in determine a modified 
  cholesky estimate. Otherwise \code{PosDefEst} returns \code{NULL}.
}

\references{
  Schnabel, R. B. and Eskow, E. 1999. "A Revised Modified Cholesky 
  Factorization Algorithm." SIAM Journal on optimization, 9(4), 1135-1148.
}

\seealso{
  \code{\link{AdjPositiveDef}}.
}

\examples{
  ##Positive Defintie Case
  M <- 0.9^apply(as.matrix(1:4),1,function(a){abs(a-(1:4))})
  modchol(M) ## Different estimate from chol(M) due to pivoting
  PosDefEst(M) ## Same as M

  ##Not Positive Definite 
  M2 <-M
  M2[2,2] <- 0.7
  eigen(M2)$values  ##has a negative eigenvalue
  L <- modchol(M2)
  L
  M2_est <- PosDefEst(M2)
  M2_est
  M2_est - M2 ## == diag(L$E)
}

\keyword{modchol,PosDefEst}