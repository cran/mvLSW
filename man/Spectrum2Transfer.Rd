% File src/library/base/man/Spectrum2Transfer.Rd
\name{Spectrum2Transfer}
\alias{Spectrum2Transfer}
\title{Convert between EWS and Transfer function matrices}
\description{
  Convert between multivariate evolutionary wavelet spectrum 
  and the set of transfer function matrices.
}
\usage{
  Spectrum2Transfer(object, S2V = TRUE)
}

\arguments{
  \item{object}{A \code{mvLSW} object containing either the EWS 
    the set of transfer function matrices.}
  \item{S2V}{Logical, if \code{TRUE} (default) then \code{object} 
    is the EWS and the set of transfer function matrices are to 
	be derived. If \code{FALSE} the \code{object} is the set of 
	transfer function matrices and the converse transformation is 
	derived.}
}

\details{
  If the EWS is supplied, then the set of transfer function matrices 
  are derived by the Choleski factorization of a real symmetric 
  semi-positive definite square matrix. In the cases where the matrix 
  is semi-definite, then the Choleski factorization is applied to the 
  submatrix that is positive definite and the remaining lower triangular 
  elements are populated such that the resulting matrix is a valid 
  factorization.
   
  Conversely, if the set of transfer function matrices are supplied, 
  then the EWS are derived by squaring the matrices.
}

\value{
  A \code{mvLSW} object containing either the EWS or set of transfer 
  function matrices depending on the specified transformation 
  direction.
}

\references{
  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes. \emph{Signal Processing, 
  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
}

\seealso{
  \code{chol}, \code{\link{mvLSW}}, \code{\link{mvEWS}}.
}

\examples{
## Define evolutionary wavelet spectrum, structure only on level 2
Spec <- array(0, dim=c(3, 3, 8, 256)) ## Ensure all are positive def.
Spec[1, 1, 2, ] <- 10
Spec[2, 2, 2, ] <- c(rep(5, 64), rep(0.6, 64), rep(5, 128))
Spec[3, 3, 2, ] <- c(rep(2, 128), rep(8, 128))
Spec[2, 1, 2, ] <- Spec[1, 2, 2, ] <- c(rep(0, 64), seq(0, 1, len=128), rep(1, 64))
Spec[3, 1, 2, ] <- Spec[1, 3, 2, ] <- c(rep(-1, 128), rep(5, 128))
Spec[3, 2, 2, ] <- Spec[2, 3, 2, ] <- -0.5

## Define EWS as mvLSW object
EWS <- mvLSW(spectrum = Spec, filter.number = 1, family = "DaubExPhase",
  min.eig.val = NA)
plot(EWS, style = 2, info = 2)

## EWS to Transfer function matrices
Transfer <- Spectrum2Transfer(object = EWS, S2V = TRUE)

## Transfer function matrices to EWS
EWS2 <- Spectrum2Transfer(object = Transfer, S2V = FALSE)
plot(EWS2, style = 2, info = 2)
}

\keyword{Spectrum2Transfer}
