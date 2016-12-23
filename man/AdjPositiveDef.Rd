% File src/library/base/man/AdjPositiveDef.Rd
\name{AdjPositiveDef}
\alias{AdjPositiveDef}
\title{Regularization of the Multivariate EWS to be Positive Definite}
\description{
  Regularizes the multivariate evolutionary wavelet spectrum 
  to ensure that each matrix for all locations and levels are 
  positive definite.
}
\usage{
  AdjPositiveDef(object, tol = 1e-10)
}

\arguments{
  \item{object}{A mvLSW object containing the smoothed and bias 
    corrected multivariate evolutionary wavelet spectrum estimate.}
  \item{tol}{Tolerance. At 1e-10 by default.}
}

\details{
  Given an multivariate evolutionary wavelet spectrum estimate, 
  smoothed and bias corrected, there is no guarantee that the 
  matrix \eqn{\hat{S}_{j,k}} for any level j and location k is 
  positive definite. This command assesses from the eigenvalues 
  whether each matrix estimate is positive definite and if not 
  attempts to regularize the estimate via \code{modchol}. If this 
  fails to determine a positive definite estimate, then a threshold 
  is applied to the eigenvalues of the orignal matrix at 
  \code{tol}.
}

\value{
  An object of class \code{mvLSW}, invisibly.
}

\references{
  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes. \emph{Signal Processing, 
  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
}

\seealso{
  \code{\link{mvEWS}}, \code{\link{mvLSW}}, \code{eigen}, 
  \code{\link{modchol}}.
}

\examples{
## Sample bivariate locally stationary time series
set.seed(100)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

## Raw EWS using the Haar wavelet
RawEWS_Haar <- RawPeriodogram(X, format = TRUE)
plot(RawEWS_Haar, style = 2, info = 1)

## Smoothed EWS using kernel "daniell"
SmoothEWS_Haar <- Smooth_EWS(RawEWS_Haar, kernel.param = 4)
plot(SmoothEWS_Haar, style = 2, info = 1)

## Correct for the estimator bias
CorrectEWS_Haar <- CorrectBias(SmoothEWS_Haar)
plot(CorrectEWS_Haar, style = 2, info = 1)

## Adjust estimate for positive definite matrices
EWSest_Haar <- AdjPositiveDef(CorrectEWS_Haar)
plot(EWSest_Haar, style = 2, info = 1)

## As above, but using Daubechies least asymmetric wavelet with 10 vanishing moments
RawEWS_DLA10 <- RawPeriodogram(X, filter.number = 10, 
   family = "DaubLeAsymm", format = TRUE)
SmoothEWS_DLA10 <- Smooth_EWS(RawEWS_DLA10, kernel.param = 2)
CorrectEWS_DLA10 <- CorrectBias(SmoothEWS_DLA10)
EWSest_DLA10 <- AdjPositiveDef(CorrectEWS_DLA10)

plot(RawEWS_DLA10, style = 2, info = 1)
plot(SmoothEWS_DLA10, style = 2, info = 1)
plot(CorrectEWS_DLA10, style = 2, info = 1)
plot(EWSest_DLA10, style = 2, info = 1)
}

\keyword{AdjPositiveDef}
