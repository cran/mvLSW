% File src/library/base/man/CorrectBias.Rd
\name{CorrectBias}
\alias{CorrectBias}
\title{Bias Correction to Multivariate Evolutionary Wavelet Spectrum Estimate}
\description{
  Bias correction to the multivariate evolutionary wavelet 
  spectrum estimate.
}
\usage{
  CorrectBias(object)
}

\arguments{
  \item{object}{A mvLSW object containing the biased multivariate 
    evolutionary wavelet spectrum estimate.}
}

\details{
  Let \eqn{\hat{I}_{l,k}} define smoothed evolutionary wavelet spectrum
  estimate (see \code{mvEWS} for more details) which is a biased. The 
  correction to the estimator is calculated as follows:
  
  \deqn{\hat{S}_{j,k} =\sum_{l=1}^{J} (A^{-1})_{j,l} \hat{I}_{l,k}}

  Here, \eqn{A_{j,l}} defines the wavelet autocorrelation inner product.
}

\value{
  A \code{mvLSW} object, invisibly.
}

\references{
  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes. \emph{Signal Processing, 
  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
}

\seealso{
  \code{\link{mvEWS}}, \code{\link{mvLSW}}, \code{ipndacw}.
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
plot(CorrectEWS_DLA10, style=2, info = 1)
plot(EWSest_DLA10, style = 2, info = 1)
}

\keyword{CorrectBias}