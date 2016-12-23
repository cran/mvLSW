% File src/library/base/man/RawPeriodogram.Rd
\name{RawPeriodogram}
\alias{RawPeriodogram}
\title{Raw Multivariate Wavelet Periodogram}
\description{
  Calculates the raw evolutionary wavelet periodogram of a 
  multivariate locally stationary time series.
}
\usage{
  RawPeriodogram (X, filter.number = 1, family = "DaubExPhase", 
    format = TRUE)
}

\arguments{
  \item{X}{A multivariate time series object of class \code{ts}. 
    The length of the time series must \eqn{2^J} for positive integer J.}
  \item{family}{Specifies the family of wavelet. Only two 
    options are available, either \code{"DaubExPhase"} (default) 
	or \code{"DaubLeAsymm"}.}
  \item{filter.number}{Integer number defining the number of 
    vanishing moments of the wavelet function. By default, 
	\code{filter.number=1} and so defining the Haar wavelet.}
  \item{format}{Logical, should the raw periodogram be returned as
    a MvLSW object? Otherwise, the raw periodogram is returned as 
	a 4D array. Set to \code{TRUE} by default.}
}

\details{
  First, each signal is transformed to obtain the non-decimated wavelet 
  coefficients \eqn{\{d^{p}_{j,k}\}} for levels j = 1,\ldots,J, 
  locations k = 0,ldots,T-1 (T=\eqn{2^J}) and signals p = 1,\ldots,P(=\code{ncol(X)}). 
  The raw periodogram is then evaluated between any signal pair p & q as:
  
  \deqn{I^{(p,q)}_{j,k} = d^{p}_{j,k}d^{q}_{j,k}} 

  Note that the estimator is inconsistent and biased.
}

\value{
  If \code{format=TRUE} then a \code{MvLSW} class object object is returned, 
  otherwise the raw periodogram is returned as a 4D array. In both cases, the 
  result is returned invisibly.
}

\references{
  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes. \emph{Signal Processing, 
  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
}

\seealso{
  \code{wd}, \code{\link{mvEWS}}, \code{\link{mvLSW}}.
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

\keyword{RawPeriodogram}
