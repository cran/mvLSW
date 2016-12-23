% File src/library/base/man/SmoothEvWaveSpec.Rd
\name{SmoothEvWaveSpec}
\alias{SmoothEvWaveSpec}
\alias{Smooth_EWS}
\alias{Smooth_GCV}
\title{Smooth the Raw Evolutionary Wavelet Periodogram}
\description{
  Smoothing of the multivariate raw evolutionary wavelet periodogram 
  estimate and determine the optimal smoothing kernel according 
  to the generalized cross-validation gamma deviance criterion.
}
\usage{
  Smooth_EWS(object, kernel.name = "daniell", kernel.param = 1, 
    optimize = FALSE, type = "all", level = NA)
	
  Smooth_GCV(object, kernel.name = "daniell", kernel.param = 1, 
    type = "all", level = NA, GCVonly = FALSE)
}

\arguments{
  \item{object}{A mvLSW object containing the raw periodogram estimate.}
  \item{kernel.name}{Name of the smoothing kernel to be supplied 
    to \code{kernel()}. Kernel \code{"daniell"} is defined by 
	default.}
  \item{kernel.param}{Parameters to be passed to \code{kernel()}. 
    This argument must be supplied as a vector that must have 
	length 2 if \code{kernel.name} is \code{"dirichlet"} or 
	\code{"fejer"}; respectively defining arguments \code{m} and 
	\code{r} for command \code{kernel}. Note that the width of the 
	kernel cannot be larger than the time series length.}
  \item{optimize}{Logical, should the smoothing be optimized? 
    If \code{TRUE} then \code{kernel.param} defines the upper bounds 
	for the parameter space. Set as \code{FALSE} by default.}
  \item{type}{The type of smoothing to apply, either \code{"all"} 
    (default) or \code{"by.level"}.}
  \item{level}{Defines which levels should be used in calculating the 
    generalized cross-validation gamma deviance criterion. If 
	\code{type="all"} and \code{level=NA} (default) then all levels 
	are used in the calculation. Otherwise, \code{level} must be 
	an integer vector indicating the required levels. If 
	\code{type="by.level"} then \code{level} must be a single integer.}
  \item{GCVonly}{Logical, should only the generalized 
    cross-validation gamma deviance criterion be return? 
    If \code{FALSE} (default) then the smoothed EWS estimate 
	is returned.}
}

\details{
  The raw multivariate evolutionary wavelet periodogram is an 
  inconsistent estimator and so must be smoothed. Notate the 
  periodogram estimate between signal pair p & q, level j and 
  location k as \eqn{I^{(p,q)}_{j,k}}. The smoothed estimate 
  is derived by:
  
  \deqn{\tilde{I}^{(p,q)}_{j,k} = \sum_i W_i I^{(p,q)}_{j,k+i}}
  
  The kernel weights \eqn{W_i} are derived from the \code{kernel} command 
  and satisfy \eqn{W_i=W_{-i}} and \eqn{\sum_i W_i = 1}. The optimal 
  parameter for the smoothing kernel is determined by minimising the 
  generalized cross-validation gamma deviance criterion:
  
  \deqn{GCV = \sum_{p,j} GCV(p,j)}
  \deqn{GCV(p,j) = (T(1-W_0)^2)^{-1} \sum_{k=0}^{T-1} q_k D^{(p,p)}_{j,k}}
  \deqn{D^{(p,p)}_{j,k} = R^{(p,p)}_{j,k} - log(R^{(p,p)}_{j,k}) - 1}
  \deqn{R^{(p,p)}_{j,k} = \tilde{I}^{(p,p)}_{j,k} / I^{(p,p)}_{j,k}}
  
  where \eqn{q_0=q_{T-1}=0.5} and \eqn{q_k=1} otherwise. Note that the 
  criterion is not applicable on the wavelet cross-spectrum. Furthermore,
  the first summation above is only applied over the signals if 
  smoothing is to be applied on a by-level basis.
  
  The command \code{Smooth_EWS} is an envelope function for 
  \code{Smooth_GCV}. The latter command can be used to evaluate 
  the generalized cross-validation gamma deviance criterion as 
  a function of the kernel parameters.
  
}

\value{
  \code{Smooth_EWS} returns a \code{mvLSW} object, invisibly, 
  containing the smoothed EWS.
  
  If \code{GCVonly=TRUE} then \code{Smooth_GCV} returns the 
  generalized cross-validation gamma deviance criterion. 
  Otherwise, a \code{mvLSW} object is returned invisibly 
  containing the newly smoothed periodogram for the specified 
  levels. The periodogram for all other levels  that are not 
  requested to be smoothed are unchanged from \code{object}.
}

\references{
  Ombao, H.C., Raz, J.A. Strawderman, R.L. and von Sachs, R. (2001) 
  A Simple Generalised Crossvalidation Method of Span Selection 
  for Periodogram Smoothing. \emph{Biometrika} \strong{88}(4) pp. 1186-1192.

  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes. \emph{Signal Processing, 
  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
}

\seealso{
  \code{\link{mvEWS}}, \code{kernel}, \code{\link{mvLSW}}, 
  \code{\link{RawPeriodogram}}.
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

## GCV as a function of kernel parameters (Level 1 only)
GCV <- rep(NA, 10)
for(m in 1:10){
  GCV[m] <- Smooth_GCV(RawEWS_Haar, kernel.param = m, level = 1, GCVonly = TRUE)
}
plot(1:10, GCV)
}

\keyword{SmoothEvWaveSpec, SmoothEWS, Smooth_GCV}
