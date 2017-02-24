% File src/library/base/man/mvEWS.Rd
\name{mvEWS}
\alias{mvEWS}
\title{Multivariate Evolutionary Wavelet Spectrum}
\description{
  Calculates the Evolutionary Wavelet Spectrum (EWS) of a 
  multivariate locally stationary time series.
}
\usage{
  mvEWS(X, filter.number = 1, family = "DaubExPhase", 
    smooth = TRUE, type = "all", kernel.name = "daniell", 
    kernel.param = floor(sqrt(nrow(X))), optimize = FALSE, 
	smooth.Jset = NA, bias.correct = TRUE, tol = 1e-10, verbose = FALSE)
}

\arguments{
  \item{X}{A multivariate time series object of class \code{ts}, 
    \code{zoo}, \code{xts} or \code{matrix}. The length of the 
    time series must \eqn{2^J} for positive integer J.}
  \item{filter.number}{Integer number defining the number of 
    vanishing moments of the wavelet function. By default, 
	\code{filter.number=1} and so defining the Haar wavelet.}
  \item{family}{Specifies the family of wavelet. Only two 
    options are available, either \code{"DaubExPhase"} (default) 
	or \code{"DaubLeAsymm"}.}
  \item{smooth}{Logical, should EWS should be smoothed?}
  \item{type}{How should the smoothing be performed. If \code{"all"}
    (default) then the same smoothing kernel is applied to all 
	levels, else if \code{"by.level"} then a different smoothing 
	kernel is applied to each level.}
  \item{kernel.name}{Name of the smoothing kernel to be supplied 
    to \code{kernel()}. Kernel \code{"daniell"} is defined by 
	default.}
  \item{kernel.param}{Parameters to be passed to \code{kernel()}. 
    This argument must be a vector if \code{type="all"}, otherwise
	it must be a matrix with each column defining the kernel parameters 
	for each \code{log2(nrow(X))} levels from coarse to fine. If 
	the kernel name is either \code{"dirichlet"} or \code{"fejer"} 
	then \code{kernel.param} must have length 2 (or a matrix with 
	2 rows) which are supplied to \code{kernel()} as arguments 
	\code{m} and \code{r} respectively. Note that the width of the 
	kernel cannot be larger than the time series length. This is set by 
	default as the square root of the length of the time series.}
  \item{optimize}{Logical, should the smoothing be optimized? 
    If \code{FALSE} (default) then the smoothing is performed as 
	specified with \code{kernel.name} and \code{kernel.param}. 
	Otherwise, \code{kernel.param} defines the upper parameter bound 
	and the optimal kernel is determined by minimising the generalized 
	cross-validation gamma deviance criterion.}
  \item{smooth.Jset}{Integer vector indicating with levels to be 
    used in calculating the generalized cross-validation gamma deviance 
	criterion. This argument is only used if \code{type="all"} and is 
	set as \code{NA} by default, implying that all levels should be used.}
  \item{bias.correct}{Logical, should the correction be applied 
    to address the bias in the EWS estimator.}
  \item{tol}{Tolerance in applying matrix regularisation 
    to ensure each EWS matrix to be strictly positive definite.
	If \code{NA} then the threshold is not applied. This is 1e-10 by default.}
  \item{verbose}{Logical. Controls the printing of messages whist 
    the computations progress. Set as \code{FALSE} as default.}
}

\details{
  This command evaluates the evolutionary wavelet spectrum of a 
  multivariate time series. The order of operations are as follows:
  
  Calculate the non-decimated wavelet coefficients \eqn{\{d^{p}_{j,k}\}} 
  for levels j = 1,\ldots,J, locations k = 0,ldots,T-1 (T=\eqn{2^J}) and signals
  p = 1,\ldots,P(=\code{ncol(X)}). The raw periodogram matrices are 
  then evaluated by \eqn{I^{(p,q)}_{j,k} = d^{p}_{j,k}d^{q}_{j,k}} 
  between any signal pair p & q.
  
  The above estimator is inconsistent and so the matrix sequence is 
  smoothed: \eqn{\tilde{I}^{(p,q)}_{j,k} = \sum_i W_i I^{(p,q)}_{j,k+i}}. 
  The kernel weights \eqn{W_i} are derived from the \code{kernel} command 
  and satisfy \eqn{W_i=W_{-i}} and \eqn{\sum_i W_i = 1}. The optimal 
  parameter for the smoothing kernel is determined by minimising the 
  generalized cross-validation gamma deviance criterion:
  
  \deqn{GCV = \sum_{p,j} GCV(p,j)}
  \deqn{GCV(p,j) = (T(1-W_0)^2)^{-1} \sum_{k=0}^{T-1} q_k D^{(p,p)}_{j,k}}
  \deqn{D^{(p,p)}_{j,k} = R^{(p,p)}_{j,k} - log(R^{(p,p)}_{j,k}) - 1}
  \deqn{R^{(p,p)}_{j,k} = \tilde{I}^{(p,p)}_{j,k} / I^{(p,p)}_{j,k}}
  
  where \eqn{q_0=q_{T-1}=0.5} and \eqn{q_k=1} otherwise. Note that the 
  criterion is not applicable on the wavelet cross-spectrum. A theshold 
  is applied to \eqn{\tilde{I}} and \eqn{I} in order to produce a valid 
  criterion estimate. In addition, the first summation is only applied 
  over the signals if smoothing is to be applied on a by-level basis.
  
  The raw wavelet periodogram is also a biased estimator. A correction is 
  subsequently applied to the smoothed estimate as follows:
  
  \deqn{\hat{S}_{j,k} =\sum_{l=1}^{J} A_{j,l}^{-1} \hat{I}_{l,k}}

  Here, \eqn{A_{j,l}} denotes the wavelet autocorrelation inner product.
  
  Finally, a threshold is applied to the eigenvalues of the EWS 
  \eqn{\hat{S}_{j,k}} to ensure that the matrices are positive definite.

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
  \code{ts}, \code{wd}, \code{kernel}, \code{\link{as.mvLSW}}, \code{ipndacw}.
}

\examples{
## Sample bivariate locally stationary time series
set.seed(100)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

## Haar wavelet, apply same smoothing to all levels & optimize
EWS <- mvEWS(X, kernel.name = "daniell", kernel.param = 20, optimize = TRUE)
summary(EWS)
plot(EWS, style = 2, info = 1)

## Over smoothed EWS
EWS_smooth <- mvEWS(X, filter.number = 10, family = "DaubLeAsymm",
  kernel.name = "modified.daniell", kernel.param = c(5, 5), optimize = FALSE)
summary(EWS_smooth)
plot(EWS_smooth, style = 2, info = 1)
}

\keyword{mvEWS}
