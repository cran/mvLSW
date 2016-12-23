% File src/library/base/man/Asymp_Quantile.Rd
\name{Asymp_Quantile}
\alias{Asymp_Quantile}
\title{Evaluate the asymptotic quantile of a mvEWS estimate}
\description{
  Evaluate the asymptotic quantile of a multivariate evolutionary 
  wavelet spectrum.
}
\usage{
  Asymp_Quantile(object, var = NULL, prob = 0.5, ...)
}

\arguments{
  \item{object}{A \code{mvLSW} object containing the multivariate 
    evolutionary wavelet spectrum estimate.}
  \item{var}{A mvLSW object containing the variance estimate of 
    the wavelet spectrum. This is specified as \code{NULL} by 
	default and is evaluated based on \code{object}}
  \item{prob}{Probability of length 1 within the range [0,1].}
  \item{...}{Additional arguments to be passed to the 
    \code{\link{varEWS}} command in calculating the spectral variance.}
}

\details{
  This command evaluates the quantile of the asymptotic Gaussian 
  distribution for the mvEWS estimate.
}

\value{
  An object of class \code{mvLSW}, invisibly.
}

\references{
  Park, T. (2014) Wavelet Methods for Multivariate Nonstationary 
  Time Series, PhD thesis, Lancaster University, pp. 91-111.
}

\seealso{
  \code{\link{mvEWS}}, \code{\link{mvLSW}}, \code{\link{varEWS}}.
}

\examples{
## Define evolutionary wavelet spectrum, structure only on level 2
Spec <- array(0, dim = c(3, 3, 8, 256))
Spec[1, 1, 2, ] <- 10
Spec[2, 2, 2, ] <- c(rep(5, 64), rep(0.6, 64), rep(5, 128))
Spec[3, 3, 2, ] <- c(rep(2, 128), rep(8, 128))
Spec[2, 1, 2, ] <- Spec[1, 2, 2, ] <- c(rep(0, 64), seq(0, 1, len = 128), rep(1, 64))
Spec[3, 1, 2, ] <- Spec[1, 3, 2, ] <- c(rep(-1, 128), rep(5, 128))
Spec[3, 2, 2, ] <- Spec[2, 3, 2, ] <- -0.5
EWS <- mvLSW(spectrum = Spec, filter.number = 1, family = "DaubExPhase", 
  min.eig.val = NA)

## Sample time series and estimate the EWS.
set.seed(10)
X <- rmvLSW(Spectrum = EWS)
EWS_X <- mvEWS(X, kernel.name = "daniell", kernel.param = 20)

## Evaluate asymptotic spectral variance 
SpecVar <- varEWS(EWS_X)

## Plot Estimate & 95/% confidence interval
Q025 <- Asymp_Quantile(object = EWS_X, var = SpecVar, prob = 0.025)
Q975 <- Asymp_Quantile(object = EWS_X, var = SpecVar, prob = 0.975)
plot(x = EWS_X, style = 2, info = 2, Int.lower = Q025, Int.upper = Q975)
}

\keyword{Asymp_Quantile}
