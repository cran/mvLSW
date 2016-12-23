% File src/library/base/man/rmvLSW.Rd
\name{rmvLSW}
\alias{rmvLSW}
\title{Sample a Multivariate Locally Stationary Wavelet Process}
\description{
  Sample a Multivariate Locally Stationary Wavelet Process.
}
\usage{
  rmvLSW(Transfer = NULL, Spectrum = NULL, noiseFN = rnorm, ...)
}

\arguments{
  \item{Transfer}{A \code{mvLSW} object containing the set of 
    transfer function matrices of the process.}
  \item{Spectrum}{A \code{mvLSW} object containing the evolutionary 
    wavelet spectrum of the process. This argument is only used if 
	\code{Transfer} is not supplied.}
  \item{noiseFN}{The function for sampling the innovations.}
  \item{...}{Optional arguments to be passed to \code{noiseFN}.}
}

\details{
  Samples a single multivariate locally stationary wavelet time 
  series for the given set of transfer function matrices. These 
  are assumed to be lower-triangular (including diagonal) matrices.
  If the evolutionary wavelet spectrum is supplied instead, then 
  this is pre-processed by \code{Spectrum2Transfer()} to obtain 
  the transfer function matrices.
  
  The \code{Transfer} and \code{Spectrum} are both \code{mvLSW} 
  objects and therefore contain information about defining the 
  wavelet function.
  
  The innovation process is assumed to be second order stationary 
  with expectation zero, orthogonal and unit variance. The first argument 
  of \code{noiseFN} must be \code{n} and define the number of samples 
  to generate. The function must also return a numerical vector 
  of length \code{n}.

}

\value{
  A \code{ts} matrix object of a locally stationary time series. The 
  columns of the matrix correspond to different signals and the rows 
  identify the time location.
}

\references{
  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes. \emph{Signal Processing, 
  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
}

\seealso{
  \code{\link{mvLSW}}, \code{\link{Spectrum2Transfer}}, 
  \code{rnorm}, \code{AvBasis}, \code{ts}.
}

\examples{
## Define evolutionary wavelet spectrum, structure only on level 2
Spec <- array(0, dim = c(3, 3, 8, 256))
Spec[1, 1, 2, ] <- 10
Spec[2, 2, 2, ] <- c(rep(5, 64), rep(0.6, 64), rep(5, 128))
Spec[3, 3, 2, ] <- c(rep(2, 128), rep(8, 128))
Spec[2, 1, 2, ] <- Spec[1, 2, 2, ] <- c(rep(0, 64), seq(0, 1, len=128), rep(1, 64))
Spec[3, 1, 2, ] <- Spec[1, 3, 2, ] <- c(rep(-1, 128), rep(5, 128))
Spec[3, 2, 2, ] <- Spec[2, 3, 2, ] <- -0.5

## Define Haar wavelet function and create mvLSW object
EWS <- mvLSW(spectrum = Spec, filter.number = 1, family = "DaubExPhase",
  min.eig.val = NA)
plot(EWS, style = 2, info = 2)

## Sample with Gaussian innovations
set.seed(10)
X <- rmvLSW(Spectrum = EWS)
plot(X)

## Define smoother wavelet function and create mvLSW object
EWS2 <- mvLSW(spectrum = Spec, filter.number = 10, family = "DaubExPhase")

## Sample with logisitc innovations
set.seed(10)
X2 <- rmvLSW(Spectrum = EWS2, noiseFN = rlogis, scale = sqrt(3)/pi)
plot(X2)
}

\keyword{rmvLSW}
