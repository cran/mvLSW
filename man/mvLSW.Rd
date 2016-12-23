% File src/library/base/man/mvLSW.Rd
\name{mvLSW}
\alias{mvLSW}
\alias{is.mvLSW}
\title{Multivariate Locally Stationary Wavelet Object}
\description{
  Constructs a multivariate locally stationary wavelet (mvLSW) object.
}
\usage{
  mvLSW(spectrum, filter.number = 1, family = "DaubExPhase",
    smooth.type = "all", smooth.kernel = kernel("daniell", 0), 
    bias.correct = FALSE, min.eig.val = -Inf, names = NULL)
  
  \method{is}{mvLSW}(object)
}

\arguments{
  \item{spectrum}{4D array of order PxPxJxT where P is the number of 
    signals of the time series of length T such that T=\eqn{2^J} for 
    some positive integer J defining the levels of the mvLSW object.}
  \item{family}{Specifies the family of wavelet. Only two 
    options are available, either \code{"DaubExPhase"} (default) 
	or \code{"DaubLeAsymm"}.}
  \item{filter.number}{Integer number defining the number of 
    vanishing moments of the wavelet function. By default, 
	\code{filter.number=1} and so defining the Haar wavelet.}
  \item{smooth.type}{What type of smoothing regime has been 
    applied. Either \code{"all"} (default) if the smoothing 
	method been applied to all levels. Otherwise \code{"by.level"}, 
	a different smoothing method is applied to each level.}
  \item{smooth.kernel}{Definition of the smoothing kernel from 
    \code{kernel()}. By default, the identity kernel is defined.}
  \item{bias.correct}{Logical, has a bias correction been applied 
    to the data. \code{FALSE} by default.}
  \item{min.eig.val}{Minimum eigenvalue from spectral matrices across 
    all levels and locations, set at \code{-Inf} by default. If \code{NA},
	then the minimum eigenvalue is calculated.}
  \item{names}{Character vector containing the signal names of the 
    multivariate time series.}
  \item{object}{An R object.}
}

\details{
  Construction of a multivariate locally stationary classed object 
  that contains all information about the constructions of various 
  multivariate wavelet estimates.
  
  The command \code{is.mvLSW} checks that the supplied R object 
  is a valid \code{mvLSW} object in that its structure and 
  contents are as expected.  
}

\value{
  The \code{mvLSW} command invisibly returns a list with the 
  following items:
  \item{spectrum}{A 4D array containing the data relating to the 
    estimate of interest.}
  \item{Information}{List containing information about the 
    estimate relating to the methodology.}
  
  The list \code{Information} contains:
  \item{names}{Character vector containing the signal names.}
  \item{dimensions}{A list containing items \code{P} - the number 
    of signals forming of the time series, \code{T} - the length 
	of the time series and \code{J} - the number of levels in the 
	wavelet transform of the data.}
  \item{wavelet}{A list containing the \code{filter.number} and 
    \code{family} of the wavelet used in the transformation.}
  \item{smooth}{A list containing information about the smoothing 
    applied to the estimate. Items of the list are: 
	\code{smooth.type} - name of the smoothing regime, 
	\code{smooth.kernels} - a \code{tskernel} class from the 
	command \code{kernel}, \code{GCV} - generalized cross-validation 
	gamma deviance criterion of the smoothing and \code{smooth.eps} - 
	smoothing threshold. If \code{smooth.type="by.level"}, then 
	\code{smooth.kernels} is a list containing the \code{tskernel} 
	object for each level from fine to coarse. In addition, 
	\code{GCV} is a length J vector containing the criterion 
	estimate for each level from fine to coarse.}
  \item{correction}{A list containing \code{bias.correction} and 
    \code{min.eig.val}.}
  
  The command \code{is.mvLSW} returns \code{TRUE} if the supplied 
  object is a \code{mvLSW} object as described above. Otherwise, the 
  command returns \code{FALSE}.
}

\seealso{
  \code{\link{mvEWS}}, \code{\link{varEWS}}, \code{kernel}.
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

## Define EWS as mvLSW object
EWS <- mvLSW(spectrum = Spec, filter.number = 1, family = "DaubExPhase", 
  names = c("A", "B", "C"), min.eig.val = NA)
is.mvLSW(EWS)
plot(EWS, style = 2, info = 2)
}

\keyword{mvLSW}