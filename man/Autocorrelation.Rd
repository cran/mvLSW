% File src/library/base/man/Autocorrelation.Rd
\name{Autocorrelation}
\alias{Autocorrelation}
\alias{AutoCorrIP}
\title{Wavelet Autocorrelation Functions}
\description{
  Wavelet autocorrelation functions between two wavelet functions 
  associated to level j & l for all pairs j, l = 1,\ldots,J and 
  associated inner product.
}
\usage{
  AutoCorrIP(J, filter.number = 1, family = "DaubExPhase", crop = TRUE)
}

\arguments{
  \item{J}{Number of levels.}
  \item{filter.number}{Number of vanishing moments for wavelet.}
  \item{family}{The family of wavelet, either \code{"DaubExPhase"} 
    or \code{"DaubLeAsymm"}. The Haar wavelet is defined as default.}
  \item{crop}{Logical, should the output of \code{AutoCorrIP} be 
    cropped such that the first dimension is limited to -\eqn{2^J}:\eqn{2^J}.
	This is set at \code{TRUE} by default.}
}

\details{
  This function returns the wavelet autocorrelation function 
  between wavelet functions at level pair j & l.
  
  Let \eqn{\psi(x)} denote the mother wavelet and the wavelet 
  defined for level j as \eqn{\psi_{j,k}(x) = 2^{j/2}\psi(2^{j}x-k)}.
  The wavelet autocorrelation function between levels j & l 
  is therefore:
  
  \deqn{\Psi_{j,l}(\tau) = \sum_\tau \psi_{j,k}(0)\psi_{l,k_\tau}(0)}
  
  Here, integer \eqn{\tau} defines the offset of the latter 
  wavelet function relative to the first.
  
  The inner product of this wavelet autocorrelation function is 
  defined as follows for level indices j, l & h and offset \eqn{\lambda}:
  
  \deqn{A^{\lambda}_{j,l,h} = \sum_{\tau} \Psi_{j,l}(\lambda - \tau) \Psi_{h,h}(\tau)}
  
}

\value{
    A 4D array (invisibly returned) of order 
    LxJxJxJ where L depends on the specified wavelet function. 
	If \code{crop=TRUE} then L=\eqn{2^{J+1}}+1. The first dimension 
	defines the offset \eqn{\lambda}, whilst the second to 
	fourth dimensions identify the levels indexed by j, l & h 
	respectively.
}

\references{
  Fryzlewicz, P. and Nason, G. (2006) HaarFisz estimation of 
  evolutionary wavelet spectra. \emph{Journal of the Royal 
  Statistical Society. Series B}, \strong{68}(4) pp. 611-634.
}

\seealso{
  \code{ipndacw}.
}

\examples{
## Plot Haar autocorrelation wavelet functions inner product
AInnProd <- AutoCorrIP(J = 8, filter.number = 1, family = "DaubExPhase")
\dontrun{
MaxOffset <- 2^8
for(h in 6:8){
  x11()
  par(mfrow = c(3, 3))
  for(l in 6:8){
    for(j in 6:8){
      plot(-MaxOffset:MaxOffset, AInnProd[, j, l, h], type = "l", 
        xlab = "lambda", ylab = "Autocorr Inner Prod", 
        main = paste("j :", j, "- l :", l, "- h :", h))
    }
  }
}
}

## Special case relating to ipndacw function from wavethresh package
Amat <- matrix(NA, ncol = 8, nrow = 8)
for(j in 1:8) Amat[, j] <- AInnProd[2^8 + 1, j, j, ]
round(Amat, 5)
round(ipndacw(J = -8, filter.number = 1, family = "DaubExPhase"), 5)
}

\keyword{PsiJLmat, AutoCorrIP, Autocorrelation}
