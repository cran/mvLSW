% File src/library/base/man/plot.mvLSW.Rd
\name{plot.mvLSW}
\alias{plot.mvLSW}

\title{Plot mvLSW Object}
\description{
  Plot the data contained within a mvLSW object based on the 
  requested format.
}
\usage{
  \method{plot}{mvLSW}(x, style = 1, info = NULL, Int.lower = NULL, 
    Int.upper = NULL, diag = TRUE, sub = "Spectrum", ...)
}

\arguments{
  \item{x}{A \code{mvLSW} object.}
  \item{style}{Index stating the type of plotting format for 
    the mvLSW object. (See details.)}
  \item{info}{Vector containing the channel and/or level indices 
    defineing the slice throught the mvEWS according to the 
	requested plotting \code{style}. (See details.)}
  \item{Int.lower, Int.upper}{\code{mvLSW} objects respectively 
    containing the lower and upper values for the interval 
	to be drawn. Both arguments must be supplied to be drawn. 
	By default, both arguments are \code{NULL}.}
  \item{diag}{Logical, should the diagonal pannels be drawn when 
    \code{style=2}. Ideally this should be \code{FALSE} if 
	\code{object} contains the coherence. Set to \code{TRUE} by default.}
  \item{sub}{Plot subtitle. Set to \code{"Spectrum"} by default.}
  \item{...}{Additional graphical parameters.}
}

\details{
  This command plots the data contained within the \code{mvLSW} based 
  on requested plotting style.
  
  Plotting style \code{style=1} with information \code{info=c(p,q,j)} 
  generates a single plot for a 
  specified signal pair \code{p} & \code{q} and level \code{j}.
  
  Plotting style \code{style=2} with information \code{info=j} 
  creates a set of plots from \code{x} for all channel pairs in a 
  lower-triangular pannel corresponding to the specified level $j$.
  If \code{diag=FALSE} then the plots along the diagonal are suppressed,
  which is ideal when \code{x} contain coherence estimates.

  Plotting style \code{style=3} with information \code{info=c(p,q)}
  creates a set of plots from \code{x} for all levels (from fine 
  to coarse) for channel pair \code{p} and \code{q}.
  
  Finally, the plotting style \code{style=4} with information 
  \code{info=j} presents the same infromation from \code{x} as
  for the previous case, bit in a compact matrix format. Please 
  refer to \code{image.plot} from the \code{fields} library for 
  additional information on this plotting style.
  
  Both arguments \code{Int.lower} and \code{Int.upper} must be supplied 
  in order to draw a polygon to indicate the interval estimate. These 
  arguments are ignored in the case \code{style=4}.
}

\value{
  Generates a plotting window. No data is returned.
}

\seealso{
  \code{plot.default}, \code{image.plot}, \code{\link{as.mvLSW}}, 
  \code{\link{mvEWS}}, \code{\link{coherence}}.
}

\examples{
## Define evolutionary wavelet spectrum, structure only on level 2
Spec <- array(0, dim=c(3, 3, 8, 256))
Spec[1, 1, 2, ] <- 10
Spec[2, 2, 2, ] <- c(rep(5, 64), rep(0.6, 64), rep(5, 128))
Spec[3, 3, 2, ] <- c(rep(2, 128), rep(8, 128))
Spec[2, 1, 2, ] <- Spec[1, 2, 2, ] <- c(rep(0, 64), seq(0, 1, len = 128), rep(1, 64))
Spec[3, 1, 2, ] <- Spec[1, 3, 2, ] <- c(rep(-1, 128), rep(5, 128))
Spec[3, 2, 2, ] <- Spec[2, 3, 2, ] <- -0.5
EWS <- as.mvLSW(x = Spec, filter.number = 1, family = "DaubExPhase",
  min.eig.val = NA)

## Sample time series and estimate the EWS and coherence.
set.seed(10)
X <- rmvLSW(Spectrum = EWS)
EWS_X <- mvEWS(X, kernel.name = "daniell", kernel.param = 20)
RHO_X <- coherence(EWS_X, partial = FALSE)

## Evaluate asymptotic spectral variance & 95\% confidence interval
SpecVar <- varEWS(EWS_X)
Q025 <- Asymp_Quantile(object = EWS_X, var = SpecVar, prob = 0.025)
Q975 <- Asymp_Quantile(object = EWS_X, var = SpecVar, prob = 0.975)

## Plot evolutionary wavelet spectrum between signals 1 & 3 at level 2
plot(x = EWS_X, style = 1, info = c(1, 3, 2), Int.lower = Q025, Int.upper = Q975)

## Plot coherence between signals 1 & 3 at level 2
plot(x = RHO_X, style = 1, info = c(1, 3, 2), ylab = "Coherence")

## Evolutionary wavelet spectrum panel plot for level 2
plot(x = EWS_X, style = 2, info = 2, Int.lower = Q025, Int.upper = Q975)

## Panel plot of coherence for level 2
plot(x = RHO_X, style = 2, info = 2, diag = FALSE, ylab = "Coherence")

## Plot evolutionary wavelet spectrum for signal pair 1 & 3 at all levels
plot(x = EWS_X, style = 3, info = c(1, 3), Int.lower = Q025, Int.upper = Q975)

## Plot coherence for signal pair 1 & 3 at all levels
plot(x = RHO_X, style = 3, info = c(1, 3), ylab = "Coherence")

## Image plot for coherence between signals 1 & 3
plot(x = RHO_X, style = 4, info = c(1, 3), sub = "Coherence")
}

\keyword{plot.mvLSW}
