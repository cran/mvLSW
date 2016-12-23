% File src/library/base/man/plot.mvLSW.Rd
\name{plot.mvLSW}
\alias{plot.mvLSW}
\alias{single_pqj_plot}
\alias{panel_pq_plot}
\alias{panel_j_plot}
\alias{matrix_j_plot}

\title{Plot mvLSW Object}
\description{
  Plot the data contained within a mvLSW object based on the 
  requested format.
}
\usage{
  \method{plot}{mvLSW}(x, style = 1, info = NULL, Int.lower = NULL, 
    Int.upper = NULL, diag = TRUE, sub = "Spectrum", ...)

  single_pqj_plot(object, Int.lower = NULL, Int.upper = NULL, p = 1, 
    q = 1, j = 1, ...)
  
  panel_j_plot(object, Int.lower = NULL, Int.upper = NULL, j = 1, 
    diag = TRUE, ...)

  panel_pq_plot(object, Int.lower = NULL, Int.upper = NULL, p = 1, 
    q = 1, ...)

  matrix_j_plot(object, p = 1, q = 1, sub = "Spectrum", ...)
}

\arguments{
  \item{x, object}{A \code{mvLSW} object.}
  \item{style}{Index stating the type of plotting format for 
    the mvLSW object. (See details.)}
  \item{info}{Vector containing the channel and level indices 
    to be passed to the requested plotting \code{style}. 
	(See details.)}
  \item{Int.lower, Int.upper}{\code{mvLSW} objects respectively 
    containing the lower and upper values for the interval 
	to be drawn. Both arguments must be supplied to be drawn. 
	By default, both arguments are \code{NULL}.}
  \item{p,q}{Signal indices.}
  \item{j}{Level index.}
  \item{diag}{Logical, should the diagonal pannels be drawn when calling 
    \code{panel_j_plot}. Ideally this should be \code{FALSE} if 
	\code{object} contains the coherence. Set to \code{TRUE} by default.}
  \item{sub}{Subtitle that is supplied to \code{matrix_j_plot}. Set to 
    \code{"Spectrum"} by default.}
  \item{...}{Additional graphical parameters.}
}

\details{
  These commands plot the data contained within the \code{mvLSW} based 
  on requested plotting style.
  
  Plotting style \code{style=1} invokes the command 
  \code{single_pqj_plot} that generates a single plot for a 
  specified signal pair \code{p} & \code{q} and level \code{j}; 
  supplied via the argument \code{info=c(p,q,j)} if called via the 
  \code{plot.mvLSW} command.
  
  Plotting style \code{style=2} invokes the command 
  \code{pannel_j_plot} that creates a panel lower-triangular, 
  including diagonal, where the data in \code{x} is plotted for 
  each channel pair for the specified level \code{j}. If called 
  via the command \code{plot.mvLSW} then the information argument 
  must be \code{info=j}.

  Plotting style \code{style=3} invokes the command 
  \code{pannel_pq_plot} that generates a panel of plots that 
  relate to the channel pair \code{p} and \code{q} for all levels from 
  fine to coarse. Id called via \code{plot.mvLSW} then the information 
  argument must be \code{info=c(p,q)}.
  
  Finally, the plotting style \code{style=4} involves the command 
  \code{matrix_j_plot} that subsequently calls \code{image.plot}
  from the \code{fields} library. A plot of the data 
  from the previous case (\code{style=3}) is generated but 
  in a compact matrix format.
  
  Both arguments \code{Int.lower} and \code{Int.upper} must be supplied 
  in order to draw a polygon to indicate the interval estimate.
}

\value{
  Generates a plotting window. No data is returned.
}

\seealso{
  \code{hsv}, \code{layout}, \code{par}, \code{plot.default}, 
  \code{image.plot}, \code{\link{mvLSW}}, \code{\link{mvEWS}}, 
  \code{\link{coherence}}, \code{tim.colors}.
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
EWS <- mvLSW(spectrum = Spec, filter.number = 1, family = "DaubExPhase",
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
