% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_local_influence.R
\name{local_influence_plot.mixpoissonreg}
\alias{local_influence_plot.mixpoissonreg}
\title{Local Influence Plot Diagnostics for Mixed Poisson Regression Models}
\usage{
\method{local_influence_plot}{mixpoissonreg}(
  model,
  which = c(1, 2, 3, 4),
  caption = list("Case Weights Perturbation", "Hidden Variable Perturbation",
    "Mean Explanatory Perturbation", "Precision Explanatory Perturbation",
    "Simultaneous Explanatory Perturbation"),
  sub.caption = NULL,
  detect.influential = TRUE,
  n.influential = 5,
  draw.benchmark = FALSE,
  lty.benchmark = 2,
  type_plot = "h",
  curvature = c("conformal", "normal"),
  direction = c("canonical", "max.eigen"),
  parameters = c("all", "mean", "precision"),
  mean.covariates = NULL,
  precision.covariates = NULL,
  main = "",
  ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive(),
  labels.id = names(stats::residuals(model)),
  cex.id = 0.75,
  cex.oma.main = 1.25,
  cex.caption = 1,
  include.modeltype = TRUE,
  ...
)
}
\arguments{
\item{model}{a \code{mixpoissonreg} object.}

\item{which}{a list or vector indicating which plots should be displayed.     If a subset of the plots is required, specify a subset of the numbers 1:5, see caption below (and the 'Details') for the different kinds.}

\item{caption}{captions to appear above the plots; character vector or list of valid graphics annotations. Can be set to "" or NA to suppress all captions.}

\item{sub.caption}{common title-above the figures if there are more than one. If NULL, as by default, a possible abbreviated version of deparse(x$call) is used.}

\item{detect.influential}{logical. Indicates whether the benchmark should be used to detect influential observations and identify them on the plot. If there is no benchmark available,
the top 'n.influential' observations will be identified in the plot by their indexes.}

\item{n.influential}{interger. The maximum number of influential observations to be identified on the plot.}

\item{draw.benchmark}{logical. Indicates whether a horizontal line identifying the benchmark should be drawn.}

\item{lty.benchmark}{the line type of the benchmark if drawn.}

\item{type_plot}{what type of plot should be drawn. The default is 'h'.}

\item{curvature}{the curvature to be returned, 'conformal' for the conformal normal curvature (see Zhu and Lee, 2001 and Poon and Poon, 1999) or
'normal' (see Zhu and Lee, 2001 and Cook, 1986).}

\item{direction}{the 'max.eigen' returns the eigenvector associated to the largest eigenvalue of the perturbation matrix. The 'canonical' considers
the curvatures under the canonical directions, which is known as "total local curvature" (see Lesaffre and Verbeke, 1998). For conformal
normal curvatures both of them coincide. The default is 'canonical'.}

\item{parameters}{the parameter to which the local influence will be computed. The options are 'all', 'mean' and 'precision'.
This argument affects the 'case_weights' and 'hidden_variable' perturbation schemes. The default is 'all'.}

\item{mean.covariates}{a list or vector of characters containing the mean-explanatory variables to be used in the 'mean-explanatory' and 'simultaneous-explanatory'
perturbation schemes. If NULL, the 'mean-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
mean-related covariates. The default is NULL.}

\item{precision.covariates}{a list or vector of characters containing the precision-explanatory variables to be used in the 'precision-explanatory'
and 'simultaneous-explanatory'
perturbation schemes. If NULL, the 'precision-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
precision-related covariates. The default is NULL.}

\item{main}{character; title to be placed at each plot additionally (and above) all captions.}

\item{ask}{logical; if \code{TRUE}, the user is asked before each plot.}

\item{labels.id}{vector of labels, from which the labels for extreme points will be chosen. The default uses the observation numbers.}

\item{cex.id}{magnification of point labels.}

\item{cex.oma.main}{controls the size of the sub.caption only if that is above the figures when there is more than one.}

\item{cex.caption}{controls the size of caption.}

\item{include.modeltype}{logical. Indicates whether the model type ('NB' or 'PIG') should be displayed on the captions.}

\item{...}{other graphical arguments to be passed.}
}
\value{
Called for its side effects.
}
\description{
Local influence plots for mixed Poisson regression models. Currently the conformal normal and normal curvatures are available
under several perturbation schemes. The default is the conformal normal curvature since
it takes values on \eqn{[0,1]} and other nice properties (see Zhu and Lee, 2001 and Poon and Poon, 1999 for further details).
}
\details{
\code{local_influence.mixpoissonreg} provides local influence diagnostics for mixed Poisson regression models for all perturbation schemes considered in
Barreto-Souza and Simas (2016), for normal and conformal normal curvatures. Further, it is also provides results for the canonical directions, which is called
the total local influence (see Lesaffre and Verbeke, 1998), as well as for the direction of largest curvature, which is the direction of the eigenvector of the
perturbation matrix associated to the largest eigenvalue.

\code{local_influence_plot.mixpoissonreg} provides a plot of the local influence diagnostics. Each plot corresponds to a perturbation scheme. The first plot considers
the 'case-weights' perturbation; the second plot considers the 'hidden-variable' perturbation (which was introduced in Barreto-Souza and Simas, 2016); the third plot
considers the mean-explanatory perturbation; the fourth plot considers the precision-explanatory perturbation; the fifth plot considers the simultanous-explanatory perturbation.

For both \code{local_influence.mixpoissonreg} and \code{local_influence_plot.mixpoissonreg}, one can select which covariates will be perturbed in the 'mean-explanatory',
'precision-explanatory' and 'simultaneous-explanatory' perturbation schemes. These are chosen in the 'mean.covariates' and 'precision.covariates' arguments.

If one considers the total local influence, then Zhu and Lee (2001) provides benchmark for influential observations for all perturbation schemes. These are returned as
attributes in the returned list from \code{local_influence.mixpoissonreg}. When using the \code{local_influence_plot.mixpoissonreg}, only points above the benchmark
will be displayed. One can also set the option 'draw_benchmark' to TRUE to plot the benchmark line.
}
\examples{
\donttest{
data("Attendance", package = "mixpoissonreg")

daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
prog | gender + math + prog, data = Attendance)
local_influence_plot(daysabs_fit)

daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
prog | gender + math + prog, data = Attendance, envelope = 100)
local_influence_plot(daysabs_fit_ml, which = 2)
}

daysabs_progML <- mixpoissonregML(daysabs ~ prog | prog, data = Attendance)
local_influence_plot(daysabs_progML)

}
\references{
DOI:10.1007/s11222-015-9601-6 \doi{10.1007/s11222-015-9601-6}(Barreto-Souza and Simas; 2016)

Cook, R. D. (1986) \emph{Assessment of Local Influence.} Journal of the Royal Statistical Society. Series B (Methodological), Vol. 48, pp.133-169. \href{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}

Lesaffre, E. and Verbeke, G. (1998) \emph{Local Influence in Linear Mixed Models}. Biometrics, 54, pp. 570-582.

Poon, W.-Y. and Poon, Y.S. (1999) \emph{Conformal normal curvature and assessment of local influence.}  Journal of the Royal Statistical Society. Series B (Methodological), Vol. 61, pp.51-61. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}

Zhu, H.-T. and Lee, S.-Y. (2001) \emph{Local influence for incomplete data models.} Journal of the Royal Statistical Society. Series B (Methodological), Vol. 63, pp.111-126. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}
}
