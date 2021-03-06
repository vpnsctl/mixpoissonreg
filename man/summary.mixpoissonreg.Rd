% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_diagnostic_functions.R
\name{summary.mixpoissonreg}
\alias{summary.mixpoissonreg}
\title{Summary Method for \code{mixpoissonreg} Objects.}
\usage{
\method{summary}{mixpoissonreg}(object, ...)
}
\arguments{
\item{object}{an object of class "mixpoissonreg" containing results from the fitted model.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class \code{summary_mixpoissonreg} containing several
informations of a \emph{mixpoissonreg} object.
}
\description{
Function providing a summary of results related to mixed Poisson regression models.
}
\examples{
\donttest{
data("Attendance", package = "mixpoissonreg")

daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
prog | gender + math + prog, data = Attendance)
summary(daysabs_fit)

daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
prog | gender + math + prog, data = Attendance)
summary(daysabs_fit_ml)
}

daysabs_prog <- mixpoissonreg(daysabs ~ prog, data = Attendance)
summary(daysabs_prog)

}
\seealso{
\code{\link{plot.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}},
\code{\link{local_influence_plot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
}
