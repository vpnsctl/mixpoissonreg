% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_tidy_mixpoissonreg.R
\name{augment.mixpoissonreg}
\alias{augment.mixpoissonreg}
\alias{augment}
\title{Augment data with information from a \code{mixpoissonreg} object}
\usage{
\method{augment}{mixpoissonreg}(
  x,
  data = stats::model.frame(x),
  newdata = NULL,
  type.predict = c("response", "link", "precision", "variance"),
  type.residuals = c("pearson", "score"),
  se_fit = FALSE,
  conf_int = TRUE,
  pred_int = FALSE,
  ...
)
}
\arguments{
\item{x}{A \code{mixpoissonreg} object.}

\item{data}{A \code{\link[base:data.frame]{base::data.frame}} or \code{\link[tibble:tibble]{tibble::tibble()}} containing the original data that was used to produce the object x.}

\item{newdata}{A \code{\link[base:data.frame]{base::data.frame}} or \code{\link[tibble:tibble]{tibble::tibble()}} containing all the original predictors used to create x.
Defaults to \code{NULL}, indicating that nothing has been passed to \code{newdata}. If \code{newdata} is specified, the data argument will be ignored.}

\item{type.predict}{Type of prediction. The options are 'response', 'link', 'precision' and 'variance'. The default is "response".}

\item{type.residuals}{Type of residuals. The options are 'pearson' and 'score'. The default is 'pearson'.}

\item{se_fit}{Logical indicating whether or not a .se.fit column should be added to the augmented output. If TRUE, it only returns a non-NA value if type of prediction is 'link'.}

\item{conf_int}{Logical indicating whether or not confidence intervals for the fitted variable with type chosen from type.predict should be built. The available type options are
'response' and 'link'.}

\item{pred_int}{Logical indicating whether or not prediction intervals for future observations should be built. It only works with type.predict = 'response'. The arguments
\code{level}, \code{nsim_pred}, \code{nsim_pred_y} are passed through additional arguments '...'. Notice that this can be computationally intensive.}

\item{...}{Additional arguments. Possible additional arguments are \code{level}, \code{nsim_pred}, \code{nsim_pred_y}, that are passed to \code{predict.mixpoissonreg} function.}
}
\value{
A \code{\link[tibble:tibble]{tibble::tibble()}} with columns:
\itemize{
\item \code{.cooksd} Cook's distance.
\item \code{.fitted} Fitted or predicted value.
\item \code{.fittedlwrconf} Lower bound of the confidence interval, if conf_int = TRUE
\item \code{.fitteduprconf} Upper bound of the confidence interval, if conf_int = TRUE
\item \code{.fittedlwrpred} Lower bound of the prediction interval, if pred_int = TRUE
\item \code{.fitteduprpred} Upper bound of the prediction interval, if pred_int = TRUE
\item \code{.hat} Diagonal of the hat matrix.
\item \code{.resid} The chosen residual.
\item \code{.resfit} The chosen residual of the fitted object.
\item \code{.se.fit} Standard errors of fitted values, if se_fit = TRUE.
\item \code{.gencooksd} Generalized Cook's distance.
\item \code{.lwrenv} Lower bound of the simulated envelope, if the fitted \code{mixpoissonreg} object, was fitted with envelopes > 0.
\item \code{.mdnenv} Median of the simulated envelope, if the fitted \code{mixpoissonreg} object, was fitted with envelopes > 0.
\item \code{.uprenv} Upper bound of the simulated envelope, if the fitted \code{mixpoissonreg} object, was fitted with envelopes > 0.
}
}
\description{
Augment accepts a model object and a dataset and adds information about each observation in the dataset. It includes
predicted values in the \code{.fitted} column, residuals in the \code{.resid} column, and standard errors for the fitted values in a \code{.se.fit} column, if
the type of prediction is 'link'. New columns always begin with a . prefix to avoid overwriting columns in the original dataset.
}
\seealso{
\code{\link{glance.mixpoissonreg}}, \code{\link{tidy.mixpoissonreg}}, \code{\link{tidy_local_influence.mixpoissonreg}},
\code{\link{autoplot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
}
