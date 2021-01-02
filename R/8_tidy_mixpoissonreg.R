#' @importFrom generics augment tidy glance
#' @importFrom ggplot2 ggplot geom_text geom_hline geom_abline geom_linerange geom_point aes autoplot
#' @importFrom tibble tibble as_tibble add_column
#' @importFrom magrittr `%>%`
#' @importFrom dplyr arrange desc bind_rows mutate filter select left_join
#' @importFrom gridExtra grid.arrange
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang arg_match
#' @importFrom rlang :=
#' @importFrom methods new
#' @importClassesFrom ggfortify ggmultiplot
#' @export autoplot
#' @export augment
#' @export tidy
#' @export glance

#############################################################################################
#' @name augment.mixpoissonreg
#' @title Augment data with information from a \code{mixpoissonreg} object
#' @aliases augment augment.mixpoissonreg
#' @description Augment accepts a model object and a dataset and adds information about each observation in the dataset. It includes
#' predicted values in the \code{.fitted} column, residuals in the \code{.resid} column, and standard errors for the fitted values in a \code{.se.fit} column, if
#' the type of prediction is 'link'. New columns always begin with a . prefix to avoid overwriting columns in the original dataset.
#' @param x A \code{mixpoissonreg} object.
#' @param data A \code{\link[base:data.frame]{base::data.frame}} or \code{\link[tibble:tibble]{tibble::tibble()}} containing the original data that was used to produce the object x.
#' @param newdata A \code{\link[base:data.frame]{base::data.frame}} or \code{\link[tibble:tibble]{tibble::tibble()}} containing all the original predictors used to create x.
#' Defaults to \code{NULL}, indicating that nothing has been passed to \code{newdata}. If \code{newdata} is specified, the data argument will be ignored.
#' @param type.predict Type of prediction. The options are 'response', 'link', 'precision' and 'variance'. The default is "response".
#' @param type.residuals Type of residuals. The options are 'pearson' and 'score'. The default is 'pearson'.
#' @param se_fit Logical indicating whether or not a .se.fit column should be added to the augmented output. If TRUE, it only returns a non-NA value if type of prediction is 'link'.
#' @param conf_int Logical indicating whether or not confidence intervals for the fitted variable with type chosen from type.predict should be built. The available type options are
#' 'response' and 'link'.
#' @param pred_int Logical indicating whether or not prediction intervals for future observations should be built. It only works with type.predict = 'response'. The arguments
#' \code{level}, \code{nsim_pred}, \code{nsim_pred_y} are passed through additional arguments '...'. Notice that this can be computationally intensive.
#' @param ... Additional arguments. Possible additional arguments are \code{level}, \code{nsim_pred}, \code{nsim_pred_y}, that are passed to \code{predict.mixpoissonreg} function.
#'
#' @return A \code{\link[tibble:tibble]{tibble::tibble()}} with columns:
#' \itemize{
#'   \item \code{.cooksd} Cook's distance.
#'   \item \code{.fitted} Fitted or predicted value.
#'   \item \code{.fittedlwrconf} Lower bound of the confidence interval, if conf_int = TRUE
#'   \item \code{.fitteduprconf} Upper bound of the confidence interval, if conf_int = TRUE
#'   \item \code{.fittedlwrpred} Lower bound of the prediction interval, if pred_int = TRUE
#'   \item \code{.fitteduprpred} Upper bound of the prediction interval, if pred_int = TRUE
#'   \item \code{.hat} Diagonal of the hat matrix.
#'   \item \code{.resid} The chosen residual.
#'   \item \code{.resfit} The chosen residual of the fitted object.
#'   \item \code{.se.fit} Standard errors of fitted values, if se_fit = TRUE.
#'   \item \code{.gencooksd} Generalized Cook's distance.
#'   \item \code{.lwrenv} Lower bound of the simulated envelope, if the fitted \code{mixpoissonreg} object, was fitted with envelopes > 0.
#'   \item \code{.mdnenv} Median of the simulated envelope, if the fitted \code{mixpoissonreg} object, was fitted with envelopes > 0.
#'   \item \code{.uprenv} Upper bound of the simulated envelope, if the fitted \code{mixpoissonreg} object, was fitted with envelopes > 0.
#'   }
#'
#' @seealso \code{\link{glance.mixpoissonreg}}, \code{\link{tidy.mixpoissonreg}}, \code{\link{tidy_local_influence.mixpoissonreg}},
#' \code{\link{autoplot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
#' @export
augment.mixpoissonreg <- function(x, data = stats::model.frame(x), newdata = NULL, type.predict = c("response", "link", "precision", "variance"),
                                  type.residuals = c("pearson", "score"), se_fit = FALSE, conf_int = TRUE, pred_int = FALSE, ...) {
  .index <- .resid <- .resfit <-  NULL
  type.predict <- rlang::arg_match(type.predict)
  type.residuals <- rlang::arg_match(type.residuals)

  if(is.null(newdata)){
    newdata = stats::model.frame(x)
    res <- stats::residuals(x, type = type.residuals)
  } else{
    res <- NULL
  }

  df <- as_tibble(newdata)

  pred <- stats::predict(x, newdata = newdata, type = type.predict, se.fit = se_fit)

  if(se_fit){
    df$.fitted <- pred$fit %>% unname()
    df$.se.fit <- pred$se.fit %>% unname()
  } else{
    df$.fitted <- pred %>% unname()
  }


  if(conf_int){
    conf_int <- stats::predict(x, newdata = newdata, type = type.predict, interval = "confidence")
    df$.fittedlwrconf <- conf_int[, "lwr"] %>% unname()
    df$.fitteduprconf <- conf_int[, "upr"] %>% unname()
  }

  if(pred_int){
    pred_int <- stats::predict(x, newdata = newdata, type = "response", interval = "prediction", ...)
    df$.fittedlwrpred <- pred_int[, "lwr"] %>% unname()
    df$.fitteduprpred <-pred_int[, "upr"] %>% unname()
  }

  if(!is.null(res)){
    df$.resid <- res
    df$.resfit <- stats::residuals(x, type = x$residualname)
    df$.hat <- stats::hatvalues(x, parameters = "mean") %>% unname()
    df$.cooksd <- stats::cooks.distance(x, type = "CD", hat = "mean") %>% unname()
    df$.gencooksd <- stats::cooks.distance(x, type = "GCD", hat = "mean") %>% unname()

    env <- x$envelope

    if (!is.null(env)) {
      df$.index <- 1:nrow(df)
      temp_tbl <- tibble(.resfit = df$.resfit, .index = df$.index)
      temp_tbl <- temp_tbl %>% dplyr::arrange(.resfit)
      temp_tbl$.lwrenv <- env[3,]
      temp_tbl$.mdnenv <- env[2,]
      temp_tbl$.uprenv = env[1,]
      df <- dplyr::left_join(df, temp_tbl, by = c(".index", ".resfit")) %>% select(-.index)
    }
  }

  df
}

#############################################################################################
#' @name glance.mixpoissonreg
#' @title Glance at a \code{mixpoissonreg} object
#' @aliases glance glance.mixpoissonreg
#' @description Glance accepts a \code{mixpoissonreg} object and returns a
#' \code{\link[tibble:tibble]{tibble::tibble()}} with exactly one row of model summaries.
#' The summaries are Efron's pseudo-\eqn{R^2}, degrees of freedom, AIC, BIC, log-likelihood,
#' the type of model used in the fit ('NB' or 'PIG'), the total number of observations and the estimation method.
#' @param x A \code{mixpoissonreg} object.
#' @param ... Additional arguments. Currently not used.
#' @return A \code{\link[tibble:tibble]{tibble::tibble()}} with exactly one row and columns:
#' \itemize{
#'   \item \code{efron.pseudo.r2} Efron's pseudo-\eqn{R^2}, that is, the squared correlation between the fitted values and the response values.
#'   \item \code{df.null} Degrees of freedom used by the null model.
#'   \item \code{logLik} The log-likelihood of the model.
#'   \item \code{AIC} Akaike's Information Criterion for the model.
#'   \item \code{BIC} Bayesian Information Criterion for the model.
#'   \item \code{df.residual} Residual degrees of freedom.
#'   \item \code{nobs} Number of observations used.
#'   \item \code{model.type} Type of model fitted, "NB" or "PIG".
#'   \item \code{est.method} The estimation method of the fitted model, "EM" or "ML".
#'   }
#' @seealso \code{\link{augment.mixpoissonreg}}, \code{\link{tidy.mixpoissonreg}}, \code{\link{tidy_local_influence.mixpoissonreg}},
#' \code{\link{autoplot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
#' @export

glance.mixpoissonreg <- function(x, ...){
  tibble(efron.pseudo.r2 = as.numeric(x$efron.pseudo.r2), df.null = x$df.null,
                   logLik = as.numeric(stats::logLik(x)), AIC = stats::AIC(x),
                   BIC = stats::BIC(x), df.residual = stats::df.residual(x),
                   nobs = stats::nobs(x), model.type = x$modeltype, est.method = x$estimation_method)
}

#############################################################################################
#' @name tidy.mixpoissonreg
#' @title Tidy a \code{mixpoissonreg} object
#' @aliases tidy tidy.mixpoissonreg
#' @description Tidy returns a \code{\link[tibble:tibble]{tibble::tibble()}} containing
#' informations on the coefficients of the model, such as the estimated parameters,
#' standard errors, z-statistics and p-values. Additionally, it may return confidence
#' intervals for the model parameters.
#' @param x A \code{mixpoissonreg} object.
#' @param conf.int Logical indicating whether or not to include a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level The confidence level to use for the confidence interval if conf.int = TRUE.
#' Must be strictly greater than 0 and less than 1. Defaults to 0.95, which corresponds to a 95 percent confidence interval.
#' @param ... Additional arguments. Currently not used.
#' @seealso \code{\link{glance.mixpoissonreg}}, \code{\link{augment.mixpoissonreg}}, \code{\link{tidy_local_influence.mixpoissonreg}},
#' \code{\link{autoplot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
#' @export

tidy.mixpoissonreg <- function(x, conf.int = FALSE, conf.level = 0.95, ...){
  join_term <- NULL
  retmean <- as_tibble(summary(x)$coefficients$mean, rownames = "term")
  colnames(retmean) <- c("term", "estimate", "std.error", "statistic",
                     "p.value")
  retmean <- retmean %>% add_column(component = "mean", .before = "term")
  retprec <- as_tibble(summary(x)$coefficients$precision, rownames = "term")
  colnames(retprec) <- c("term", "estimate", "std.error", "statistic",
                     "p.value")
  retprec <- retprec %>% add_column(component = "precision", .before = "term")

  ret <- dplyr::bind_rows(retmean,retprec)

  if (conf.int) {
    ret$join_term <- names(stats::coef(x, parameters = "all"))
    ci <- as_tibble(stats::confint(x, level = conf.level), rownames = "term")
    names(ci) <- c("join_term", "conf.low", "conf.high")
    ret <- dplyr::left_join(ret, ci, by = "join_term") %>% select(-join_term)
  }
  ret
}


#############################################################################################
#' @name autoplot.mixpoissonreg
#' @title Autoplot Method for \code{mixpoissonreg} Objects
#' @aliases autoplot autoplot.mixpoissonreg
#' @description This function provides \pkg{ggplot2}-based counterparts to the plots produced by \code{\link{plot.mixpoissonreg}}.
#' Currently there are six plots available. They contain residual analysis and global influence diagnostics. The plots are selectable by
#' the \code{which} argument. The plots are: Residuals vs. obs. numbers; Normal Q-Q plots, which may contain simulated envelopes, if the fitted object
#' has simulated envelopes; Cook's distances vs. obs. numbers; Generalized Cook's distances vs. obs. numbers; Cook's distances vs. Generalized Cook's distances;
#' Response variables vs. fitted means. By default, the first two plots and the last two plots are provided.
#'
#' If both ncol and nrow are \code{NULL}, the plots will be placed one at a time. To place multiple plots, set the values for \code{nrow} or \code{ncol}.
#' @param object A \code{mixpoissonreg} object.
#' @param which a list or vector indicating which plots should be displayed. 	If a subset of the plots is required, specify a subset of the numbers 1:6,
#' see caption below for the different kinds. In
#' plot number 2, 'Normal Q-Q', if the \code{mixpoissonreg} object was fitted with envelopes, a quantile-quantile plot with simulated envelopes will be displayed.
#' @param title titles to appear above the plots; character vector or list of valid graphics annotations. Can be set to "" to suppress all captions.
#' @param sub.caption common title-above the figures if there are more than one. If NULL, as by default, a possible abbreviated version of \code{deparse(x$call)} is used.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param label.label vector of labels. If \code{NULL}, rownames will be used as labels.
#' @param env_alpha alpha of the envelope region (when the fitted model has envelopes)
#' @param env_fill the colour of the filling in the envelopes.
#' @param gpar_sub.caption list of gpar parameters to be used as common title in the case of multiple plots. The title will be given in sub.caption argument. See
#' the help of \code{\link[grid]{gpar}} function from the \pkg{grid} package for all the available options.
#' @param include.modeltype logical. Indicates whether the model type ('NB' or 'PIG') should be displayed on the captions.
#' @param include.residualtype local. Indicates whether the name of the residual ('Pearson' or 'Score') should be displayed on the caption of plot 1 (Residuals vs. Index).
#' @param label.repel Logical flag indicating whether to use \pkg{ggrepel} to place the labels.
#' @param nrow Number of facet/subplot rows. If both \code{nrow} and \code{ncol} are \code{NULL}, the plots will be placed one at a time. For multiple plots, set values for \code{nrow}
#' or \code{ncol}.
#' @param ncol Number of facet/subplot columns. If both \code{nrow} and \code{ncol} are \code{NULL}, the plots will be placed one at a time. For multiple plots, set values for \code{nrow}
#' or \code{ncol}.
#' @param qqline logical; if \code{TRUE} and the fit does *not* contain simulated
#' envelopes, a qqline passing through the first and third quartiles of a standard normal distribution will be added to the normal Q-Q plot.
#' @param colour line colour.
#' @param size	point size.
#' @param linetype	line type.
#' @param alpha	alpha of the plot.
#' @param fill fill colour.
#' @param shape	point shape.
#' @param label Logical value whether to display labels.
#' @param label.colour Colour for text labels.
#' @param label.alpha	Alpha for text labels.
#' @param label.size Size for text labels.
#' @param label.angle	Angle for text labels.
#' @param label.family Font family for text labels.
#' @param label.fontface Fontface for text labels.
#' @param label.lineheight Lineheight for text labels.
#' @param label.hjust	Horizontal adjustment for text labels.
#' @param label.vjust	Vertical adjustment for text labels.
#' @param label.n	Number of points to be laeled in each plot, starting with the most extreme.
#' @param ad.colour	Line colour for additional lines.
#' @param ad.linetype	Line type for additional lines.
#' @param ad.size	Fill colour for additional lines.
#' @param ... other arguments passed to methods.
#' @details Based on \code{autoplot.lm} from the excellent \pkg{ggfortify} package, \href{https://github.com/sinhrks/ggfortify/}{ggfortify}.
#'
#' sub.caption - by default the function call - is shown as a subtitle (under the x-axis title) on each plot when plots are on separate pages, or as a subtitle
#' in the outer margin when there are multiple plots per page.
#'
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' autoplot(daysabs_fit, which = 1:6)
#'
#' autoplot(daysabs_fit, nrow = 2)
#'
#' daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance, envelope = 100)
#' autoplot(daysabs_fit_ml, which = 2)
#' }
#' @export



autoplot.mixpoissonreg <- function(object, which = c(1,2,5,6), title = list("Residuals vs Obs. number",
                                                                            "Normal Q-Q",
                                                                            "Cook's distance",
                                                                            "Generalized Cook's distance",
                                                                            "Cook's dist vs Generalized Cook's dist",
                                                                            "Response vs Fitted means"),
                                   label.repel = TRUE,
                                   nrow = NULL, ncol = NULL,
                                   qqline = TRUE, ask = prod(graphics::par("mfcol")) <
                                     length(which) && grDevices::dev.interactive(), include.modeltype = TRUE,
                                   include.residualtype = FALSE, sub.caption = NULL,
                                   env_alpha = 0.5, env_fill = "grey70", gpar_sub.caption = list(fontface = "bold"),
                                   colour = "#444444", size = NULL, linetype = NULL, alpha = NULL, fill = NULL,
                                   shape = NULL, label = TRUE, label.label = NULL, label.colour = "#000000",
                                   label.alpha = NULL, label.size = NULL, label.angle = NULL,
                                   label.family = NULL, label.fontface = NULL, label.lineheight = NULL,
                                   label.hjust = NULL, label.vjust = NULL,
                                   label.n = 3, ad.colour = "#888888", ad.linetype = "dashed", ad.size = 0.2, ...){
  p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- NULL
  .resid <- .cooksd <- .gencooksd <- .obs <- .fitted <-.lwrenv <- .uprenv <- .mdnenv <- NULL

  if (is.null(sub.caption)) {
    cal <- object$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2L] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1L], "c")
    abbr <- length(cc) > 1 || nc > 75
    sub.caption <- if (abbr)
      paste(substr(cc[1L], 1L, min(75L, nc)), "...")
    else cc[1L]
  }

  plot.data <- augment(object, type.residuals = object$residualname, type.predict = "response")
  n <- stats::nobs(object)
  plot.data$.index <- 1:n
  if(is.null(label.label)){
    plot.data$.label <- rownames(plot.data)
  } else{
    plot.data$.label <- as.vector(label.label)
  }

  if(is.null(object$y)){
    y = object$residuals + object$fitted.values
  } else{
    y = object$y
  }

  if(2 %in% which){
    ylim <- range(plot.data$.resid, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    qx <- stats::qqnorm(plot.data$.resid, ylim = ylim,
                        plot.it = FALSE)$x
    plot.data$.qqx <- qx

    qprobs <- c(0.25, 0.75)
    qy <- stats::quantile(plot.data$.resid, probs = qprobs,
                          names = FALSE, type = 7, na.rm = TRUE)
    qx <- stats::qnorm(qprobs)
    slope <- diff(qy)/diff(qx)
    int <- qy[1L] - slope * qx[1L]
  }

  plot.data$.obs <- as.vector(y)
  
  # Internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  plot_label <- function (p, data, x = NULL, y = NULL, label = TRUE, label.label = "rownames",
                          label.colour = NULL, label.alpha = NULL, label.size = NULL,
                          label.angle = NULL, label.family = NULL, label.fontface = NULL,
                          label.lineheight = NULL, label.hjust = NULL, label.vjust = NULL,
                          label.repel = FALSE, label.show.legend = NA)
  {
    if (!is.data.frame(data)) {
      stop(paste0("Unsupported class: ", class(data)))
    }
    if (!missing(label.colour) && !is.null(label.colour) && missing(label)) {
      label <- TRUE
    }
    if (label || label.repel) {
      if (is.null(label.colour)) {
        label.colour <- "#000000"
      }
      if (label.repel) {
        textfunc <- ggrepel::geom_text_repel
      }
      else {
        textfunc <- ggplot2::geom_text
      }
      p <- p + geom_factory(textfunc, data, x = x, y = y, label = ".label",
                            colour = label.colour, alpha = label.alpha, size = label.size,
                            angle = label.angle, family = label.family, fontface = label.fontface,
                            lineheight = label.lineheight, hjust = label.hjust,
                            vjust = label.vjust, show.legend = label.show.legend)
    }
    p
  }

  # Internal function from ggfortify package (https://github.com/sinhrks/ggfortify)
  
  flatten <- function (df)
  {
    ismatrix <- vapply(df, is.matrix, logical(1))
    if (any(ismatrix)) {
      return(data.frame(c(df[!ismatrix], do.call(data.frame,
                                                 df[ismatrix])), stringsAsFactors = FALSE))
    }
    else {
      return(df)
    }
  }

  # Internal function from ggfortify package (https://github.com/sinhrks/ggfortify)
  
  geom_factory <- function (geomfunc, data = NULL, ...)
  {
    mapping <- list()
    option <- list()
    columns <- colnames(data)
    for (key in names(list(...))) {
      value <- list(...)[[key]]
      if (is.null(value)) {
      }
      else if (value %in% columns) {
        mapping[[key]] <- value
      }
      else {
        option[[key]] <- value
      }
    }
    if (!is.null(data)) {
      option[["data"]] <- data
    }
    option[["mapping"]] <- do.call(ggplot2::aes_string, mapping)
    return(do.call(geomfunc, option))
  }

  if (is.logical(shape) && !shape) {
    if (missing(label)) {
      label <- TRUE
    }
    if (missing(label.n)) {
      label.n <- nrow(plot.data)
    }
  }

  dev_ask <- is.null(nrow) & is.null(ncol)

  if(dev_ask){

    if (ask) {
      oask <- grDevices::devAskNewPage(TRUE)
      on.exit(grDevices::devAskNewPage(oask))
    }
  }

  plot.data <- flatten(plot.data)
  if(any(c(1,2) %in% which)){
    residualname <- paste0(toupper(substring(object$residualname, 1, 1)), substring(object$residualname, 2))
  }

  if (label.n > 0L) {
    if (any(c(1,2,6) %in% which)) {
      r.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.resid)))
      r.data <- utils::head(r.data, label.n)
    }
    if (3 %in% which) {
      cd.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.cooksd)))
      cd.data <- utils::head(cd.data, label.n)
    }
    if (4 %in% which){
      gcd.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.gencooksd)))
      gcd.data <- utils::head(gcd.data, label.n)
    }
    if(5 %in% which){
      cdgcd.data <- dplyr::arrange(plot.data, dplyr::desc((.gencooksd/sum(.gencooksd))^2 + (.cooksd/sum(.cooksd))^2))
      cdgcd.data <- utils::head(cdgcd.data, label.n)
    }
    if(6 %in% which){
      respobs.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.obs - .fitted)))
      respobs.data <- utils::head(respobs.data, label.n)
    }
  }
  
  # Internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  .decorate.label <- function(p, data) {
    if (label & label.n > 0) {
      p <- plot_label(p = p, data = data, label = label,
                      label.label = ".label", label.colour = label.colour,
                      label.alpha = label.alpha, label.size = label.size,
                      label.angle = label.angle, label.family = label.family,
                      label.fontface = label.fontface, label.lineheight = label.lineheight,
                      label.hjust = label.hjust, label.vjust = label.vjust,
                      label.repel = label.repel)
    }
    p
  }
  
  # Internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  .decorate.plot <- function(p, xlab = NULL, ylab = NULL, title = NULL) {
    p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
  }
  
  # Plots

  if (1 %in% which) {
    if(include.residualtype){
      title[[1]] = paste(residualname, title[[1]])
    }

    ylab <- paste0(residualname, " residuals")
    
    if(include.modeltype){
      title[[1]] <- paste0(title[[1]], " - ", object$modeltype, " Regression")
    }
    
    t1 <- title[[1]]
    
    mapping <- ggplot2::aes_string(x = ".index", y = ".resid")
    p1 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p1 <- p1 + geom_factory(geom_point, plot.data, colour = colour,
                              size = size, linetype = linetype, alpha = alpha,
                              fill = fill, shape = shape)
    }
    p1 <- p1 + ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype,
                                   size = ad.size, colour = ad.colour)
    p1 <- .decorate.label(p1, r.data)

    if(sub.caption == "" | !dev_ask){
      xlab = "Obs. number"
    } else {
      xlab = paste0("Obs. number\n",sub.caption)
    }
    p1 <- .decorate.plot(p1, xlab = xlab, ylab = ylab,
                         title = t1)
    if(dev_ask){
      print(p1)
    }
  }

  if(2 %in% which){
    env <- object$envelope
    if (!is.null(env)) {
      title[[2]] <- paste0(title[[2]]," with simulated envelopes")
    }
    
    if(include.modeltype){
      title[[2]] <- paste0(title[[2]], " - ", object$modeltype, " Regression")
    }
    
    t2 <- title[[2]]

    mapping <- ggplot2::aes_string(x = ".qqx", y = ".resid")
    p2 <- ggplot2::ggplot(data = plot.data, mapping = mapping)

    if(qqline & is.null(env)){
      p2 <- p2 + ggplot2::geom_abline(intercept = int, slope = slope,
                                      linetype = ad.linetype, size = ad.size, colour = ad.colour)
    }

    if(!is.null(env)){
      p2 <- p2 + ggplot2::geom_ribbon(aes(ymin=.lwrenv, ymax=.uprenv), alpha=env_alpha, fill = env_fill)
      p2 <- p2 + ggplot2::geom_line(aes(y = .mdnenv), linetype = ad.linetype, size = ad.size, colour = ad.colour)
    }

    if (!is.logical(shape) || shape) {
      p2 <- p2 + geom_factory(geom_point, plot.data, colour = colour,
                              size = size, linetype = linetype, alpha = alpha,
                              fill = fill, shape = shape)
    }

    ylab <- paste0(residualname, " residuals")

    if(sub.caption == "" | !dev_ask){
      xlab = "Theoretical Quantiles"
    } else {
      xlab = paste0("Theoretical Quantiles\n", sub.caption)
    }

    p2 <- .decorate.label(p2, r.data)
    p2 <- .decorate.plot(p2, xlab = xlab,
                         ylab = ylab, title = t2)
    if(dev_ask){
      print(p2)
    }
  }

  if(3 %in% which){
    
    if(include.modeltype){
      title[[3]] <- paste0(title[[3]], " - ", object$modeltype, " Regression")
    }
    
    t3 <- title[[3]]
    
    mapping <- ggplot2::aes_string(x = ".index", y = ".cooksd",
                                   ymin = 0, ymax = ".cooksd")
    p3 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p3 <- p3 + geom_factory(geom_linerange, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p3 <- .decorate.label(p3, cd.data)

    if(sub.caption == "" | !dev_ask){
      xlab = "Obs. Number"
    } else {
      xlab = paste0("Obs. Number\n", sub.caption)
    }

    p3 <- .decorate.plot(p3, xlab = xlab, ylab = "Cook's distance",
                         title = t3)
    if(dev_ask){
      print(p3)
    }
  }

  if(4 %in% which){
    
    if(include.modeltype){
      title[[4]] <- paste0(title[[4]], " - ", object$modeltype, " Regression")
    }
    
    t4 <- title[[4]]
    
    mapping <- ggplot2::aes_string(x = ".index", y = ".gencooksd",
                                   ymin = 0, ymax = ".gencooksd")
    p4 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p4 <- p4 + geom_factory(geom_linerange, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p4 <- .decorate.label(p4, gcd.data)

    if(sub.caption == "" | !dev_ask){
      xlab = "Obs. Number"
    } else {
      xlab = paste0("Obs. Number\n", sub.caption)
    }

    p4 <- .decorate.plot(p4, xlab = xlab, ylab = "Generalized Cook's distance",
                         title = t4)
    if(dev_ask){
      print(p4)
    }
  }

  if (5 %in% which) {
    
    if(include.modeltype){
      title[[5]] <- paste0(title[[5]], " - ", object$modeltype, " Regression")
    }
    
    t5 <- title[[5]]
    
    mapping <- ggplot2::aes_string(x = ".gencooksd", y = ".cooksd")
    p5 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p5 <- p5 + geom_factory(geom_point, plot.data, colour = colour,
                              size = size, linetype = linetype, alpha = alpha,
                              fill = fill, shape = shape)
    }
    p5 <- p5 + ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype,
                                   size = ad.size, colour = ad.colour)
    p5 <- .decorate.label(p5, cdgcd.data)

    if(sub.caption == "" | !dev_ask){
      xlab = "Generalized Cook's distance"
    } else {
      xlab = paste0("Generalized Cook's distance\n", sub.caption)
    }

    p5 <- .decorate.plot(p5, xlab = xlab, ylab = "Cook's distance",
                         title = t5)
    if(dev_ask){
      print(p5)
    }
  }

  if (6 %in% which) {
    
    if(include.modeltype){
      title[[6]] <- paste0(title[[6]], " - ", object$modeltype, " Regression")
    }
    
    t6 <- title[[6]]
    
    mapping <- ggplot2::aes_string(x = ".fitted", y = ".obs")
    p6 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p6 <- p6 + geom_factory(geom_point, plot.data, colour = colour,
                              size = size, linetype = linetype, alpha = alpha,
                              fill = fill, shape = shape)
    }
    p6 <- p6 + ggplot2::geom_abline(intercept = 0L, slope = 1L, linetype = ad.linetype,
                                    size = ad.size, colour = ad.colour)
    p6 <- .decorate.label(p6, respobs.data)

    if(sub.caption == "" | !dev_ask){
      xlab = "Predicted values"
    } else {
      xlab = paste0("Predicted values\n", sub.caption)
    }

    p6 <- .decorate.plot(p6, xlab = xlab, ylab = "Response values",
                         title = t6)
    if(dev_ask){
      print(p6)
    }
  }


  if(!dev_ask){
    grDevices::devAskNewPage(ask = FALSE)
    plot.list <- list(p1, p2, p3, p4, p5, p6)[which]
    if(is.null(nrow))
      nrow <- 0
    if(is.null(ncol))
      ncol <- 0
    p <- methods::new("ggmultiplot", plots = plot.list, nrow = nrow, ncol = ncol)

    gpar_multi <- do.call(grid::gpar, gpar_sub.caption)
    title_multi <- grid::textGrob(sub.caption, gp=gpar_multi)
    gridExtra::grid.arrange(grobs = p@plots, top = title_multi)
    grDevices::devAskNewPage(ask = grDevices::dev.interactive())
  } else{
    invisible()
  }
}


#############################################################################################
#' @name tidy_local_influence.mixpoissonreg
#' @title Tidy Functions for Local Influence Diagnostics for \code{mixpoissonreg} Objects
#' @aliases local_influence_benchmarks.mixpoissonreg tidy_local_influence.mixpoissonreg
#' @description Functions to provide tidy outputs of local influence diagnostics. \code{tidy_local_influence.mixpoissonreg}
#' provides a \code{\link[tibble:tibble]{tibble::tibble()}} containing the local influence diagnostics under the chosen perturbation schemes.
#' \code{local_influence_benchmarks.mixpoissonreg} provides a \code{\link[tibble:tibble]{tibble::tibble()}} with a single row and one column
#' for each selected perturbation scheme containing influential benchmarks for each perturbation scheme.
#' @param model A \code{mixpoissonreg} model.
#' @param perturbation a list or vector of perturbation schemes to be returned. The currently available schemes are
#' "case_weights", "hidden_variable", "mean_explanatory", "precision_explanatory", "simultaneous_explanatory". See Barreto-Souza and Simas (2016) for further details.
#' @param curvature the curvature to be returned, 'conformal' for the conformal normal curvature (see Zhu and Lee, 2001 and Poon and Poon, 1999) or
#' 'normal' (see Zhu and Lee, 2001 and Cook, 1986).
#' @param direction the 'max.eigen' returns the eigenvector associated to the largest eigenvalue of the perturbation matrix. The 'canonical' considers
#' the curvatures under the canonical directions, which is known as "total local curvature" (see Lesaffre and Verbeke, 1998). For conformal
#' normal curvatures both of them coincide. The default is 'canonical'.
#' @param parameters the parameter to which the local influence will be computed. The options are 'all', 'mean' and 'precision'.
#' This argument affects the 'case_weights' and 'hidden_variable' perturbation schemes. The default is 'all'.
#' @param mean.covariates a list or vector of characters containing the mean-explanatory variables to be used in the 'mean-explanatory' and 'simultaneous-explanatory'
#' perturbation schemes. If NULL, the 'mean-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
#' mean-related covariates. The default is NULL.
#' @param precision.covariates a list or vector of characters containing the precision-explanatory variables to be used in the 'precision-explanatory'
#' and 'simultaneous-explanatory'
#' perturbation schemes. If NULL, the 'precision-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
#' precision-related covariates. The default is NULL.
#' @param ... Currently not used.
#' @references
#' DOI:10.1007/s11222-015-9601-6 (\href{https://doi.org/10.1007/s11222-015-9601-6}{Barreto-Souza and Simas; 2016})
#'
#' Cook, R. D. (1986) *Assessment of Local Influence.* Journal of the Royal Statistical Society. Series B (Methodological), Vol. 48, pp.133-169. \href{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}
#'
#' Lesaffre, E. and Verbeke, G. (1998) *Local Influence in Linear Mixed Models*. Biometrics, 54, pp. 570-582. \href{https://www.jstor.org/stable/3109764}{https://www.jstor.org/stable/3109764}
#'
#' Poon, W.-Y. and Poon, Y.S. (1999) *Conformal normal curvature and assessment of local influence.*  Journal of the Royal Statistical Society. Series B (Methodological), Vol. 61, pp.51-61. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}
#'
#' Zhu, H.-T. and Lee, S.-Y. (2001) *Local influence for incomplete data models.* Journal of the Royal Statistical Society. Series B (Methodological), Vol. 63, pp.111-126. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}
#' @seealso \code{\link{glance.mixpoissonreg}}, \code{\link{augment.mixpoissonreg}}, \code{\link{tidy.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' tidy_local_influence(daysabs_fit)
#'
#' daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance, envelope = 100)
#' tidy_local_influence(daysabs_fit_ml, perturbation = "case_weights")
#' }
#' @rdname tidy_local_influence.mixpoissonreg
#' @export
tidy_local_influence.mixpoissonreg <- function(model, perturbation = c("case_weights", "hidden_variable",
                                                     "mean_explanatory", "precision_explanatory",
                                                     "simultaneous_explanatory"), curvature = c("conformal", "normal"),
                                 direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                 mean.covariates = NULL, precision.covariates = NULL, ...){
    loc_infl <- suppressWarnings(local_influence(model, perturbation = perturbation, curvature = curvature, direction = direction,
                                parameters = parameters, mean.covariates = mean.covariates, precision.covariates = precision.covariates))

    tidy_loc_infl <- tibble(.rows = stats::nobs(model))

    for(pert in perturbation){
      tidy_loc_infl = tidy_loc_infl %>% add_column(!!pert := loc_infl[[pert]])
    }
    tidy_loc_infl
}

#' @rdname tidy_local_influence.mixpoissonreg
#' @export
local_influence_benchmarks.mixpoissonreg <- function(model, perturbation = c("case_weights", "hidden_variable",
                                                           "mean_explanatory", "precision_explanatory",
                                                           "simultaneous_explanatory"), curvature = c("conformal", "normal"),
                                       direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                       mean.covariates = NULL, precision.covariates = NULL, ...){
  loc_infl <- local_influence(model, perturbation = perturbation, curvature = curvature, direction = direction,
                              parameters = parameters, mean.covariates = mean.covariates, precision.covariates = precision.covariates)
  benchmarks <- c()
  for(pert in perturbation){
    benchmarks <- c(benchmarks, attr(loc_infl[[pert]], "benchmark"))
  }
  benchmarks <- matrix(benchmarks, nrow = 1)
  colnames(benchmarks) <- perturbation
  benchmarks_tbl <- as_tibble(benchmarks)
  benchmarks_tbl
}

#############################################################################################
#' @name local_influence_autoplot.mixpoissonreg
#' @title Local Influence Autoplots for \code{mixpoissonreg} Objects
#' @description Function to provide customizable ggplot2-based plots of local influence diagnostics.
#' @param model A \code{mixpoissonreg} model.
#' @param which a list or vector indicating which plots should be displayed. 	If a subset of the plots is required, specify a subset of the numbers 1:5, see caption below (and the 'Details') for the different kinds.
#' @param title titles to appear above the plots; character vector or list of valid graphics annotations. Can be set to "" to suppress all titles.
#' @param sub.caption	common title-above the figures if there are more than one. If NULL, as by default, a possible abbreviated version of deparse(x$call) is used.
#' @param curvature the curvature to be returned, 'conformal' for the conformal normal curvature (see Zhu and Lee, 2001 and Poon and Poon, 1999) or
#' 'normal' (see Zhu and Lee, 2001 and Cook, 1986).
#' @param direction the 'max.eigen' returns the eigenvector associated to the largest eigenvalue of the perturbation matrix. The 'canonical' considers
#' the curvatures under the canonical directions, which is known as "total local curvature" (see Lesaffre and Verbeke, 1998). For conformal
#' normal curvatures both of them coincide. The default is 'canonical'.
#' @param parameters the parameter to which the local influence will be computed. The options are 'all', 'mean' and 'precision'.
#' This argument affects the 'case_weights' and 'hidden_variable' perturbation schemes. The default is 'all'.
#' @param mean.covariates a list or vector of characters containing the mean-explanatory variables to be used in the 'mean-explanatory' and 'simultaneous-explanatory'
#' perturbation schemes. If NULL, the 'mean-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
#' mean-related covariates. The default is NULL.
#' @param precision.covariates a list or vector of characters containing the precision-explanatory variables to be used in the 'precision-explanatory'
#' and 'simultaneous-explanatory'
#' perturbation schemes. If NULL, the 'precision-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
#' precision-related covariates. The default is NULL.
#' @param detect.influential logical. Indicates whether the benchmark should be used to detect influential observations and identify them on the plot. If there is no benchmark available,
#' the top 'n.influential' observations will be identified in the plot by their indexes.
#' @param n.influential interger. The maximum number of influential observations to be identified on the plot.
#' @param draw.benchmark logical. Indicates whether a horizontal line identifying the benchmark should be drawn.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param label.label vector of labels. If \code{NULL}, rownames will be used as labels.
#' @param gpar_sub.caption list of gpar parameters to be used as common title in the case of multiple plots. The title will be given in sub.caption argument. See
#' the help of \code{\link[grid]{gpar}} function from the \pkg{grid} package for all the available options.
#' @param include.modeltype logical. Indicates whether the model type ('NB' or 'PIG') should be displayed on the captions.
#' @param label.repel Logical flag indicating whether to use \pkg{ggrepel} to place the labels.
#' @param nrow Number of facet/subplot rows. If both \code{nrow} and \code{ncol} are \code{NULL}, the plots will be placed one at a time. For multiple plots, set values for \code{nrow}
#' or \code{ncol}.
#' @param ncol Number of facet/subplot columns. If both \code{nrow} and \code{ncol} are \code{NULL}, the plots will be placed one at a time. For multiple plots, set values for \code{nrow}
#' or \code{ncol}.
#' @param colour line colour.
#' @param size	point size.
#' @param linetype	line type.
#' @param alpha	alpha of the plot.
#' @param fill fill colour.
#' @param shape	point shape.
#' @param label Logical value whether to display labels.
#' @param label.colour Colour for text labels.
#' @param label.alpha	Alpha for text labels.
#' @param label.size Size for text labels.
#' @param label.angle	Angle for text labels.
#' @param label.family Font family for text labels.
#' @param label.fontface Fontface for text labels.
#' @param label.lineheight Lineheight for text labels.
#' @param label.hjust	Horizontal adjustment for text labels.
#' @param label.vjust	Vertical adjustment for text labels.
#' @param ad.colour	Line colour for additional lines.
#' @param ad.linetype	Line type for additional lines.
#' @param ad.size	Fill colour for additional lines.
#' @param ... Currently not used.
#' @references
#' DOI:10.1007/s11222-015-9601-6 (\href{https://doi.org/10.1007/s11222-015-9601-6}{Barreto-Souza and Simas; 2016})
#'
#' Cook, R. D. (1986) *Assessment of Local Influence.* Journal of the Royal Statistical Society. Series B (Methodological), Vol. 48, pp.133-169. \href{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}
#'
#' Lesaffre, E. and Verbeke, G. (1998) *Local Influence in Linear Mixed Models*. Biometrics, 54, pp. 570-582. \href{https://www.jstor.org/stable/3109764}{https://www.jstor.org/stable/3109764}
#'
#' Poon, W.-Y. and Poon, Y.S. (2002) *Conformal normal curvature and assessment of local influence.*  Journal of the Royal Statistical Society. Series B (Methodological), Vol. 61, pp.51-61. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}
#'
#' Zhu, H.-T. and Lee, S.-Y. (2001) *Local influence for incomplete data models.* Journal of the Royal Statistical Society. Series B (Methodological), Vol. 63, pp.111-126. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}
#' @seealso \code{\link{glance.mixpoissonreg}}, \code{\link{augment.mixpoissonreg}}, \code{\link{tidy.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' local_influence_autoplot(daysabs_fit)
#'
#' local_influence_autoplot(daysabs_fit, nrow = 2)
#'
#' daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance, envelope = 100)
#' local_influence_autoplot(daysabs_fit_ml, which = 2)
#' }
#' @export
local_influence_autoplot.mixpoissonreg <- function(model, which = c(1,2,3,4), title = list("Case Weights Perturbation",
                                                                                           "Hidden Variable Perturbation",
                                                                                           "Mean Explanatory Perturbation",
                                                                                           "Precision Explanatory Perturbation",
                                                                                           "Simultaneous Explanatory Perturbation"),
                                                   curvature = c("conformal", "normal"),
                                                   direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                                   mean.covariates = NULL, precision.covariates = NULL,
                                                   label.repel = TRUE,
                                                   nrow = NULL, ncol = NULL,
                                                   ask = prod(graphics::par("mfcol")) <
                                                     length(which) && grDevices::dev.interactive(), include.modeltype = TRUE,
                                                   sub.caption = NULL,
                                                   gpar_sub.caption = list(fontface = "bold"), detect.influential = TRUE, n.influential = 5,
                                                   draw.benchmark = FALSE,
                                                   colour = "#444444", size = NULL, linetype = NULL, alpha = NULL, fill = NULL,
                                                   shape = NULL, label = TRUE, label.label = NULL, label.colour = "#000000",
                                                   label.alpha = NULL, label.size = NULL, label.angle = NULL,
                                                   label.family = NULL, label.fontface = NULL, label.lineheight = NULL,
                                                   label.hjust = NULL, label.vjust = NULL,
                                                   ad.colour = "#888888", ad.linetype = "dashed", ad.size = 0.2, ...){
  p <- list()
  tmp <- NULL
  p[[1]] <- p[[2]] <- p[[3]] <- p[[4]] <- p[[5]] <- NULL

  n_beta = length(model$coefficients$mean)
  n_alpha = length(model$coefficients$precision)

  if(model$intercept[1] & n_beta == 1 & any(c(3,5)%in% which)){
    warning("Removing mean explanatory and simultaneous explanatory perturbations since
          there is only the intercept for the mean.")
    which <- setdiff(which, c(3,5))
  }

  if(model$intercept[2] & n_alpha == 1 & any(c(4,5)%in% which)){
    warning("Removing precision explanatory and simultaneous explanatory perturbations since
          there is only the intercept for the precision parameter.")
    which <- setdiff(which, c(4,5))
  }

  if (is.null(sub.caption)) {
    cal <- model$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2L] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1L], "c")
    abbr <- length(cc) > 1 || nc > 75
    sub.caption <- if (abbr)
      paste(substr(cc[1L], 1L, min(75L, nc)), "...")
    else cc[1L]
  }

  direction <- rlang::arg_match(direction)
  curvature <- rlang::arg_match(curvature)

  pert <- c("case_weights", "hidden_variable",
            "mean_explanatory", "precision_explanatory",
            "simultaneous_explanatory")

  plot.data <- tidy_local_influence(model, perturbation = pert[which], curvature = curvature,
                                    direction = direction, parameters = parameters,
                                    mean.covariates = mean.covariates,
                                    precision.covariates = precision.covariates)
  n <- stats::nobs(model)
  plot.data$.index <- 1:n

  if(is.null(label.label)){
    plot.data$.label <- rownames(plot.data)
  } else{
    plot.data$.label <- as.vector(label.label)
  }

  # Based on internal function from ggfortify package (https://github.com/sinhrks/ggfortify)
  
  plot_label_influential <- function (p, data, x = NULL, y = NULL, label = TRUE, label.label = "rownames",
                                      label.colour = NULL, label.alpha = NULL, label.size = NULL,
                                      label.angle = NULL, label.family = NULL, label.fontface = NULL,
                                      label.lineheight = NULL, label.hjust = NULL, label.vjust = NULL,
                                      label.repel = FALSE, label.show.legend = NA)
  {
    if (!is.data.frame(data)) {
      stop(paste0("Unsupported class: ", class(data)))
    }
    if (!missing(label.colour) && !is.null(label.colour) && missing(label)) {
      label <- TRUE
    }
    if (label || label.repel) {
      if (is.null(label.colour)) {
        label.colour <- "#000000"
      }
      if (label.repel) {
        textfunc <- ggrepel::geom_text_repel
      }
      else {
        textfunc <- ggplot2::geom_text
      }
      p <- p + geom_factory_influential(textfunc, data, x = x, y = y, label = ".label",
                                        colour = label.colour, alpha = label.alpha, size = label.size,
                                        angle = label.angle, family = label.family, fontface = label.fontface,
                                        lineheight = label.lineheight, hjust = label.hjust,
                                        vjust = label.vjust, show.legend = label.show.legend)
    }
    p
  }
  
  # Internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  flatten <- function (df)
  {
    ismatrix <- vapply(df, is.matrix, logical(1))
    if (any(ismatrix)) {
      return(data.frame(c(df[!ismatrix], do.call(data.frame,
                                                 df[ismatrix])), stringsAsFactors = FALSE))
    }
    else {
      return(df)
    }
  }
  
  # Based on internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  geom_factory_influential <- function (geomfunc, data = NULL, ...)
  {
    mapping <- list()
    option <- list()
    columns <- colnames(data)
    for (key in names(list(...))) {
      value <- list(...)[[key]]
      if (is.null(value)) {
      }
      else if (value %in% columns) {
        mapping[[key]] <- value
      }
      else {
        option[[key]] <- value
      }
    }
    if (!is.null(data)) {
      option[["data"]] <- data
    }
    option[["mapping"]] <- do.call(ggplot2::aes_string, mapping)
    return(do.call(geomfunc, option))
  }

  if (is.logical(shape) && !shape) {
    if (missing(label)) {
      label <- TRUE
    }
    if (missing(n.influential)) {
      n.influential <- nrow(plot.data)
    }
  }

  dev_ask <- is.null(nrow) & is.null(ncol)

  if(dev_ask){
    if (ask) {
      oask <- grDevices::devAskNewPage(TRUE)
      on.exit(grDevices::devAskNewPage(oask))
    }
  }

  plot.data <- flatten(plot.data)

  if (detect.influential) {
    bm <- local_influence_benchmarks(model, perturbation = pert[which], curvature = curvature,
                                     direction = direction, parameters = parameters,
                                     mean.covariates = mean.covariates,
                                     precision.covariates = precision.covariates)

    p.data <- list()

    for(i in 1:length(pert)){
      if(i %in% which){
        p.data[[i]] <- plot.data %>% dplyr::mutate(tmp = !!as.name(pert[[i]]), tmp = abs(tmp)) %>% dplyr::arrange(dplyr::desc(tmp))
        if(!is.na(bm[[pert[[i]]]])){
          p.data[[i]] <- p.data[[i]] %>% dplyr::filter(tmp > bm[[pert[[i]]]])
        }
        p.data[[i]] <- utils::head(p.data[[i]], n.influential) %>% dplyr::select(-tmp)
      }
    }
  }
  
  # Based on internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  .decorate.label.influential <- function(p, data) {
    if (label & n.influential > 0) {
      p <- plot_label_influential(p = p, data = data, label = label,
                                  label.label = ".label", label.colour = label.colour,
                                  label.alpha = label.alpha, label.size = label.size,
                                  label.angle = label.angle, label.family = label.family,
                                  label.fontface = label.fontface, label.lineheight = label.lineheight,
                                  label.hjust = label.hjust, label.vjust = label.vjust,
                                  label.repel = label.repel)
    }
    p
  }
  
  # Based on internal function from ggfortify package (https://github.com/sinhrks/ggfortify)

  .decorate.plot <- function(p, xlab = NULL, ylab = NULL, title = NULL) {
    p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
  }

  ylab_infl <- switch(curvature,
                      "conformal" = {
                        yl <- switch(direction,
                                     "canonical" = "Total Local Influence (Conformal)",
                                     "max.eigen" = "Largest Curvature Direction (Conformal)")
                        yl
                      },
                      "normal" = {
                        yl <- switch(direction,
                                     "canonical" = "Total Local Influence (Normal)",
                                     "max.eigen" = "Largest Curvature Direction (Normal)")
                      }
  )

  # Plots
  
  for(i in 1:length(pert)){
    if(i %in% which){

      if(include.modeltype){
        title[[i]] <- paste0(title[[i]], " - ", model$modeltype, " Regression")
      }

      t <- title[[i]]
      mapping <- ggplot2::aes_string(x = ".index", y = pert[[i]],
                                     ymin = 0, ymax = pert[[i]])
      p[[i]] <- ggplot2::ggplot(data = plot.data, mapping = mapping)
      if (!is.logical(shape) || shape) {
        p[[i]] <- p[[i]] + geom_factory_influential(geom_linerange, plot.data,
                                                    colour = colour, size = size, linetype = linetype,
                                                    alpha = alpha, fill = fill, shape = shape)
      }
      if(detect.influential){
        p[[i]] <- .decorate.label.influential(p[[i]], p.data[[i]])
      }

      if(sub.caption == "" | !dev_ask){
        xlab = "Obs. Number"
      } else {
        xlab = paste0("Obs. Number\n", sub.caption)
      }

      p[[i]] <- .decorate.plot(p[[i]], xlab = xlab, ylab = ylab_infl,
                               title = t)

      bm_i <- bm[[pert[[i]]]]

      if(draw.benchmark){
        if(!is.na(bm_i)){
          p[[i]] <- p[[i]] + ggplot2::geom_abline(intercept = bm_i, slope = 0,
                                                  linetype = ad.linetype, size = ad.size, colour = ad.colour)
        }
      }

      if(dev_ask){
        print(p[[i]])
      }

    }
  }


  if(!dev_ask){
    grDevices::devAskNewPage(ask = FALSE)
    plot.list <- lapply(which, function(i){p[[i]]})
    if(is.null(nrow))
      nrow <- 0
    if(is.null(ncol))
      ncol <- 0
    p <- methods::new("ggmultiplot", plots = plot.list, nrow = nrow, ncol = ncol)

    gpar_multi <- do.call(grid::gpar, gpar_sub.caption)
    title_multi <- grid::textGrob(sub.caption, gp=gpar_multi)
    gridExtra::grid.arrange(grobs = p@plots, top = title_multi)
    grDevices::devAskNewPage(ask = grDevices::dev.interactive())
  } else{
    invisible()
  }
}


#############################################################################################
#' @name tidy_local_influence
#' @title Tidy Functions for Local Influence Diagnostics
#' @aliases local_influence_autoplot tidy_local_influence local_influence_benchmarks
#' @description Functions to provide tidy outputs or ggplot2-based plots of local influence diagnostics.
#' @param model A model object for which local influence diagnostics are desired.
#' @param ... additional arguments to be passed.
#' @details
#' Local influence diagnostics were first introduced by Cook (1986), where several perturbation schemes were introduced and normal curvatures were obtained. Poon and Poon (1999)
#' introduced the conformal normal curvature, which has nice properties and takes values on the unit interval \eqn{[0,1]}. Zhu and Lee (2001) following Cook (1986) and Poon and Poon (1999)
#' introduced normal and conformal normal curvatures for EM-based models.
#' @references
#' Cook, R. D. (1986) *Assessment of Local Influence.* Journal of the Royal Statistical Society. Series B (Methodological), Vol. 48, pp.133-169. \href{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1986.tb01398.x}
#'
#' Poon, W.-Y. and Poon, Y.S. (1999) *Conformal normal curvature and assessment of local influence.*  Journal of the Royal Statistical Society. Series B (Methodological), Vol. 61, pp.51-61. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00162}
#'
#' Zhu, H.-T. and Lee, S.-Y. (2001) *Local influence for incomplete data models.* Journal of the Royal Statistical Society. Series B (Methodological), Vol. 63, pp.111-126. \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9868.00279}

#' @rdname tidy_local_influence
#' @export
local_influence_autoplot <- function(model, ...){
  UseMethod("local_influence_autoplot", model)
}

#' @rdname tidy_local_influence
#' @export
tidy_local_influence <- function(model, ...){
  UseMethod("tidy_local_influence", model)
}

#' @rdname tidy_local_influence
#' @export
local_influence_benchmarks <- function(model, ...){
  UseMethod("local_influence_benchmarks", model)
}

