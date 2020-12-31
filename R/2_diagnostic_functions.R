#' @importFrom lmtest coeftest coefci lrtest waldtest coefci.default coeftest.default

#############################################################################################
#' @name plot.mixpoissonreg
#' @title Plot Diagnostics for \code{mixpoissonreg} Objects
#' @description Currently there are six plots available. They contain residual analysis and global influence diagnostics. The plots are selectable by
#' the \code{which} argument. The plots are: Residuals vs. obs. numbers; Normal Q-Q plots, which may contain simulated envelopes, if the fitted object
#' has simulated envelopes; Cook's distances vs. obs. numbers; Generalized Cook's distances vs. obs. numbers; Cook's distances vs. Generalized Cook's distances;
#' Response variables vs. fitted means. By default, the first two plots and the last two plots are provided.
#' @param x object of class "mixpoissonreg" containing results from the fitted model.
#' If the model was fitted with envelope = 0, the Q-Q plot will be produced without envelopes.
#' @param which a list or vector indicating which plots should be displayed. 	If a subset of the plots is required, specify a subset of the numbers 1:6,
#' see caption below for the different kinds. In
#' plot number 2, 'Normal Q-Q', if the \code{mixpoissonreg} object was fitted with envelopes, a quantile-quantile plot with simulated envelopes will be displayed.
#' @param caption captions to appear above the plots; character vector or list of valid graphics annotations. Can be set to "" or NA to suppress all captions.
#' @param sub.caption common title-above the figures if there are more than one. If NULL, as by default, a possible abbreviated version of \code{deparse(x$call)} is used.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param id.n number of points to be labelled in each plot, starting with the most extreme.
#' @param main character; title to be placed at each plot additionally (and above) all captions.
#' @param labels.id	 vector of labels, from which the labels for extreme points will be chosen. The default uses the observation numbers.
#' @param label.pos positioning of labels, for the left half and right half of the graph respectively, for plots 2 and 6.
#' @param type.cookplot character; what type of plot should be drawn for Cook's and Generalized Cook's distances (plots 3 and 4). The default is 'h'.
#' @param cex.id magnification of point labels.
#' @param cex.caption	controls the size of caption.
#' @param cex.oma.main controls the size of the sub.caption only if that is above the figures when there is more than one.
#' @param include.modeltype logical. Indicates whether the model type ('NB' or 'PIG') should be displayed on the captions.
#' @param include.residualtype local. Indicates whether the name of the residual ('Pearson' or 'Score') should be displayed on the caption of plot 1 (Residuals vs. Index).
#' @param qqline logical; if \code{TRUE} and the fit does *not* contain simulated
#' envelopes, a qqline passing through the first and third quartiles of a standard normal distribution will be added to the normal Q-Q plot.
#' @param alpha_env alpha channel for envelope shading when the \code{mixpoissonreg} object was fitted with envelopes.
#' @param ... graphical parameters to be passed.
#' @details
#' The \code{plot} method is implemented following the same structure as the \link[stats]{plot.lm}, so it will be easy to be used by practitioners that
#' are familiar with \code{glm} objects.
#'
#' These plots allows one to perform residuals analsysis and influence diagnostics. There are other global influence functions, see \code{\link{influence.mixpoissonreg}}.
#'
#' See Barreto-Souza and Simas (2016), Cook and Weisberg (1982) and Zhu et al. (2001).
#'
#' @references
#' DOI:10.1007/s11222-015-9601-6 (\href{https://doi.org/10.1007/s11222-015-9601-6}{Barreto-Souza and Simas; 2016})
#'
#' Cook, D.R. and Weisberg, S. (1982) *Residuals and Influence in Regression*. (New York: Chapman and Hall, 1982)
#'
#' Zhu, H.T., Lee, S.Y., Wei, B.C., Zhu, J. (2001) *Case-deletion measures formodels with incomplete data.* Biometrika, 88, 727–737. \href{https://www.jstor.org/stable/2673442?seq=1}{https://www.jstor.org/stable/2673442?seq=1}
#'
#' @seealso
#' \code{\link{autoplot.mixpoissonreg}}, \code{\link{local_influence_plot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}},
#' \code{\link{summary.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}}, \code{\link{influence.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' plot(daysabs_fit, which = 1:6)
#'
#' par(mfrow = c(2,2))
#' plot(daysabs_fit)
#'
#' par(mfrow = c(1,1))
#' daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance, envelope = 100)
#' plot(daysabs_fit_ml, which = 2)
#' }
#' @export
plot.mixpoissonreg <- function(x, which = c(1,2,5,6),
                               caption = list("Residuals vs Obs. number",
                                              "Normal Q-Q",
                                              "Cook's distance",
                                              "Generalized Cook's distance",
                                              "Cook's dist vs Generalized Cook's dist",
                                              "Response vs Fitted means"
                                              ),
                               sub.caption = NULL, qqline = TRUE, alpha_env = 0.7,
                               main = "",
                               ask = prod(graphics::par("mfcol")) <
                                 length(which) && grDevices::dev.interactive(),
                               labels.id = names(stats::residuals(x)),
                               label.pos = c(4,2),
                               type.cookplot = 'h',
                               id.n = 3,
                               cex.id = 0.75,
                               cex.oma.main = 1.25,
                               cex.caption = 1,
                               include.modeltype = TRUE,
                               include.residualtype = FALSE,
                               ...) {

  getCaption <- function(k) if (length(caption) < k)
    NA_character_
  else {
    if(include.modeltype){
      grDevices::as.graphicsAnnot(paste0(caption[[k]], " - ", x$modeltype, " Regression"))
    } else {
      grDevices::as.graphicsAnnot(caption[[k]])
    }
  }

  if (is.null(sub.caption)) {
    cal <- x$call
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

  place_ids <- function(x_coord, y_coord, offset, dif_pos_neg){
    extreme_points <- as.vector(Rfast::nth(abs(y_coord), k = id.n,
                                    num.of.nths = id.n,
                                    index.return = TRUE, descending = TRUE))

    if(dif_pos_neg){
      idx_x_pos <- extreme_points[which(y_coord[extreme_points] >= 0)]
      idx_x_neg <- setdiff(extreme_points, idx_x_pos)
      idx_y_pos <- y_coord[idx_x_pos]
      idx_y_neg <- y_coord[idx_x_neg]
      idx_x_pos_id <- x_coord[idx_x_pos]
      idx_x_neg_id <- x_coord[idx_x_neg]
      if(length(idx_x_pos)>0){
        graphics::text(idx_x_pos_id, idx_y_pos, labels = labels.id[idx_x_pos], cex = cex.id, xpd = TRUE, pos = 3, offset = offset)
      }
      if(length(idx_x_neg)>0){
        graphics::text(idx_x_neg_id, idx_y_neg, labels = labels.id[idx_x_neg], cex = cex.id, xpd = TRUE, pos = 1, offset = offset)
      }
    } else{
      idx_x <- extreme_points
      idx_y <- y_coord[idx_x]
      idx_x_id <- x_coord[idx_x]
      labpos <- label.pos[1 + as.numeric(idx_x_id > mean(range(x_coord)))]
      graphics::text(idx_x_id, idx_y, labels = labels.id[idx_x], cex = cex.id, pos = labpos, xpd = TRUE, offset = offset)
    }
  }

  one.fig <- prod(graphics::par("mfcol")) == 1

  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }


  if (1 %in% which) {
    # First plot (residuals vs index)
    res <- stats::residuals(x, type = x$residualname)
    ylim <- range(res, na.rm = TRUE)
    if (id.n > 0)
      ylim <- grDevices::extendrange(r = ylim, f = 0.08)
    grDevices::dev.hold()
    residualname <- paste0(toupper(substring(x$residualname, 1, 1)), substring(x$residualname, 2))

    if(include.residualtype){
      caption[[1]] = paste(residualname, caption[[1]])
    }

    ylab <- paste0(residualname, " residuals")
    graphics::plot(res, ylab = ylab, xlab = "Obs. number", main = main, ylim = ylim, ...)
    graphics::abline(0, 0, lty = 3)

    if (one.fig)
      graphics::title(sub = sub.caption, ...)

    graphics::mtext(getCaption(1), side = 3, cex = cex.caption)

    place_ids(1:length(res), res, 0.5, TRUE)
    grDevices::dev.flush()
  }

  # Second plot (QQ-plot)
  if (2 %in% which) {
    env <- x$envelope
    res <- stats::residuals(x, type = x$residualname)
    ylim <- range(res, env, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    grDevices::dev.hold()
    n <- length(x$residuals)
    residualname <- paste0(toupper(substring(x$residualname, 1, 1)), substring(x$residualname, 2))
    ylab <- paste0(residualname, " residuals")
    if (!is.null(env)) {
      caption[[2]] <- paste0(caption[[2]]," with simulated envelopes")
    }
    qq <- stats::qqnorm(res, xlab = "Theoretical quantiles", ylab = ylab, ylim = ylim, main = main, ...)

    if (one.fig)
      graphics::title(sub = sub.caption, ...)

    if (is.null(env) == FALSE) {
      aux <- sort(qq$x)
      graphics::lines(aux, env[1, ], col = grDevices::rgb(0.7, 0.7, 0.7))
      graphics::lines(aux, env[3, ], col = grDevices::rgb(0.7, 0.7, 0.7))
      graphics::polygon(c(aux, rev(aux)), c(env[3, ], rev(env[1, ])), col = grDevices::rgb(0.7, 0.7, 0.7, alpha_env), border = NA)
      graphics::lines(aux, env[2, ], lty = 2, lwd = 2)
    } else {
      if (qqline) {
        stats::qqline(res, lty = 3)
      }
    }
    graphics::points(qq$x, qq$y, ...)
    if (id.n > 0)
      place_ids(qq$x, qq$y, 0.5, FALSE)
    graphics::mtext(getCaption(2), side = 3, cex = cex.caption)
    grDevices::dev.flush()
  }

  # Third plot (Cook's distance)

  if (3 %in% which) {
    CD <- stats::cooks.distance(x, type = "CD")

    ylab = "Cook's distance"
    ylim <- range(CD, na.rm = TRUE)

    grDevices::dev.hold()

    if(id.n >0)
      ylim <- grDevices::extendrange(r = ylim, f = 0.08)

    graphics::plot(CD, type = type.cookplot, main = main, xlab = "Obs. number", ylab = ylab, ylim=ylim, ...)
    if (one.fig)
      graphics::title(sub = sub.caption, ...)

    graphics::mtext(getCaption(3), side = 3, cex = cex.caption)

    if (id.n > 0)
      place_ids(1:length(CD), CD, 0.2, TRUE)
    grDevices::dev.flush()
  }

  # Fourth plot (Generalized Cook's distance)

  if (4 %in% which) {
    GCD <- stats::cooks.distance(x, type = "GCD")

    ylab = "Generalized Cook's distance"
    ylim <- range(GCD, na.rm = TRUE)

    grDevices::dev.hold()

    if(id.n >0)
      ylim <- grDevices::extendrange(r = ylim, f = 0.08)

    graphics::plot(GCD, type = type.cookplot, main = main, ylab = ylab,
                   xlab = "Obs. number", ylim=ylim, ...)
    if (one.fig)
      graphics::title(sub = sub.caption, ...)

    graphics::mtext(getCaption(4), side = 3, cex = cex.caption)

    if (id.n > 0)
      place_ids(1:length(GCD), GCD, 0.2, TRUE)
    grDevices::dev.flush()
  }

  # Fifth plot (Cook's dist vs Generalized Cook's dist)

  if(5 %in% which) {
    CD <- stats::cooks.distance(x, type = "CD")
    GCD <- stats::cooks.distance(x, type = "GCD")

    ylim <- range(CD, na.rm = TRUE)

    if (id.n > 0)
      ylim <- grDevices::extendrange(r = ylim, f = 0.08)

    grDevices::dev.hold()

    graphics::plot(GCD, CD, main = main, ylab = "Cook's distance", xlab = "Generalized Cook's distance", ylim=ylim, ...)
    if (one.fig)
      graphics::title(sub = sub.caption, ...)

    graphics::mtext(getCaption(5), side = 3, cex = cex.caption)

    if (id.n > 0){
      extreme_points <- as.vector(Rfast::nth((CD/sum(CD))^2 + (GCD/sum(GCD))^2, k = id.n,
                                      num.of.nths = id.n,
                                      index.return = TRUE, descending = TRUE))
      labpos <- label.pos[1 + as.numeric(GCD[extreme_points] > mean(range(GCD)))]
      graphics::text(GCD[extreme_points], CD[extreme_points], labels = labels.id[extreme_points], cex = cex.id, pos = labpos,
                     xpd = TRUE, offset = 0.5)
    }

    grDevices::dev.flush()
  }

  # Sixth plot (Response vs Fitted means)

  if (6 %in% which) {
    mu_est <- stats::fitted(x, type = "response")
    if(is.null(x$y)){
      y = x$residuals + x$fitted.values
    } else{
      y = x$y
    }

    graphics::plot(mu_est, y, xlab = "Predicted values", ylab = "Response values", main = main, ...)
    if (one.fig)
      graphics::title(sub = sub.caption, ...)
    graphics::abline(0, 1, lty = 3)
    graphics::mtext(getCaption(6), side = 3, cex = cex.caption)

    if (id.n > 0){
      extreme_points <- as.vector(Rfast::nth(abs(y - mu_est), k = id.n,
                                      num.of.nths = id.n,
                                      index.return = TRUE, descending = TRUE))
      labpos <- label.pos[1 + as.numeric(mu_est[extreme_points] > mean(range(mu_est)))]
      graphics::text(mu_est[extreme_points], y[extreme_points], labels = labels.id[extreme_points], cex = cex.id, pos = labpos,
                     xpd = TRUE, offset = 0.5)
    }

    grDevices::dev.flush()
  }

  if (!one.fig && graphics::par("oma")[3L] >= 1)
    graphics::mtext(sub.caption, outer = TRUE, cex = 1.25)

  invisible()
}


#############################################################################################
#' @name fitted.mixpoissonreg
#' @title Fitted Method for \code{mixpoissonreg} Objects
#' @description Function providing the fitted means, linear predictors, precisions or variances for mixed Poisson regression models.
#' @param object object of class "mixpoissonreg" containing results from the fitted model.
#' @param type the type of variable to get the fitted values. The default is the "response" type, which provided the estimated values for the means.
#' The type "link" provides the estimates for the linear predictor of the mean. The type "precision" provides estimates for the precision parameters
#' whereas the type "variance" provides estimates for the variances.
#' @param ... Currently not used.
#' @seealso
#' \code{\link{predict.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}},
#' \code{\link{coef.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' fitted(daysabs_fit)
#' fitted(daysabs_fit, type = "precision")
#' }
#' @export
fitted.mixpoissonreg <- function(object, type = c("response", "link", "precision", "variance"), ...) {
  fit <- object
  if (length(type) > 1) {
    type <- type[1]
  }
  possible_types <- c("response", "link", "precision", "variance")
  if (!(type %in% c("response", "link", "precision", "variance"))) {
    stop(paste0("type must be one of ", possible_types))
  }

  model <- fit$modeltype

  ddB <- switch(model,
         "PIG" = {function(x) {1/( (-2*x)^(3/2))}},
         "NB" = {function(x) {1/(x^2)}})

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  link_mean <- build_links_mpreg(fit$link.mean)

  link_precision <- build_links_mpreg(fit$link.precision)

  fitted_values <- switch(type,
                          "response" = {
                            mu <- c(fit$fitted.values)
                            names(mu) <- 1:length(mu)
                            mu
                          },
                          "link" = {
                            link_fitted <- c(link_mean$linkfun(fit$fitted.values))
                            names(link_fitted) <- 1:length(link_fitted)
                            link_fitted
                          },
                          "precision" = {
                            fitted_prec <- c(fit$fitted.precisions)
                            names(fitted_prec) <- 1:length(fitted_prec)
                            fitted_prec
                          },
                          "variance" = {
                            fitted_prec <- c(fit$fitted.precisions)
                            variance_fitted <- c(fit$fitted.values * (1 + ddB(qsi0) * fit$fitted.values / fitted_prec))
                            names(variance_fitted) <- 1:length(variance_fitted)
                            variance_fitted
                          }
  )
  return(fitted_values)
}


#############################################################################################
#' @name predict.mixpoissonreg
#' @title Predict Method for \code{mixpoissonreg} Objects
#' @description Function to obtain various predictions based on the fitted mixed Poisson regression models.
#' @param object object of class "mixpoissonreg" containing results from the fitted model.
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted response values will be provided.
#' @param type the type of prediction. The default is the "response" type, which provided the estimated values for the means. The type "link" provides the estimates for the linear predictor. The type "precision" provides estimates for the precision parameters whereas the type "variance" provides estimates for the variances.
#' @param se.fit logical switch indicating if standard errors on the scale of linear predictors should be returned. If \code{TRUE}, it only returns the standard deviations of
#' the linear predictors when type = 'link', otherwise returns NA and a warning indicating that the type must be 'link'.  When using \code{mixpoissonreg} objects, the fit must be done with \code{x = TRUE}.
#' @param interval Type of interval calculation for the response variables, 'none', 'confidence' or 'prediction'. If 'confidence', the confidence intervals for the means are returned.
#' If 'prediction', prediction intervals for future response variables are reported. For confidence intervals, the type of the prediction must be 'response' or 'link'.
#' For prediction intervals the type of prediction must be 'response'. For 'confidence' intervals, when using \code{mixpoissonreg} objects, the fit must be done with \code{x = TRUE} and for
#' predictions intervals, the fit must be done with \code{x = TRUE} and \code{w = TRUE}.
#' @param level Tolerance/confidence level. The default is set to 0.95.
#' @param nsim_pred number of means and predictions to be generated in each step of the simulation. The default is set to 100.
#' @param nsim_pred_y number of response variables generated for each pair of mean and precision to compute the prediction intervals. The default is set to 100.
#' @param ... further arguments passed to or from other methods.
#' @details
#' The \code{se.fit} argument only returns a non-NA vector for type = 'link', that is, on the scale of the linear predictor for the mean parameter. For the response scale,
#' one can obtain confidence or prediction intervals. It is important to notice that confidence intervals *must not* be used for future observations as they will underestimate
#' the uncertainty. In this case prediction intervals should be used. Currently, we do not have closed-form expressions for the prediction interval and, therefore, they
#' are obtained by simulation and can be computationally-intensive.
#'
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}}, \code{\link{plot.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}},
#' \code{\link{coef.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' predict(daysabs_fit, interval = "confidence")
#' predict(daysabs_fit, type = "link", se.fit = TRUE)
#' }
#' @export
predict.mixpoissonreg <- function(object, newdata = NULL, type = c("response", "link", "precision", "variance"), se.fit = FALSE,
                                  interval = c("none", "confidence", "prediction"), level = 0.95, nsim_pred = 100, nsim_pred_y = 100, ...) {
  fit <- object
  if (length(type) > 1) {
    type <- type[1]
  }

  if(level > 1 | level < 0){
    stop("level must be a number between 0 and 1.")
  }

  interval <- rlang::arg_match(interval)

  possible_intervals <- c("none", "confidence", "prediction")

  if(!(interval %in% possible_intervals)){
    stop(paste("interval should be one of", possible_intervals))
  }

  possible_types <- c("response", "link", "precision", "variance")
  if (!(type %in% c("response", "link", "precision", "variance"))) {
    stop(paste0("type must be one of ", possible_types))
  }


  if (missing(newdata)) {
    predictions <- stats::fitted(fit, type)
    if(interval != "none"){
      if(is.null(fit$x)){
        stop('the fit must contain the matrix x, run the fit again with argument x = TRUE')
      }
      pred_matrix <- switch(interval,
                            "confidence" = {
                              x_matrix <- object$x
                              std_error <- sqrt(diag(x_matrix%*%stats::vcov(fit, parameters = "mean")%*%t(x_matrix)))
                              pred_matrix <- switch(type,
                                                    "response" = {
                                                      link_mean <- build_links_mpreg(fit$link.mean)
                                                      lwr_bound <- link_mean$linkinv(stats::fitted(fit, type = 'link') + stats::qnorm( (1-level)/2 )*std_error)
                                                      upr_bound <- link_mean$linkinv(stats::fitted(fit, type = 'link') + stats::qnorm( (1+level)/2 )*std_error)
                                                      cbind(stats::fitted(fit, type = "response"), lwr_bound, upr_bound)
                                                    },
                                                    "link" = {
                                                      lwr_bound <- stats::fitted(fit, type = 'link') + stats::qnorm( (1-level)/2 )*std_error
                                                      upr_bound <- stats::fitted(fit, type = 'link') + stats::qnorm( (1+level)/2 )*std_error
                                                      cbind(stats::fitted(fit, type = "link"), lwr_bound, upr_bound)
                                                    },
                                                    "variance" = {
                                                      warning("for confidence intervals the type must be 'response' or 'link'")
                                                      lwr_bound <-rep(NA, length(fit$fitted.values))
                                                      upr_bound <- rep(NA, length(fit$fitted.values))
                                                      cbind(stats::fitted(fit, type = "variance"), lwr_bound, upr_bound)
                                                    },
                                                    "precision" = {
                                                      warning("for confidence intervals the type must be 'response' or 'link'")
                                                      lwr_bound <-rep(NA, length(fit$fitted.values))
                                                      upr_bound <- rep(NA, length(fit$fitted.values))
                                                      cbind(stats::fitted(fit, type = "precision"), lwr_bound, upr_bound)
                                                    })
                              colnames(pred_matrix) <- c("fit", "lwr", "upr")
                              rownames(pred_matrix) <- 1:length(fit$fitted.values)
                              pred_matrix
                            },
                            "prediction" = {
                              if(is.null(fit$w)){
                                stop('the fit must contain the matrix w, run the fit again with argument w = TRUE')
                              }
                              pred_matrix <- switch(type,
                                                    "response" = {
                                                      warning("predictions on current data refer to _future_ responses")
                                                      x_matrix <- object$x
                                                      w_matrix <- object$w
                                                      alpha_est <- fit$coefficients$precision

                                                      link_mean <- build_links_mpreg(fit$link.mean)
                                                      link_precision <- build_links_mpreg(fit$link.precision)

                                                      mean_link_scale <- stats::fitted(fit, type = 'link')
                                                      precision_link_scale <- fit$w %*% alpha_est
                                                      std_error_mean <- sqrt(diag(x_matrix%*%stats::vcov(object, parameters = "mean")%*%t(x_matrix)))
                                                      std_error_precision <- sqrt(diag(w_matrix%*%stats::vcov(object, parameters = "precision")%*%t(w_matrix)))

                                                      pred_matrix <- switch(fit$modeltype,
                                                                            "NB" = {
                                                                              pred_matrix <- lapply(1:length(fit$fitted.values), function(i){
                                                                                mu_pred <- link_mean$linkinv(stats::rnorm(nsim_pred, mean = mean_link_scale[i], sd = std_error_mean[i]))
                                                                                phi_pred <- link_precision$linkinv(stats::rnorm(nsim_pred, mean = precision_link_scale[i], sd = std_error_precision[i]))
                                                                                phi_pred[phi_pred <= 0] = 10^(-5)
                                                                                y_sim <- sapply(1:nsim_pred, function(j){rNBI(nsim_pred_y, mu = mu_pred[j], sigma = 1/phi_pred[j])})
                                                                                y_sim <- sort(y_sim)
                                                                                idx_lwr <- max(1, round(nsim_pred * nsim_pred_y * (1 - level) / 2))
                                                                                idx_upr <- round(nsim_pred * nsim_pred_y * (1 + level) / 2)
                                                                                c(y_sim[idx_lwr], y_sim[idx_upr])
                                                                              })
                                                                            },
                                                                            "PIG" = {
                                                                              pred_matrix <- lapply(1:length(fit$fitted.values), function(i){
                                                                                mu_pred <- link_mean$linkinv(stats::rnorm(nsim_pred, mean = mean_link_scale[i], sd = std_error_mean[i]))
                                                                                phi_pred <- link_precision$linkinv(stats::rnorm(nsim_pred, mean = precision_link_scale[i], sd = std_error_precision[i]))
                                                                                phi_pred[phi_pred <= 0] = 10^(-5)
                                                                                y_sim  = sapply(1:nsim_pred, function(j){ ig <- rIG(nsim_pred_y,mu=1,sigma=1/sqrt(phi_pred[j]))
                                                                                stats::rpois(nsim_pred_y,ig*mu_pred[j])
                                                                                })
                                                                                y_sim <- sort(y_sim)
                                                                                idx_lwr <- max(1, round(nsim_pred * nsim_pred_y * (1 - level) / 2))
                                                                                idx_upr <- round(nsim_pred * nsim_pred_y * (1 + level) / 2)
                                                                                c(y_sim[idx_lwr], y_sim[idx_upr])
                                                                              })
                                                                            })
                                                      pred_matrix <- unlist(pred_matrix)
                                                      pred_matrix <- t(matrix(pred_matrix, nrow = 2))
                                                      cbind(stats::fitted(fit, type = "response"), pred_matrix)
                                                    },
                                                    "link" = {
                                                      warning("for prediction intervals the type must be 'response'")
                                                      lwr_bound <-rep(NA, length(fit$fitted.values))
                                                      upr_bound <- rep(NA, length(fit$fitted.values))
                                                      cbind(stats::fitted(fit, type = "link"), lwr_bound, upr_bound)
                                                    },
                                                    "precision" = {
                                                      warning("for prediction intervals the type must be 'response'")
                                                      lwr_bound <-rep(NA, length(fit$fitted.values))
                                                      upr_bound <- rep(NA, length(fit$fitted.values))
                                                      cbind(stats::fitted(fit, type = "precision"), lwr_bound, upr_bound)
                                                    },
                                                    "variance" = {
                                                      warning("for prediction intervals the type must be 'response'")
                                                      lwr_bound <-rep(NA, length(fit$fitted.values))
                                                      upr_bound <- rep(NA, length(fit$fitted.values))
                                                      cbind(stats::fitted(fit, type = "variance"), lwr_bound, upr_bound)
                                                    })
                              colnames(pred_matrix) <- c("fit", "lwr", "upr")
                              rownames(pred_matrix) <- 1:length(fit$fitted.values)
                              pred_matrix
                            })
      predictions <- pred_matrix
    }

    if(se.fit){
      if(is.null(fit$x)){
        stop('the fit must contain the matrix x, run the fit again with argument x = TRUE')
      }
      predictions_temp <- list()
      predictions_temp$fit <- predictions
      if(type != "link"){
        warning("se.fit should be TRUE only with type = 'link'.")
        std_error <- rep(NA, length(object$y))
        names(std_error) = 1:length(object$y)
      } else{
        x_matrix <- object$x
        std_error <- sqrt(diag(x_matrix%*%stats::vcov(object, parameters = "mean")%*%t(x_matrix)))
      }
      predictions_temp$se.fit <- std_error
      predictions <- predictions_temp
    }
  } else {
    formula_temp <- Formula(fit$formula)
    matrix_temp_x <- stats::model.matrix(object = formula_temp, data = newdata, rhs = 1)
    matrix_temp_w <- stats::model.matrix(object = formula_temp, data = newdata, rhs = 2)

    beta_est <- fit$coefficients$mean
    alpha_est <- fit$coefficients$precision

    link_mean <- build_links_mpreg(fit$link.mean)
    link_precision <- build_links_mpreg(fit$link.precision)

    mu_est <- link_mean$linkinv(matrix_temp_x %*% beta_est)
    mu_est <- c(mu_est)
    names(mu_est) <- 1:length(mu_est)

    phi_est <- c(link_precision$linkinv(matrix_temp_w %*% alpha_est))
    names(phi_est) <- 1:length(phi_est)

    model <- fit$modeltype

    ddB <- switch(model,
                  "PIG" = {function(x) {1/( (-2*x)^(3/2))}},
                  "NB" = {function(x) {1/(x^2)}})

    qsi0 <- switch(model,
                   "NB" = {-1},
                   "PIG" = {-1/2}
    )

    predictions <- switch(type,
                          "response" = {
                            mu_est
                          },
                          "link" = {
                            link_predict <- c(matrix_temp_x %*% beta_est)
                            names(link_predict) <- 1:length(link_predict)
                            link_predict
                          },
                          "precision" = {
                            phi_est
                          },
                          "variance" = {
                            variance_fitted <- c(mu_est * (1 + ddB(qsi0) * mu_est / phi_est))
                            names(variance_fitted) <- 1:length(variance_fitted)
                            variance_fitted
                          }
    )

    if(interval != "none"){
      if(is.null(fit$x)){
        stop('the fit must contain the matrix x, run the fit again with argument x = TRUE')
      }
      pred_matrix <- switch(interval,
                            "confidence" = {
                              std_error <- sqrt(diag(matrix_temp_x%*%stats::vcov(object, parameters = "mean")%*%t(matrix_temp_x)))
                              pred_matrix <- switch(type,
                                                    "response" = {
                                                      link_mean <- build_links_mpreg(fit$link.mean)
                                                      lwr_bound <- link_mean$linkinv(matrix_temp_x %*% beta_est + stats::qnorm( (1-level)/2 )*std_error)
                                                      upr_bound <- link_mean$linkinv(matrix_temp_x %*% beta_est + stats::qnorm( (1+level)/2 )*std_error)
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    },
                                                    "link" = {
                                                      lwr_bound <- predictions + stats::qnorm( (1-level)/2 )*std_error
                                                      upr_bound <- predictions + stats::qnorm( (1+level)/2 )*std_error
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    },
                                                    "variance" = {
                                                      warning("for confidence intervals the type must be 'response' or 'link'")
                                                      lwr_bound <-rep(NA, nrow(newdata))
                                                      upr_bound <- rep(NA, nrow(newdata))
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    },
                                                    "precision" = {
                                                      warning("for confidence intervals the type must be 'response' or 'link'")
                                                      lwr_bound <-rep(NA, nrow(newdata))
                                                      upr_bound <- rep(NA, nrow(newdata))
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    })
                              colnames(pred_matrix) <- c("fit", "lwr", "upr")
                              rownames(pred_matrix) <- 1:nrow(newdata)
                              pred_matrix
                            },
                            "prediction" = {
                              if(is.null(fit$w)){
                                stop('the fit must contain the matrix w, run the fit again with argument w = TRUE')
                              }
                              pred_matrix <- switch(type,
                                                    "response" = {
                                                      warning("predictions on current data refer to _future_ responses")
                                                      beta_est <- fit$coefficients$mean
                                                      alpha_est <- fit$coefficients$precision

                                                      link_mean <- build_links_mpreg(fit$link.mean)
                                                      link_precision <- build_links_mpreg(fit$link.precision)

                                                      mean_link_scale <- matrix_temp_x %*% beta_est
                                                      precision_link_scale <- matrix_temp_w %*% alpha_est
                                                      std_error_mean <- sqrt(diag(matrix_temp_x%*%stats::vcov(object, parameters = "mean")%*%t(matrix_temp_x)))
                                                      std_error_precision <- sqrt(diag(matrix_temp_w%*%stats::vcov(object, parameters = "precision")%*%t(matrix_temp_w)))

                                                      pred_matrix <- switch(fit$modeltype,
                                                                            "NB" = {
                                                                              pred_matrix <- lapply(1:nrow(newdata), function(i){
                                                                                mu_pred <- link_mean$linkinv(stats::rnorm(nsim_pred, mean = mean_link_scale[i], sd = std_error_mean[i]))
                                                                                phi_pred <- link_precision$linkinv(stats::rnorm(nsim_pred, mean = precision_link_scale[i], sd = std_error_precision[i]))
                                                                                phi_pred[phi_pred <= 0] = 10^(-5)
                                                                                y_sim <- sapply(1:nsim_pred, function(j){rNBI(nsim_pred_y, mu = mu_pred[j], sigma = 1/phi_pred[j])})
                                                                                y_sim <- sort(y_sim)
                                                                                idx_lwr <- max(1, round(nsim_pred * nsim_pred_y * (1 - level) / 2))
                                                                                idx_upr <- round(nsim_pred * nsim_pred_y * (1 + level) / 2)
                                                                                c(y_sim[idx_lwr], y_sim[idx_upr])
                                                                              })
                                                                            },
                                                                            "PIG" = {
                                                                              pred_matrix <- lapply(1:nrow(newdata), function(i){
                                                                                mu_pred <- link_mean$linkinv(stats::rnorm(nsim_pred, mean = mean_link_scale[i], sd = std_error_mean[i]))
                                                                                phi_pred <- link_precision$linkinv(stats::rnorm(nsim_pred, mean = precision_link_scale[i], sd = std_error_precision[i]))
                                                                                phi_pred[phi_pred <= 0] = 10^(-5)
                                                                                y_sim  = sapply(1:nsim_pred, function(j){ ig <- rIG(nsim_pred_y,mu=1,sigma=1/sqrt(phi_pred[j]))
                                                                                stats::rpois(nsim_pred_y,ig*mu_pred[j])
                                                                                })
                                                                                y_sim <- sort(y_sim)
                                                                                idx_lwr <- max(1, round(nsim_pred * nsim_pred_y * (1 - level) / 2))
                                                                                idx_upr <- round(nsim_pred * nsim_pred_y * (1 + level) / 2)
                                                                                c(y_sim[idx_lwr], y_sim[idx_upr])
                                                                              })
                                                                            })
                                                      pred_matrix <- unlist(pred_matrix)
                                                      pred_matrix <- t(matrix(pred_matrix, nrow = 2))
                                                      cbind(predictions, pred_matrix)
                                                    },
                                                    "link" = {
                                                      warning("for prediction intervals the type must be 'response'")
                                                      lwr_bound <-rep(NA, nrow(newdata))
                                                      upr_bound <- rep(NA, nrow(newdata))
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    },
                                                    "precision" = {
                                                      warning("for prediction intervals the type must be 'response'")
                                                      lwr_bound <-rep(NA, nrow(newdata))
                                                      upr_bound <- rep(NA, nrow(newdata))
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    },
                                                    "variance" = {
                                                      warning("for prediction intervals the type must be 'response'")
                                                      lwr_bound <-rep(NA, nrow(newdata))
                                                      upr_bound <- rep(NA, nrow(newdata))
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    })
                              colnames(pred_matrix) <- c("fit", "lwr", "upr")
                              rownames(pred_matrix) <- 1:nrow(newdata)
                              pred_matrix
                            })
      predictions <- pred_matrix
    }

    if(se.fit){
      if(is.null(fit$x)){
        stop('the fit must contain the matrix x, run the fit again with argument x = TRUE')
      }
      predictions_temp <- list()
      predictions_temp$fit <- predictions
      if(type != "link"){
        warning("se.fit should be TRUE only with type = 'link'.")
        std_error <- rep(NA, nrow(newdata))
        names(std_error) = 1:nrow(newdata)
      } else{
        std_error <- sqrt(diag(matrix_temp_x%*%stats::vcov(object, parameters = "mean")%*%t(matrix_temp_x)))
      }
      predictions_temp$se.fit <- std_error
      predictions <- predictions_temp
    }

  }
  return(predictions)
}

#############################################################################################
#' @name vcov.mixpoissonreg
#' @title Calculate Variance-Covariance Matrix for \code{mixpoissonreg} Objects
#' @description Returns the variance-covariance matrix of the parameters for fitted mixed Poisson regression models. The \code{parameters} argument
#' indicates for which parameters the variance-covariance matrix should be computed, namely, 'mean' for mean-relatex parameters or 'precision' for precision-related parameters.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param parameters a string to determine which coefficients should be extracted: 'all' extracts all coefficients, 'mean' extracts the coefficients of the mean parameters and 'precision' extracts coefficients of the precision parameters.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{coef.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' vcov(daysabs_fit)
#' vcov(daysabs_fit, parameters = "mean")
#' }
#' @export
vcov.mixpoissonreg <- function(object, parameters = c("all", "mean", "precision"), ...) {
  if (length(parameters) > 1) {
    parameters <- parameters[1]
  }
  nbeta <- length(object$coefficients$mean)
  nalpha <- length(object$coefficients$precision)
  theta <- c(object$coefficients$mean, object$coefficients$precision)

  coeff_beta <- object$coefficients$mean
  coeff_alpha <- object$coefficients$precision
  coeff_names = c(names(coeff_beta),paste(names(coeff_alpha), ".precision", sep = ""))

  vcov_mp_complete <- object$vcov
  vcov_mp <- switch(parameters,
                    "all" = {
                      colnames(vcov_mp_complete) <- rownames(vcov_mp_complete) <- coeff_names
                      vcov_mp_complete
                    },
                    "mean" = {
                      vcov_mp_complete[1:nbeta, 1:nbeta]
                    },
                    "precision" = {
                      vcov_mp_complete[(nbeta + 1):(nbeta + nalpha), (nbeta + 1):(nbeta + nalpha)]
                    }
  )

  return(vcov_mp)
}

#############################################################################################
#' @name coef.mixpoissonreg
#' @title Coef Method for \code{mixpoissonreg} Objects.
#' @description Extract model coefficients of fitted mixed Poisson regression models. The parameters arguments allows one to chose if all coefficients should be extracted,
#' with \code{parameters = 'all'}; if the coefficients of the mean-related parameters should be extracted, with \code{parameters = 'mean'}; if the coefficients of the
#' precision-related parameters should be extracted, with \code{parameters = 'precision'}.
#' @param object object of class "mixpoissonreg" containing results from the fitted model.
#' @param parameters a string to determine which coefficients should be extracted: 'all' extracts all coefficients, 'mean' extracts the coefficients of the mean parameters and 'precision' extracts coefficients of the precision parameters.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{vcov.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' coef(daysabs_fit)
#' coef(daysabs_fit, parameters = "precision")
#' }
#' @export
coef.mixpoissonreg <- function(object, parameters = c("all", "mean", "precision"), ...) {
  parameters <- rlang::arg_match(parameters)

  coef_ext <- switch(parameters,
                     "all" = {
                       coeff_beta <- object$coefficients$mean
                       coeff_alpha <- object$coefficients$precision
                       coeff_names = c(names(coeff_beta),paste(names(coeff_alpha), ".precision", sep = ""))
                       coef_complete <- c(coeff_beta, coeff_alpha)
                       names(coef_complete) <- coeff_names
                       coef_complete
                     }, "mean" = {
                       object$coefficients$mean
                     }, "precision" = {
                       object$coefficients$precision
                     }
  )
  return(coef_ext)
}



#############################################################################################
#' @name terms.mixpoissonreg
#' @title Terms Method for \code{mixpoissonreg} Objects.
#' @description Function to extract the terms of a fitted mixpoissonreg object.
#' @param x an object of class "mixpoissonreg" containing results from the fitted model.
#' @param parameters characters the parameters to be chosen. The options are 'mean' and 'precision'.
#' @param ... Currently not used.
#' @noRd
#' @export

terms.mixpoissonreg <- function(x, parameters = c("mean", "precision"), ...){
  parameters <- rlang::arg_match(parameters)
  x$terms[[parameters]]
}


#############################################################################################
#' @name summary.mixpoissonreg
#' @title Summary Method for \code{mixpoissonreg} Objects.
#' @description Function providing a summary of results related to mixed Poisson regression models.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{plot.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}},
#' \code{\link{local_influence_plot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' summary(daysabs_fit)
#'
#' daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' summary(daysabs_fit_ml)
#' }
#' @export
summary.mixpoissonreg <- function(object, ...) {
  ans <- list()

  nbeta <- length(object$coefficients$mean)
  nalpha <- length(object$coefficients$precision)
  #
  coeff <- object$coefficients
  coeff <- c(coeff$mean, coeff$precision)
  SEr <- object$std_errors
  tab <- cbind(coeff, SEr, coeff / SEr, 2 * stats::pnorm(-abs(coeff / SEr)))
  colnames(tab) <- c("Estimate", "Std.error", "z-value", "Pr(>|z|)")
  rownames(tab) <- names(coeff)
  tab <- list(mean = tab[seq.int(length.out = nbeta), , drop = FALSE], precision = tab[seq.int(length.out = nalpha) + nbeta, , drop = FALSE])

  ans$coefficients <- tab

  ans$estimation_method <- object$estimation_method

  ans$modelname <- object$modelname

  ans$call <- object$call

  ans$residualname <- paste0(toupper(substring(object$residualname, 1, 1)), substring(object$residualname, 2))

  ans$RSS <- sum(stats::residuals(object, type = object$residualname)^2)

  ans$res_quantiles <- stats::quantile(stats::residuals(object, type = object$residualname))

  ans$efron.pseudo.r2 <- object$efron.pseudo.r2

  ans$envelope_prop <- object$envelope_prop

  ans$envelope <- object$envelope

  ans$niter <- object$niter

  class(ans) <- "summary_mixpoissonreg"
  ans
}


#############################################################################################
#' @name print.mixpoissonreg
#' @title Print Method for \code{mixpoissonreg} Objects
#' @description Provides a brief description of results related to mixed Poisson regression models.
#' @param x object of class "mixpoissonreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @noRd
#' @export

print.mixpoissonreg <- function(x, ...) {
  nbeta <- length(x$coefficients$mean)
  nalpha <- length(x$coefficients$precision)
  #
  call_name <- switch(x$estimation_method,
                      "EM" = {"mixpoissonreg"},
                      "ML" = {"mixpoissonregML"}
  )

  coeff_beta <- x$coefficients$mean
  coeff_alpha <- x$coefficients$precision
  optim_algo <- switch(x$estimation_method,
                       "EM" = {"Expectation-Maximization Algorithm"},
                       "ML" = {"Maximum-Likelihood Estimation"}
  )
  cat("\n")
  cat(paste0(x$modelname, " - ", optim_algo))
  cat("\n\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(paste0("Coefficients modeling the mean (with ", x$link.mean, " link):", "\n"))
  print(coeff_beta)
  cat(paste0("Coefficients modeling the precision (with ", x$link.precision, " link):", "\n"))
  print(coeff_alpha)
}



#############################################################################################
#' @name print.summary_mixpoissonreg
#' @title Print Method for \code{summary_mixpoissonreg} Objects
#' @description Provides a brief description of results related to mixed Poisson regression models.
#' @param x object of class "summary_mixpoissonreg" containing results of summary method applied to a fitted model.
#' @param ... further arguments passed to or from other methods.
#' @noRd
#' @export
print.summary_mixpoissonreg <- function(x, ...) {
  tab <- x$coefficients

  #
  digits <- max(3, getOption("digits") - 3)
  #

  call_name <- switch(x$estimation_method,
                      "EM" = {"mixpoissonreg"},
                      "ML" = {"mixpoissonregML"}
  )

  optim_algo <- switch(x$estimation_method,
                       "EM" = {"Expectation-Maximization Algorithm"},
                       "ML" = {"Maximum-Likelihood Estimation"}
  )
  cat("\n")
  cat(paste0(x$modelname, " - ", optim_algo))

  cat("\n\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")


  #
  RSS <- x$RSS
  residualname <- x$residualname

  cat(sprintf("\n%s:\n", paste0(residualname, " residuals")))
  print(structure(round(c(RSS, as.vector(x$res_quantiles)), digits = digits), .Names = c("RSS", "Min", "1Q", "Median", "3Q", "Max")))
  #
  if (NROW(tab$mean)) {
    cat(paste0("\nCoefficients modeling the mean (with ", x$link.mean, " link):\n"))
    stats::printCoefmat(tab$mean, digits = digits, signif.legend = FALSE)
  } else {
    message("\nNo coefficients modeling the mean. \n")
  }
  #
  if (NROW(tab$precision)) {
    cat(paste0("\nCoefficients modeling the precision (with ", x$link.precision, " link):\n"))
    stats::printCoefmat(tab$precision, digits = digits, signif.legend = FALSE)
  } else {
    message("\nNo coefficients modeling the precision. \n")
  }
  #
  if (getOption("show.signif.stars")) {
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n\n")
  }
  #

  cat("Efron's pseudo R-squared: ", x$efron.pseudo.r2,"\n")

  if (is.null(x$envelope) == FALSE) {
    cat(sprintf("%s\n", paste0("Percentage of residuals within the envelope = ", round(x$envelope_prop, digits = digits))))
  }
  if(x$estimation_method == "EM"){
    cat(paste0("Number of iterations of the EM algorithm = ", x$niter,"\n"))
  } else{
    cat(paste0("Number of function calls by 'optim' = ", x$niter[1],"\n"))
  }
}

#############################################################################################
#' @name residuals.mixpoissonreg
#' @title Residuals Method for \code{mixpoissonreg} Objects
#' @description Function to return 'pearson' or 'score' residuals for mixed Poisson regression models.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param type the type of residual to be returned. Currently, the options are 'pearson' or 'score'. The default is set to 'pearson'. Notice that these
#' residuals coincide for Negative-Binomial models.
#' @param ... Currently not used.
#' @seealso
#' \code{\link{plot.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}},
#' \code{\link{autoplot.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' residuals(daysabs_fit)
#'
#' #Score residuals:
#' residuals(daysabs_fit, type = "score")
#' }
#' @export
residuals.mixpoissonreg <- function(object, type = c("pearson", "score"), ...) {
  type <- rlang::arg_match(type)

  if(!(type%in%c("pearson", "score"))){
    stop("the type must be 'pearson' or 'score'")
  }
if(is.null(object$y)){
  y = object$residuals + object$fitted.values
} else{
  y = object$y
}
  mu <- object$fitted.values
  phi <- object$fitted.precisions

  ddB <- switch(object$modeltype,
                "NB" = {function(x) {1/(x^2)}},
                "PIG" = {function(x) {1/( (-2*x)^(3/2))}}
  )

  qsi0 <- switch(object$modeltype,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  lambda <- lambda_r(y, mu, phi, object$modeltype)


  residuals <- switch(type,
                      "pearson" = {
                      (y-mu)/sqrt(mu*(1+ddB(qsi0)*mu/phi))
                        },
                      "score" = {
                        switch(object$modeltype,
                               "NB" = {(y-mu*lambda)/sqrt(mu -mu^2*ddB(qsi0)/phi + mu^3*ddB(qsi0)/(phi*(mu+phi)))},
                               "PIG" = {y - mu*lambda}
                        )
                      }
                      )
  residuals <- c(residuals)
  names(residuals) <- 1:length(residuals)
  residuals
}


#############################################################################################
#' @name logLik.mixpoissonreg
#' @title logLik Method for \code{mixpoissonreg} Objects
#' @description Function to compute the log-likelihood at the estimated parameters for mixed Poisson regression models.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{vcov.mixpoissonreg}}
#' @examples
#' \donttest{
#' data("Attendance", package = "mixpoissonreg")
#'
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' logLik(daysabs_fit)
#' }
#' @export
logLik.mixpoissonreg <- function(object, ...){
  logLik <- object$logLik
  attr(logLik,"df") = length(object$coefficients$mean) + length(object$coefficients$precision)
  class(logLik) = "logLik"
  logLik
}

#############################################################################################
#' @name coeftest.mixpoissonreg
#' @title coeftest and coefci methods for \model{mixpoissonreg} models
#' @aliases coeftest.mixpoissonreg coefci.mixpoissonreg
#' @description Functions to perform z Wald tests of estimated coefficients and to compute the corresponding Wald confidence intervals.
#' @param x an object of class "mixpoissonreg" containing results from the fitted model.
#' @param vcov. \code{NULL}. The vcov matrix from the \code{mixpoissonreg} object.
#' @param df \code{Inf}, which means that \eqn{z} tests are performed.
#' @param level the confidence level required.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param ... further arguments passed to or from other methods.
#' @noRd
#' @export
coeftest.mixpoissonreg <- function(x, vcov. = NULL, df = Inf, ...){
  coeftest.default(x = x, vcov. = vcov., df = df, ...)
}

#' @noRd
#' @export

coefci.mixpoissonreg <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = Inf, ...){
  coefci.default(x = x, parm = parm, level = level, vcov. = vcov., df = df, ...)
}
