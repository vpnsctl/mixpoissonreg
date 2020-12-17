#############################################################################################
#' @title plot.mixpoissonreg
#' @description Function to build useful plots for mixed Poisson regression models.
#' @param x object of class "mixpoissonreg" containing results from the fitted model.
#' If the model is fitted with envelope = 0, the Q-Q plot will be produced without envelopes.
#' @param which a number of a vector of numbers between 1 and 4. Plot 1:
#' Residuals vs. Index; Plot 2: Q-Q Plot (if the fit contains simulated envelopes,
#' the plot will be with the simulated envelopes); Plot 3: Fitted means vs. Response;
#' Plot 4: Residuals vs. Fitted means.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param main character; title to be placed at each plot additionally (and above) all captions.
#' @param qqline logical; if \code{TRUE} and the fit does *not* contain simulated
#' envelopes, a qqline will be added to the normal Q-Q plot.
#' @param ... graphical parameters to be passed.
#' @seealso
#' \code{\link{summary.mixpoissonreg}}, \code{\link{coef.mixpoissonreg}},
#' \code{\link{vcov.mixpoissonreg}}, \code{\link{fitted.mixpoissonreg}},
#' \code{\link{predict.mixpoissonreg}}
#' @examples
#' \donttest{
#' n <- 100
#' x <- cbind(rbinom(n, 1, 0.5), runif(n, -1, 1))
#' v <- runif(n, -1, 1)
#' z <- simdata_bes(
#'   kap = c(1, 1, -0.5), lam = c(0.5, -0.5), x, v, repetitions = 1,
#'   link.mean = "logit", link.precision = "log"
#' )
#' z <- unlist(z)
#' fit <- bbreg(z ~ x | v, envelope = 10)
#' plot(fit)
#' plot(fit, which = 2)
#' plot(fit, which = c(1, 4), ask = FALSE)
#' }
#' @export
plot.mixpoissonreg <- function(x, which = c(1, 2, 3, 4), ask = TRUE, main = "", qqline = TRUE,
                                include.modeltype = TRUE, ...) {
  if (length(which) == 1) {
    ask <- FALSE
  }

  call_name <- switch(x$estimation_method,
                      "EM" = {"mixpoissonreg"},
                      "ML" = {"mixpoissonregML"}
)

  current_ask = grDevices::devAskNewPage()

  grDevices::devAskNewPage(ask = ask)
  res <- residuals(x, type = x$residualname)
  call_mod <- deparse(x$call)
  name <- x$modelname
  residualname <- paste0(toupper(substring(x$residualname, 1, 1)), substring(x$residualname, 2))
  residualname <- paste0(residualname, " residuals")
  mu_est <- stats::fitted(x, type = "response")

  if (1 %in% which) {
    # First plot (residuals vs index)
    ylab <- residualname
    xlab <- paste0("Index\n ",call_name,"(", call_mod, ", model = \"",x$modeltype, "\")")
    title_1 <- paste0(residualname, " vs Index - ", name)
    graphics::plot(res, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(0, 0, lty = 3)
    graphics::mtext(title_1, side = 3)
  }

  # Second plot (QQ-plot)
  if (2 %in% which) {
    env <- x$envelope
    n <- length(res)
    residualname <- paste0(toupper(substring(x$residualname, 1, 1)), substring(x$residualname, 2))
    ylim <- range(res, env)
    xlim <- c(stats::qnorm(0.5 / n), stats::qnorm(1 - 0.5 / n))
    xlab <- paste0("Theoretical quantiles\n ",call_name,"(", call_mod, ", model = \"",x$modeltype, "\")")
    ylab <- paste0(residualname, " residuals")
    if (is.null(env) == FALSE) {
      title_2 <- paste0("Q-Q Plot with simulated envelopes - ", name)
    } else {
      title_2 <- paste0("Q-Q Plot - ", name)
    }
    RR <- stats::qqnorm(res, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main = main, ...)
    if (is.null(env) == FALSE) {
      aux <- sort(RR$x)
      graphics::lines(aux, env[1, ], col = grDevices::rgb(0.7, 0.7, 0.7))
      graphics::lines(aux, env[3, ], col = grDevices::rgb(0.7, 0.7, 0.7))
      graphics::polygon(c(aux, rev(aux)), c(env[3, ], rev(env[1, ])), col = grDevices::rgb(0.7, 0.7, 0.7), border = NA)
      graphics::lines(aux, env[2, ], lty = 2, lwd = 2)
    } else {
      if (qqline) {
        stats::qqline(res, lty = 3)
      }
    }
    graphics::points(RR$x, RR$y, ...)
    graphics::mtext(title_2, side = 3)
  }

  # Third plot (fitted vs response)

  if (3 %in% which) {
    obs <- x$y

    title_3 <- paste0("Response vs Fitted means - ", name)
    ylab <- "Response"
    xlab <- paste0("Predicted values\n ",call_name,"(", call_mod, ", model = \"",x$modeltype, "\")")
    graphics::plot(mu_est, obs, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(0, 1, lty = 3)
    graphics::mtext(title_3, side = 3)
  }


  # Fourth plot

  if (4 %in% which) {
    title_4 <- paste0(residualname, " vs Fitted means - ", name)
    ylab <- paste0(residualname, " residuals")
    xlab <- paste0("Predicted values\n ",call_name,"(", call_mod, ", model = \"",x$modeltype, "\")")
    graphics::plot(mu_est, res, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(0, 0, lty = 3)
    graphics::mtext(title_4, side = 3)
  }

  grDevices::devAskNewPage(ask = current_ask)
}


#############################################################################################
#' @title fitted.mixpoissonreg
#' @description Function providing the fitted means for the mixed Poisson regression model.
#' @param object object of class "mixpoissonreg" containing results from the fitted model.
#' @param type the type of variable to get the fitted values. The default is the "response" type, which provided the estimated values for the means. The type "link" provides the estimates for the linear predictor of the mean. The type "precision" provides estimates for the precision parameters whereas the type "variance" provides estimates for the variances.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{predict.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}},
#' \code{\link{coef.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' fitted(fit)
#' fitted(fit, type = "precision")
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
#' @title predict.mixpoissonreg
#' @description Function to obtain various predictions based on the fitted mixed Poisson regression model.
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
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}},
#' \code{\link{coef.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' predict(fit)
#' new_data_example <- data.frame(priming = c(0, 0, 1), eliciting = c(0, 1, 1))
#' predict(fit, new_data = new_data_example)
#' predict(fit, new_data = new_data_example, type = "precision")
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

  if(length(interval) > 1){
    interval = interval[1]
  }

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
                              std_error <- sqrt(diag(x_matrix%*%vcov(fit, parameters = "mean")%*%t(x_matrix)))
                              pred_matrix <- switch(type,
                                                    "response" = {
                                                      link_mean <- build_links_mpreg(fit$link.mean)
                                                      lwr_bound <- link_mean$linkinv(stats::fitted(fit, type = 'link') + qnorm( (1-level)/2 )*std_error)
                                                      upr_bound <- link_mean$linkinv(stats::fitted(fit, type = 'link') + qnorm( (1+level)/2 )*std_error)
                                                      cbind(stats::fitted(fit, type = "response"), lwr_bound, upr_bound)
                                                    },
                                                    "link" = {
                                                      lwr_bound <- stats::fitted(fit, type = 'link') + qnorm( (1-level)/2 )*std_error
                                                      upr_bound <- stats::fitted(fit, type = 'link') + qnorm( (1+level)/2 )*std_error
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
                                                      std_error_mean <- sqrt(diag(x_matrix%*%vcov(object, parameters = "mean")%*%t(x_matrix)))
                                                      std_error_precision <- sqrt(diag(w_matrix%*%vcov(object, parameters = "precision")%*%t(w_matrix)))

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
                                                                                                  rpois(nsim_pred_y,ig*mu_pred[j])
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
        std_error <- sqrt(diag(x_matrix%*%vcov(object, parameters = "mean")%*%t(x_matrix)))
      }
      predictions_temp$se.fit <- std_error
      predictions <- predictions_temp
    }
  } else {
    formula_temp <- Formula(fit$call)
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
                              std_error <- sqrt(diag(matrix_temp_x%*%vcov(object, parameters = "mean")%*%t(matrix_temp_x)))
                              pred_matrix <- switch(type,
                                                    "response" = {
                                                      link_mean <- build_links_mpreg(fit$link.mean)
                                                      lwr_bound <- link_mean$linkinv(matrix_temp_x %*% beta_est + qnorm( (1-level)/2 )*std_error)
                                                      upr_bound <- link_mean$linkinv(matrix_temp_x %*% beta_est + qnorm( (1+level)/2 )*std_error)
                                                      cbind(predictions, lwr_bound, upr_bound)
                                                    },
                                                    "link" = {
                                                      lwr_bound <- predictions + qnorm( (1-level)/2 )*std_error
                                                      upr_bound <- predictions + qnorm( (1+level)/2 )*std_error
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
                                                      std_error_mean <- sqrt(diag(matrix_temp_x%*%vcov(object, parameters = "mean")%*%t(matrix_temp_x)))
                                                      std_error_precision <- sqrt(diag(matrix_temp_w%*%vcov(object, parameters = "precision")%*%t(matrix_temp_w)))

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
                                                                                rpois(nsim_pred_y,ig*mu_pred[j])
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
        std_error <- sqrt(diag(matrix_temp_x%*%vcov(object, parameters = "mean")%*%t(matrix_temp_x)))
      }
      predictions_temp$se.fit <- std_error
      predictions <- predictions_temp
    }

  }
  return(predictions)
}


#############################################################################################
#' @title print.mixpoissonreg
#' @description Function providing a brief description of results related to the mixed Poisson regression model.
#' @param x object of class "mixpoissonreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}},
#' \code{\link{coef.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' fit
#' }
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
  cat("Call:", "\n")
  call_model <- deparse(x$call)
  cat(call_model)
  cat("\n\n")
  cat(paste0("Coefficients modeling the mean (with ", x$link.mean, " link):", "\n"))
  print(coeff_beta)
  cat(paste0("Coefficients modeling the precision (with ", x$link.precision, " link):", "\n"))
  print(coeff_alpha)
}


#############################################################################################
#' @title coef.mixpoissonreg
#' @description Function to extract the coefficients of a fitted mixed Poisson regression model.
#' @param object object of class "mixpoissonreg" containing results from the fitted model.
#' @param parameters a string to determine which coefficients should be extracted: 'all' extracts all coefficients, 'mean' extracts the coefficients of the mean parameters and 'precision' extracts coefficients of the precision parameters.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{summary.mixpoissonreg}},
#' \code{\link{print.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' coef(fit)
#' coef(fit, parameters = "precision")
#' }
#' @export
coef.mixpoissonreg <- function(object, parameters = c("all", "mean", "precision"), ...) {
  if(length(parameters)>1){
    parameters = parameters[1]
  }
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
#' @title vcov.mixpoissonreg
#' @description Function to extract the variance-covariance matrix of the parameters of the fitted mixed Poisson regression model.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param parameters a string to determine which coefficients should be extracted: 'all' extracts all coefficients, 'mean' extracts the coefficients of the mean parameters and 'precision' extracts coefficients of the precision parameters.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{obs_fisher_information_mixpoisson}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
#' vcov(fit)
#' vcov(fit, parameters = "precision")
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
#' @title summary.mixpoissonreg
#' @description Function providing a summary of results related to the mixed Poisson regression model.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{coef.mixpoissonreg}},
#' \code{\link{print.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
#' summary(fit)
#' }
#' @export
summary.mixpoissonreg <- function(object, ...) {
  nbeta <- length(object$coefficients$mean)
  nalpha <- length(object$coefficients$precision)
  #
  coeff <- c(coef(object)$mean, coef(object)$precision)
  SEr <- object$std_errors
  tab <- cbind(coeff, SEr, coeff / SEr, 2 * stats::pnorm(-abs(coeff / SEr)))
  colnames(tab) <- c("Estimate", "Std.error", "z-value", "Pr(>|z|)")
  rownames(tab) <- names(coeff)
  tab <- list(mean = tab[seq.int(length.out = nbeta), , drop = FALSE], precision = tab[seq.int(length.out = nalpha) + nbeta, , drop = FALSE])
  #
  digits <- max(3, getOption("digits") - 3)
  #

  call_name <- switch(object$estimation_method,
                      "EM" = {"mixpoissonreg"},
                      "ML" = {"mixpoissonregML"}
  )

  optim_algo <- switch(object$estimation_method,
                       "EM" = {"Expectation-Maximization Algorithm"},
                       "ML" = {"Maximum-Likelihood Estimation"}
  )
  cat("\n")
  cat(paste0(object$modelname, " - ", optim_algo))

  cat("\n\n")
  cat("Call:", "\n")
  call_model <- deparse(object$call)
  cat(call_model)
  cat("\n")


  #
  RSS <- sum(object$residuals^2)
  residualname <- paste0(toupper(substring(object$residualname, 1, 1)), substring(object$residualname, 2))
  cat(sprintf("\n%s:\n", paste0(residualname, " residuals")))
  print(structure(round(c(RSS, as.vector(stats::quantile(object$residuals))), digits = digits), .Names = c("RSS", "Min", "1Q", "Median", "3Q", "Max")))
  #
  if (NROW(tab$mean)) {
    cat(paste0("\nCoefficients modeling the mean (with ", object$link.mean, " link):\n"))
    stats::printCoefmat(tab$mean, digits = digits, signif.legend = FALSE)
  } else {
    message("\nNo coefficients modeling the mean. \n")
  }
  #
  if (NROW(tab$precision)) {
    cat(paste0("\nCoefficients modeling the precision (with ", object$link.precision, " link):\n"))
    stats::printCoefmat(tab$precision, digits = digits, signif.legend = FALSE)
  } else {
    message("\nNo coefficients modeling the precision. \n")
  }
  #
  if (getOption("show.signif.stars")) {
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n\n")
  }
  #

  cat("Efron's pseudo R-squared: ", object$efron.pseudo.r2,"\n")

  if (is.null(object$envelope) == FALSE) {
    message(sprintf("%s\n", paste0("Percentage of residual within the envelope = ", round(object$envelope_prop, digits = digits))))
  }
  if(object$estimation_method == "EM"){
    cat(paste0("Number of iterations of the EM algorithm = ", object$niter,"\n"))
  } else{
    cat(paste0("Number of function calls by 'optim' = ", object$niter[1],"\n"))
  }
}

#############################################################################################
#' @title residuals.mixpoissonreg
#' @description Function to return 'pearson' or 'score' residuals for mixed Poisson regression model.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param type the type of residual to be returned. Currently, the options are 'pearson' or 'score'. The default is set to 'pearson'. Notice that these
#' residuals coincide for Negative-Binomial models.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{coef.mixpoissonreg}},
#' \code{\link{print.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
#' summary(fit)
#' }
#' @export
residuals.mixpoissonreg <- function(object, type = c("pearson", "score")) {
  if(length(type)>1){
    type = type[1]
  }
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
#' @title logLik.mixpoissonreg
#' @description Function to compute the log-likelihood at the estimated parameters for mixed Poisson regression model.
#' @param object an object of class "mixpoissonreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.mixpoissonreg}}, \code{\link{coef.mixpoissonreg}},
#' \code{\link{print.mixpoissonreg}}, \code{\link{vcov.mixpoissonreg}},
#' \code{\link{plot.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
#' summary(fit)
#' }
#' @export
logLik.mixpoissonreg <- function(object){
  logLik <- object$logLik
  attr(logLik,"df") = length(object$coefficients$mean) + length(object$coefficients$precision)
  class(logLik) = "logLik"
  logLik
}
