#' @import gamlss.dist

#############################################################################################
#' @title ML_mixpoisson
#' @description Function to run the obtain maximum-likelihood estimates of a mixed Poisson regression model.
#' @param beta initial values for the mean-related coefficients.
#' @param alpha initial values for the precision-related coefficients.
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param epsilon tolerance to control the convergence criterion.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function.
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return a list containing
#' \itemize{
#'     \item niter - number of \code{optim} iterations;
#'     \item coefficients - a list containing estimated vectors 'beta' (mean-related coefficients) and 'alpha' (precision-related coefficients);
#'     \item fitted.values - the estimated means.
#' }
ML_mixpoisson <- function(beta, alpha, y, x, w,
                          link.mean, link.precision, model, optim_method, optim_controls){
  nbeta <- length(beta)

  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  theta <- c(beta, alpha)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  M <- tryCatch(
    stats::optim(
    par = theta,
    fn = loglik_mixpoisson,
    gr = score_mixpoisson,
    y = y,
    x = x,
    w = w,
    link.mean = link.mean,
    link.precision = link.precision,
    model = model,
    control = c(fnscale = -1, optim_controls),
    method = optim_method
  ), error = function(e) {
    "Error"
  })

  if (length(M) == 1) {
    warning("Trying with numerical derivatives")
    M <- tryCatch(
      stats::optim(
        par = theta,
        fn = loglik_mixpoisson,
        gr = NULL,
        y = y,
        x = x,
        w = w,
        link.mean = link.mean,
        link.precision = link.precision,
        model = model,
        control = c(fnscale = -1, optim_controls),
        method = optim_method
      ), error = function(e) {
        "Error"
      })
  }

  if (length(M) == 1) {
    warning("Trying with another optimization algorithm")
    if(optim_method == "L-BFGS-B"){
      optim_temp = "Nelder-Mead"
    } else if(optim_method == "Nelder-Mead"){
      optim_temp = "L-BFGS-B"
    } else {
      optim_temp = "Nelder-Mead"
    }
    M <- tryCatch(
      stats::optim(
        par = theta,
        fn = loglik_mixpoisson,
        gr = NULL,
        y = y,
        x = x,
        w = w,
        link.mean = link.mean,
        link.precision = link.precision,
        model = model,
        control = c(fnscale = -1, optim_controls),
        method = optim_temp
      ), error = function(e) {
        "Error"
      })
  }

  if (length(M) == 1) {
    stop("The algorithm did not converge")
  }

  theta <- M$par

  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  if (sum(phi <= 0) > 0) {
    warning("one or more estimates of precision parameters were negative. Please,
         consider using another link function for the precision parameter.")
  }

  if (sum(is.nan(phi)) > 0) {
    warning("one or more estimates of precision parameters were not a number. Please,
         consider using another link function for the precision parameter.")
  }

  out <- list()
  out$coefficients <- list(mean = beta, precision = alpha)
  out$fitted.values <- mu
  out$fitted.precisions <- phi
  out$niter <- M$counts
  return(out)
}


#############################################################################################
#' @title loglik_mixpoisson
#' @description log-likelihood function of the mixed Poisson regression model.
#' @param theta vector of parameters (all coefficients).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return scalar representing the value of the log-likelihood function at 'theta'.
loglik_mixpoisson <- function(theta, y, x, w,
                              link.mean, link.precision,
                              model) {
  nbeta <- ncol(x)
  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]
  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  loglik = switch(model,
                  "PIG" = {sum(dPIG(y, mu = mu, sigma = 1/phi, log = TRUE))},
                  "NB" = {sum(dNBI(y, mu = mu, sigma = 1/phi, log = TRUE))}
                  )
  loglik
}

#############################################################################################
#' @title score_mixpoisson
#' @description Function to calculate the score vector of a mixed Poisson regression model.
#' @param theta vector of parameters (all coefficients).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return vector containing the gradient of the Q function at theta.
score_mixpoisson <- function(theta, y, x, w, link.mean, link.precision, model) {
  nbeta <- ncol(x)
  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]
  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  dmudeta <- link_mean$mu.eta(x %*% beta)
  dphideta <- link_precision$mu.eta(w %*% alpha)

  grad1 <- switch(model,
                  "PIG" = {t(x)%*% ( dmudeta*(y/mu - lambda_r(y,mu,phi,"PIG")))},
                  "NB" = {t(x) %*% (dmudeta * phi * ( (y-mu)/(mu*(mu + phi))))}
                  )

  grad2 <- switch(model,
                  "PIG" = {t(w)%*% (dphideta * ( (1+y/phi) - (1 + mu/phi) * lambda_r(y,mu,phi,"PIG") ))},
                  "NB" = {
                    t(w) %*% (dphideta * (mu/(mu+phi) - log(mu/phi+1) + mu*y/(phi*mu + phi^2) -y/phi +digamma(y+phi) -digamma(phi)))}
                  )
  c(grad1, grad2)
}
