#' @import statmod

#############################################################################################
#' @name startvalues_mpreg
#' @title Initial guesses for mixed Poisson regression models
#' @description Function providing initial values for the optimization algorithms for mixed Poisson regression models.
#' @param y response vector with z_i >= 0 and integers.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" or "inverse.sqrt".
#' @return a list containing
#' \itemize{
#'   \item beta - the initial guesses for the mean-related parameters;
#'   \item alpha - the initial guesses for the precision-relared parameters.
#'   }
#'   @noRd
startvalues_mpreg <- function(y, x, w, link.mean, link.precision, model) {
  nbeta <- ncol(x)
  nalpha <- ncol(w)
  n <- length(y)
  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  fit_aux <- stats::glm.fit(x = x, y = y, family = stats::quasipoisson(link = link.mean))
  beta_start <- fit_aux$coefficients

  mu_est <- link_mean$linkinv(x %*% beta_start)

  disp_temp <- sum( ((y - mu_est)^2 / mu_est ) )/(n - nbeta)



  phi_est <- 1 / (disp_temp - 1)

  if(phi_est <= 0){
    phi_est <- 1
  }

  alpha_start <- c(link_precision$linkfun(phi_est), rep(0, nalpha - 1))

  out <- list(beta = beta_start, alpha = alpha_start)

  return(out)
}


#############################################################################################
#' @name build_links_mpreg
#' @title Build link functions
#' @description Function to construct the link functions.
#' @param link the name of the link function.
#' @return A link function object of the class \code{link-glm}.
#' @noRd
build_links_mpreg <- function(link) {
  if(link %in% c("log", "sqrt", "identity")){
    return(stats::make.link(link))
  } else if(link == "inverse.sqrt"){
    linkfun <- function(y){1/sqrt(y)}
    ## inverse link
    linkinv <- function(eta){1/(eta^2)}
    ## derivative of invlink wrt eta
    mu.eta <- function(eta) { -2/eta^3 }
    valideta <- function(eta) TRUE
    name.link <- "inverse.sqrt"
    inverse_sqrt = structure(list(linkfun = linkfun,
                                  linkinv = linkinv,
                                  mu.eta = mu.eta,
                                  valideta = valideta,
                                  name = name.link),
              class = "link-glm")
    return(inverse_sqrt)
  } else{
    stop("link must be 'log', 'sqrt', 'identity' or 'inverse.sqrt'.")
  }
}



#############################################################################################
#' @name d2mudeta2
#' @title Second derivative of the mean with respect to the linear predictor.
#' @description Function to obtain the second derivatives of the mean parameter with respect to the linear predictor.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" or "sqrt".
#' @param mu mean parameter.
#' @noRd

d2mudeta2 <- function(link.mean, mu) {
  d2mu <- switch(link.mean,
                 "log" = {
                   mu
                 },
                 "sqrt" = {
                   2
                 }
  )
  return(d2mu)
}

#############################################################################################
#' @name d2phideta2
#' @title Second derivative of the precision parameter with respect to the linear predictor.
#' @description Function to obtain the second derivatives of the precision parameter with respect to the linear predictor.
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" or "inverse.sqrt".
#' @param phi precision parameter.
#' @noRd

d2phideta2 <- function(link.precision, phi) {
  d2phi <- switch(link.precision,
                  "identity" = {
                    0
                  },
                  "log" = {
                    phi
                  },
                  "inverse.sqrt" = {
                    6/phi^4
                  }
  )
  return(d2phi)
}





#############################################################################################
#' @name generate_data_mixpoisson
#' @title Auxiliar function to generate mixed Poisson random numbers
#' @description Function to generate synthetic data from mixed Poisson regression models.
#' @param coefficients a list containing elements 'beta' and 'alpha'. 'beta' and 'alpha' are vectors of estimated parameters of the model.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param repetitions the number of random draws to be made.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return a list of response vectors y, with y being nonnegative integers.
#' @noRd

generate_data_mixpoisson <- function(coefficients, x, w,
                      repetitions, link.mean,
                      link.precision,model){
  x = as.matrix(x)
  w = as.matrix(w)

  ncolx <- ncol(x)
  ncolw <- ncol(w)


  beta <- coefficients$mean
  alpha <- coefficients$precision

  nbeta <- length(beta)
  nalpha <- length(alpha)

  if ((nbeta - ncolx) == 1) {
    x <- cbind(1, x)
  }
  if ((nalpha - ncolw) == 1) {
    w <- cbind(1, w)
  }

  if (abs(nbeta - ncolx) > 1) {
    stop("check dimension of beta and x")
  }
  if (abs(nalpha - ncolw) > 1) {
    stop("check dimension of alpha and w")
  }
  #
  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)
  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)
  n <- length(mu)

  y <- lapply(1:repetitions, function(x) {
    y_temp <- switch(model,
                "NB" = {stats::rnbinom(n = n, mu = mu, size = phi)},
                "PIG" = {ig = rinvgauss(n = n,mean=1,dispersion=1/phi)
                  y  = rpois(n = n, lambda = ig*mu)
                  y}
                )
    return(y_temp) #
  })
  return(y)
}


#############################################################################################
#' @name envelope_mixpoisson
#' @title Function to compute simulated envelopes.
#' @description Function to calculate envelopes based on residuals for the mixed Poisson regression models.
#' @param residual character indicating the type of residual ("pearson" or "score").
#' @param estimation_method character indicating the estimation method ("EM" or "ML").
#' @param coefficients a list containing elements 'mean' and 'precision'. 'mean' and 'precision' are vectors of estimated parameters of the model.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_env number of synthetic data sets to be generated.
#' @param prob confidence level of the envelope (number between 0 and 1).
#' @param n sample size.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @param em_controls only used with the 'EM' method. A list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000;
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5); \code{em_tolgrad} that defines the tolerance value
#' of the maximum-norm of the the gradient of the Q-function, the default is set to 10^(-2).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function.
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @noRd
envelope_mixpoisson <- function(residual, estimation_method,
                    coefficients, x, w, nsim_env,
                    prob, n, link.mean,
                    link.precision,model, em_controls, optim_method, optim_controls) {

  ysim <- generate_data_mixpoisson(coefficients, x, w, nsim_env, link.mean, link.precision, model)

  residuals_envelope <- switch(residual,
                "pearson" = {
                  residuals_env <- pblapply(ysim, function(y_env) {
                    coeff_env <- switch(estimation_method,
                           "EM" = {EM_mixpoisson(coefficients$mean, coefficients$precision, y_env, x, w, link.mean, link.precision,
                                                 model,  em_controls, optim_method, optim_controls)$coefficients},
                           "ML" = {ML_mixpoisson(coefficients$mean, coefficients$precision, y_env, x, w,
                                                 link.mean, link.precision, model, optim_method, optim_controls)$coefficients}
                           )
                    pearson_res <- pearson_residual_mixpoisson(coeff_env, y_env, x, w, link.mean, link.precision, model)
                    pearson_res
                  })
                  residuals_env
                },
                "score" = {
                  residuals_env <- pblapply(ysim, function(y_env) {
                    coeff_env <- switch(estimation_method,
                                        "EM" = {EM_mixpoisson(coefficients$mean, coefficients$precision, y_env, x, w, link.mean, link.precision,
                                                              model,  em_controls, optim_method, optim_controls)$coefficients},
                                        "ML" = {ML_mixpoisson(coefficients$mean, coefficients$precision, y_env, x, w,
                                                              link.mean, link.precision, model, optim_method, optim_controls)$coefficients}
                    )
                    score_res <- score_residual_mixpoisson(coeff_env, y_env, x, w, link.mean, link.precision, model)
                    score_res
                  })
                  residuals_env
                  })
  residuals_envelope <- t(matrix(unlist(residuals_envelope), n, nsim_env))
  residuals_envelope <- t(apply(residuals_envelope, 1, sort))
  residuals_envelope <- apply(residuals_envelope, 2, sort)
  id1 <- max(1, round(nsim_env * (1 - prob) / 2))
  id2 <- round(nsim_env * (1 + prob) / 2)
  envelopes <- rbind(residuals_envelope[id2, ], apply(residuals_envelope, 2, median), residuals_envelope[id1, ])
  rownames(envelopes) <- c("upper", "median", "lower")
  return(envelopes)
}

#############################################################################################
#' @name pearson_residual_mixpoisson
#' @title Pearson residuals for mixed Poisson regression models
#' @description Function to calculate the pearson residuals for mixed Poisson regression models.
#' @param coefficients a list containing elements 'beta' and 'alpha'. 'beta' and 'alpha' are vectors of estimated parameters of the model.
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @noRd
pearson_residual_mixpoisson <- function(coefficients, y, x, w,
                            link.mean, link.precision, model){
  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  beta_est <- coefficients$mean
  alpha_est <- coefficients$precision

  mu_hat = link_mean$linkinv(x%*%beta_est)
  phi_hat = link_precision$linkinv(w%*%alpha_est)

  ddB <- switch(model,
         "NB" = {function(x) {1/(x^2)}},
         "PIG" = {function(x) {1/( (-2*x)^(3/2))}}
         )

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
                 )

  res.pearson <- (y-mu_hat)/sqrt(mu_hat*(1+ddB(qsi0)*mu_hat/phi_hat))
  res.pearson <- c(res.pearson)
  names(res.pearson) <- 1:length(res.pearson)
  return(res.pearson)
}


#############################################################################################
#' @name score_residual_mixpoisson
#' @title Score residuals for mixed Poisson regression models
#' @description Function to calculate the score residuals for mixed Poisson regression models.
#' @param coefficients a list containing elements 'beta' and 'alpha'. 'beta' and 'alpha' are vectors of estimated parameters of the model.
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return Vector containing the score residuals. Notice that the score residuals for the PIG regression
#' models are not variance-normalized, whereas the score residuals for the NB regression models are
#' variance-normalized. See Barreto-Souza and Simas (2015) for further details.
#' @noRd


score_residual_mixpoisson <- function(coefficients, y, x, w,
                            link.mean, link.precision, model){
  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  beta_est <- coefficients$mean
  alpha_est <- coefficients$precision

  mu_hat = link_mean$linkinv(x%*%beta_est)
  phi_hat = link_precision$linkinv(w%*%alpha_est)

  ddB <- switch(model,
                "NB" = {function(x) {1/(x^2)}},
                "PIG" = {function(x) {1/( (-2*x)^(3/2))}}
  )

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  lambda <- lambda_r(y, mu_hat, phi_hat, model)

  res.score <- switch(model,
                      "NB" = {(y-mu_hat*lambda)/sqrt(mu_hat -mu_hat^2*ddB(qsi0)/phi_hat + mu_hat^3*ddB(qsi0)/(phi_hat*(mu_hat+phi_hat)))},
                      "PIG" = {y - mu_hat*lambda}
  )

  res.score <- c(res.score)
  names(res.score) <- 1:length(res.score)
  return(res.score)
}

#############################################################################################
#' @name std_error_mixpoisson
#' @title Standard errors of coefficients of mixed Poisson regression models.
#' @description Function to calculate the standard errors of the coefficients of a fitted mixed Poisson regression model.
#' @param coefficients a list containing elements 'beta' and 'alpha'. 'beta' and 'alpha' are vectors of estimated parameters of the model.
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return Vector containing the score residuals.
#' @noRd
std_error_mixpoisson <- function(coefficients, y, x, w,
                link.mean, link.precision, model){
  theta = c(coefficients$mean, coefficients$precision)
  obs_inf_fisher <- obs_fisher_information_mixpoisson(theta, y, x, w,
                                                link.mean, link.precision, model)
  inv_fisher <- tryCatch(solve(obs_inf_fisher), error = function(e) matrix(NA,
                                                                        nrow(obs_inf_fisher), col(obs_inf_fisher)))
  sqrt(diag(inv_fisher))
}


#############################################################################################
#' @name lambda_r
#' @title Auxiliary function for the E-step
#' @description Auxiliary function given by the conditional expected value E(Z_i|Y; theta), where
#' Z is the latent variable (for further details see Barreto-Souza and Simas, 2015).
#' @param y nonnegative integers response vector.
#' @param mu mean parameter (vector having the same size of z).
#' @param phi precision parameter (vector having the same size of z).
#' @param model the mixed Poisson regression model, "NB" or "PIG".
#' @return Vector of expected values.
#' @noRd


lambda_r <- function(y, mu, phi, model) {
  out <- switch(model,
                "NB" = {(y+phi)/(mu+phi)},
                "PIG" ={(y+1)*dPIG(y+1,mu,1/phi)/(mu*dPIG(y,mu,1/phi))}
  )
  return(out)
}


#############################################################################################
#' @name kappa_r
#' @title Auxiliary function for the E-step
#' @description Auxiliary function given by the conditional expected value E(g(Z_i)|Y; theta), where
#' Z is the latent variable and the function g() is obtained in the decomposition of the function
#' c( , ) of the expression of the density of Z in the exponential family (for further details see Barreto-Souza and Simas, 2015).
#' @param y nonnegative integers response vector.
#' @param mu mean parameter (vector having the same size of z).
#' @param phi precision parameter (vector having the same size of z).
#' @param model the mixed Poisson regression model, "NB" or "PIG".
#' @return Vector of expected values.
#' @noRd

kappa_r <- function(y, mu, phi, model) {
  out <- switch(model,
                "NB" = {digamma(y+phi)-log(mu+phi)},
                "PIG" ={res<-(y==0)*((1+sqrt(phi*(2*mu+phi)))/phi)+
                  (y>0)*(mu*dPIG((y>0)*(y-1),mu,1/phi)/
                           (((y==0)+y)*dPIG(y,mu,1/phi)))
                -.5*res}
  )
  return(out)
}


#############################################################################################
#' @name update.mixpoissonreg
#' @title Update method for \code{mixpoissonreg} objects
#' @description Update and re-fit a model.
#' @param object a fitted \code{mixpoissonreg} object.
#' @param formula changes to the formula;
#' @param ... Additional arguments to the call, or arguments with changed values. Use name = NULL to remove the argument name.
#' @param evaluate If true evaluate the new call else return the call.
#' @return If evaluate = TRUE the fitted object, otherwise the updated call.
#' @noRd
update.mixpoissonreg <- function (object, formula., ..., evaluate = TRUE)
{
  if (is.null(call <- getCall(object)))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- formula(update(Formula(formula(object)),
                                   formula.))
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, parent.frame())
  else call
}
