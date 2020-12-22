#' @importFrom stats cooks.distance hatvalues influence terms

#############################################################################################
#' @name influence.mixpoissonreg
#' @title Global Influence Diagnostics for Mixed Poisson Regression Models
#' @aliases cooks.distance.mixpoissonreg hatvalues.mixpoissonreg influence.mixpoissonreg
#' @description These functions provides global influence diagnostic quantities such as Cook's distance, hat values, generalized Cook's distance (through argument
#' on \code{cooks.distance.mixpoissonreg} function), likelihood displacement (through argument
#' on \code{cooks.distance.mixpoissonreg} function) and Q-displacement (through argument
#' on \code{cooks.distance.mixpoissonreg} function).
#' @param model a \code{mixpoissonreg} object.
#' @param hat hat values \code{H[i,i]}. The default is obtained through the second-derivative of the Q-function in the spirit of Zhu et al. (2001) and Pregibon (1981), see details.
#' @param type the type of Cook's distance to be used. The options are 'CD', the standard Cook's distance;
#' 'GCD', the generalized Cook's distance with respect to all parameters, Zhu et al. (2001);
#' 'GCDmean', the generalized Cook's distance with respect to the mean parameters;
#' 'GCDprecision', the generalized Cook's distance with respect to the precision parameters;
#' 'LD', the likelihood displacement (also known as likelihood distance), see Cook and Weisberg (1982); 'QD', the Q-displacement, see Zhu et al. (2001).
#' See 'details' for more informations. For 'GCD', 'GCDmean', 'LD' and 'QD', the model must be fitted with 'x' set to TRUE,
#' and for 'GCD', 'GCDprecision', 'LD' and 'QD', the model must be fitted with 'w' set to TRUE. For 'CD', if 'hat' is set to 'mean',
#' the model must be fitted with 'x' set to TRUE, whereas
#' if 'hat' is set to 'precision', the model must be fitted with 'w' set to TRUE.
#' @param parameters the parameter to which the hat values will be computed. The options are 'mean' and 'precision'. The default is 'mean'. For hatvalues with respect to the mean
#' the model must be fitted with 'x' set to TRUE, and for hatvalues with respect to the precision the model must be fitted with 'w' set to TRUE.
#' @param do.coef logical indicating if the the approximation to the change of coefficients values after case removal are desired. The model must be fitted with x = TRUE. See details for further explanations.
#' @param ... Currently not used.
#' @details
#' For hatvalues of mixed Poisson regression models, we follow Zhu et al. (2001) to consider the negative of the hessian of the Q-function as weight matrix, and follow
#' Pregibon (1981) to define the 'hat' matrix with respect to this weight matrix. We can consider the hessian of the Q-function with respect to mean-related parameters,
#' which usually considered. We can also consider the hessian of the Q-function with respect to the precision-related parameters to give rise to hatvalues related to the precision
#' parameters.
#'
#' The Generalized Cook's distance and Q-displacement for EM-based models were defined in Zhu et al. (2001) and computed for mixed Poisson regression models
#' in Barreto-Souza and Simas (2015). We implemented first-order approximation to these quantities to be computationally feasible. These first-order approximations
#' are available in Barreto-Souza and Simas (2015). We also provide versions of generalized Cook's distance for mean-related or precision-related parameters, whose
#' details can be found in Barreto-Souza and Simas (2015).
#'
#' In the influence method we provide a 'do.coef' argument that computes first-order approximations to the impact of removal of each case to each parameter, in the
#' same spirit as the 'do.coeff' argument in 'influence.lm'.
#'
#'
#'
#' @references
#' DOI:10.1007/s11222-015-9601-6 (\href{https://doi.org/10.1007/s11222-015-9601-6}{Barreto-Souza and Simas; 2015})
#'
#' Cook, D.R. and Weisberg, S. (1982) *Residuals and Influence in Regression*. (New York: Chapman and Hall, 1982)
#'
#' DOI:10.1214/aos/1176345513 (\href{https://projecteuclid.org/euclid.aos/1176345513}{Pregibon; 1981})
#'
#' Zhu, H.T., Lee, S.Y., Wei, B.C., Zhu, J. (2001) *Case-deletion measures formodels with incomplete data.* Biometrika, 88, 727–737. \href{https://www.jstor.org/stable/2673442?seq=1}{https://www.jstor.org/stable/2673442?seq=1}

#' @rdname influence.mixpoissonreg
#' @export

hatvalues.mixpoissonreg <- function(model, parameters = c("mean", "precision"), ...){
  parameters <- rlang::arg_match(parameters)

  if(!(parameters%in% c("mean", "precision"))){
    stop("the parameter must be 'mean' or 'precision'")
  }
  if(parameters == "mean" & is.null(model$x)){
    stop("x component not found. fit the model again with argument x = TRUE")
  }
  if(parameters == "precision" & is.null(model$w)){
    stop("w component not found. fit the model again with argument w = TRUE")
  }
 W = obs_fisher_weight_matrix_mixpoisson(model, parameters)
 eigen_W <- eigen(W, symmetric = TRUE)
 ev_W <- eigen_W$vectors
 sqrt_W <- ev_W %*% diag(sqrt(eigen_W$values)) %*% t(ev_W)
 H <- switch(parameters,
             "mean" = {x = model$x
               sqrt_W %*% x %*% solve(t(x)%*%W%*%x)%*%t(x)%*% sqrt_W},
             "precision" = {w = model$w
             sqrt_W %*% w %*% solve(t(w)%*%W%*%w)%*%t(w)%*% sqrt_W})
 h <- diag(H)
 names(h) <- 1:length(h)
 h
}

#' @rdname influence.mixpoissonreg
#' @export

cooks.distance.mixpoissonreg <- function(model, type = c("CD", "GCD", "GCDmean", "GCDprecision", "LD", "QD"), hat = c("mean", "precision"), ...){
type <- rlang::arg_match(type)

if(!(type%in%c("CD", "GCD", "GCDmean", "GCDprecision", "LD", "QD"))){
  stop("type must be one of 'CD', 'GCD', 'GCDmean'', 'GCDprecision', 'LD', 'QD'")
}

hat <- rlang::arg_match(hat)

if(!(hat%in%c("mean","precision"))){
  stop("hat must be one of 'mean' or 'precision'")
}

h = switch(hat,
           "mean" = {
             if(is.null(model$x)){
               stop("x component not found. fit the model again with argument x = TRUE")
             }
             stats::hatvalues(model, parameters = hat)
             },
           "precision" = {
             if(is.null(model$w)){
               stop("w component not found. fit the model again with argument w = TRUE")
             }
             stats::hatvalues(model, parameters = hat)
           })
dist <- switch(type,
               "CD" = {
                 p <- length(model$coefficients$mean)
                 pearson_residuals <- stats::residuals(model, type = "pearson")
                 h*(pearson_residuals^2)/(p*(1-h)^2)
               },
               "GCD" = {
                 if(is.null(model$x)){
                   stop("x component not found. fit the model again with argument x = TRUE")
                 }
                 if(is.null(model$w)){
                   stop("w component not found. fit the model again with argument w = TRUE")
                 }
                 x = model$x
                 w = model$w
                 mu = model$fitted.values
                 phi = model$fitted.precisions
                 if(is.null(model$y)){
                   y = model$residuals + model$fitted.values
                 } else{
                   y = model$y
                 }

                 lambda = lambda_r(y,mu,phi,model$modeltype)

                 a = mu*lambda - y

                 link_precision <- build_links_mpreg(model$link.precision)
                 eta_prec <- link_precision$linkfun(phi)
                 dphideta <- link_precision$mu.eta(eta_prec)

                 qsi0 <- switch(model$modeltype,
                                "NB" = {-1},
                                "PIG" = {-1/2}
                 )

                 B<- switch(model$modeltype,
                            "PIG" = {function(x) {-sqrt((-2*x))}},
                            "NB" = {function(x) {-log(-x)}}
                 )

                 dD<-switch(model$modeltype,
                            "PIG" = {function(x) {1/(2*x)}},
                            "NB" = {function(x) {log(x)+1-digamma(x)}}
                 )

                 kappa = kappa_r(y,mu,phi, model$modeltype)

                 b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)
                 Wbeta <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "mean")
                 Walpha <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "precision")
                 gencook = a^2 * diag(x %*% solve(t(x)%*%Wbeta%*%x)%*%t(x)) + b^2 * diag(w%*%solve(t(w)%*%Walpha%*%w)%*%t(w))
                 gencook <- c(gencook)
                 names(gencook) <- 1:length(gencook)
                 gencook
               },
               "GCDmean" = {
                 if(is.null(model$x)){
                   stop("x component not found. fit the model again with argument x = TRUE")
                 }
                 x = model$x
                 phi = model$fitted.precisions
                 if(is.null(model$y)){
                   y = model$residuals + model$fitted.values
                 } else{
                   y = model$y
                 }
                 mu = model$fitted.values
                 phi = model$fitted.precisions
                 lambda = lambda_r(y,mu,phi,model$modeltype)

                 a = mu*lambda - y

                 kappa = kappa_r(y,mu,phi, model$modeltype)

                 Wbeta <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "mean")
                 gencook = a^2 * diag(x %*% solve(t(x)%*%Wbeta%*%x)%*%t(x))
                 gencook <- c(gencook)
                 names(gencook) <- 1:length(gencook)
                 gencook
               },
               "GCDprecision" = {
                 if(is.null(model$w)){
                   stop("w component not found. fit the model again with argument w = TRUE")
                 }
                 w = model$w
                 mu = model$fitted.values
                 phi = model$fitted.precisions
                 if(is.null(model$y)){
                   y = model$residuals + model$fitted.values
                 } else{
                   y = model$y
                 }
                 mu = model$fitted.values
                 phi = model$fitted.precisions
                 lambda = lambda_r(y,mu,phi,model$modeltype)

                 a = mu*lambda - y

                 link_precision <- build_links_mpreg(model$link.precision)
                 eta_prec <- link_precision$linkfun(phi)
                 dphideta <- link_precision$mu.eta(eta_prec)

                 qsi0 <- switch(model$modeltype,
                                "NB" = {-1},
                                "PIG" = {-1/2}
                 )

                 B<- switch(model$modeltype,
                            "PIG" = {function(x) {-sqrt((-2*x))}},
                            "NB" = {function(x) {-log(-x)}}
                 )

                 dD<-switch(model$modeltype,
                            "PIG" = {function(x) {1/(2*x)}},
                            "NB" = {function(x) {log(x)+1-digamma(x)}}
                 )

                 kappa = kappa_r(y,mu,phi, model$modeltype)

                 b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)
                 Walpha <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "precision")
                 gencook = b^2 * diag(w%*%solve(t(w)%*%Walpha%*%w)%*%t(w))
                 gencook <- c(gencook)
                 names(gencook) <- 1:length(gencook)
                 gencook
               },
               "LD" = {
                 if(is.null(model$x)){
                   stop("x component not found. fit the model again with argument x = TRUE")
                 }
                 if(is.null(model$w)){
                   stop("w component not found. fit the model again with argument w = TRUE")
                 }
                 x = model$x
                 w = model$w
                 mu = model$fitted.values
                 phi = model$fitted.precisions
                 if(is.null(model$y)){
                   y = model$residuals + model$fitted.values
                 } else{
                   y = model$y
                 }

                 lambda = lambda_r(y,mu,phi,model$modeltype)

                 a = mu*lambda - y

                 link_precision <- build_links_mpreg(model$link.precision)
                 eta_prec <- link_precision$linkfun(phi)
                 dphideta <- link_precision$mu.eta(eta_prec)

                 qsi0 <- switch(model$modeltype,
                                "NB" = {-1},
                                "PIG" = {-1/2}
                 )

                 B<- switch(model$modeltype,
                            "PIG" = {function(x) {-sqrt((-2*x))}},
                            "NB" = {function(x) {-log(-x)}}
                 )

                 dD<-switch(model$modeltype,
                            "PIG" = {function(x) {1/(2*x)}},
                            "NB" = {function(x) {log(x)+1-digamma(x)}}
                 )

                 kappa = kappa_r(y,mu,phi, model$modeltype)

                 b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)

                 covariates_matrix <- cbind(x, matrix(0,nrow=nrow(x), ncol = ncol(w)))
                 covariates_matrix <- rbind(covariates_matrix, cbind(matrix(0, nrow = nrow(x), ncol = ncol(x)), w))
                 coeff_beta <- model$coefficients$mean
                 coeff_alpha <- model$coefficients$precision
                 coeff <- c(coeff_beta, coeff_alpha)
                 new_coeff <- lapply(1:nrow(x), function(i){
                   W_Q_mean <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "mean")
                   W_Q_precision <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "precision")
                   new_beta <- coeff_beta + solve(t(x)%*% W_Q_mean %*% x)%*% x[i ,] * a[i]
                   new_alpha <- coeff_alpha + solve(t(w)%*% W_Q_precision %*% w)%*% w[i ,] * b[i]
                   c(new_beta, new_alpha)
                 })
                 lik_disp <- lapply(new_coeff, function(i){2*(loglik_mixpoisson(coeff, y, x, w, model$link.mean, model$link.precision, model$modeltype) - loglik_mixpoisson(i, y, x, w, model$link.mean, model$link.precision, model$modeltype))})
                 lik_disp <- c(unlist(lik_disp))
                 names(lik_disp) <- 1:length(lik_disp)
                 lik_disp
               },
               "QD" = {
                 if(is.null(model$x)){
                   stop("x component not found. fit the model again with argument x = TRUE")
                 }
                 if(is.null(model$w)){
                   stop("w component not found. fit the model again with argument w = TRUE")
                 }
                 x = model$x
                 w = model$w
                 mu = model$fitted.values
                 phi = model$fitted.precisions
                 if(is.null(model$y)){
                   y = model$residuals + model$fitted.values
                 } else{
                   y = model$y
                 }

                 lambda = lambda_r(y,mu,phi,model$modeltype)

                 a = mu*lambda - y

                 link_precision <- build_links_mpreg(model$link.precision)
                 eta_prec <- link_precision$linkfun(phi)
                 dphideta <- link_precision$mu.eta(eta_prec)

                 qsi0 <- switch(model$modeltype,
                                "NB" = {-1},
                                "PIG" = {-1/2}
                 )

                 B<- switch(model$modeltype,
                            "PIG" = {function(x) {-sqrt((-2*x))}},
                            "NB" = {function(x) {-log(-x)}}
                 )

                 dD<-switch(model$modeltype,
                            "PIG" = {function(x) {1/(2*x)}},
                            "NB" = {function(x) {log(x)+1-digamma(x)}}
                 )

                 kappa = kappa_r(y,mu,phi, model$modeltype)

                 b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)

                 covariates_matrix <- cbind(x, matrix(0,nrow=nrow(x), ncol = ncol(w)))
                 covariates_matrix <- rbind(covariates_matrix, cbind(matrix(0, nrow = nrow(x), ncol = ncol(x)), w))
                 coeff_beta <- model$coefficients$mean
                 coeff_alpha <- model$coefficients$precision
                 coeff <- c(coeff_beta, coeff_alpha)
                 new_coeff <- lapply(1:nrow(x), function(i){
                   W_Q_mean <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "mean")
                   W_Q_precision <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "precision")
                   new_beta <- coeff_beta + solve(t(x)%*% W_Q_mean %*% x)%*% x[i ,] * a[i]
                   new_alpha <- coeff_alpha + solve(t(w)%*% W_Q_precision %*% w)%*% w[i ,] * b[i]
                   c(new_beta, new_alpha)
                 })
                 Q_disp <- lapply(new_coeff, function(i){2*(Q_function_mixpoisson(coeff, mu, phi, y, x, w, model$link.mean, model$link.precision, model$modeltype) - Q_function_mixpoisson(i, mu, phi, y, x, w, model$link.mean, model$link.precision, model$modeltype))})
                 Q_disp <- c(unlist(Q_disp))
                 names(Q_disp) <- 1:length(Q_disp)
                 Q_disp
               }
               )

dist
}

#' @rdname influence.mixpoissonreg
#' @export

influence.mixpoissonreg <- function(model, do.coef = TRUE, ...){
influence_mpreg <- list()
hat_mean <- stats::hatvalues(model, parameters = "mean")
hat_precision <- stats::hatvalues(model, parameters = "precision")
influence_mpreg$hat.mean <- hat_mean
influence_mpreg$hat.precision <- hat_precision

if(do.coef){
  if(is.null(model$x)){
    stop("x component not found. fit the model again with argument x = TRUE")
  }
  if(is.null(model$w)){
    stop("w component not found. fit the model again with argument w = TRUE")
  }
  x = model$x
  w = model$w
  mu = model$fitted.values
  phi = model$fitted.precisions
  if(is.null(model$y)){
    y = model$residuals + model$fitted.values
  } else{
    y = model$y
  }

  lambda = lambda_r(y,mu,phi,model$modeltype)


  a = mu*lambda - y

  link_precision <- build_links_mpreg(model$link.precision)
  eta_prec <- link_precision$linkfun(phi)
  dphideta <- link_precision$mu.eta(eta_prec)

  qsi0 <- switch(model$modeltype,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  B<- switch(model$modeltype,
             "PIG" = {function(x) {-sqrt((-2*x))}},
             "NB" = {function(x) {-log(-x)}}
  )

  dD<-switch(model$modeltype,
             "PIG" = {function(x) {1/(2*x)}},
             "NB" = {function(x) {log(x)+1-digamma(x)}}
  )

  kappa = kappa_r(y,mu,phi, model$modeltype)

  b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)

  n_beta <- length(model$coefficients$mean)
  n_alpha <- length(model$coefficients$precision)

  coeff_beta <- model$coefficients$mean
  coeff_alpha <- model$coefficients$precision
  new_coeff <- lapply(1:nrow(x), function(i){
    W_Q_mean <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "mean")
    W_Q_precision <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "precision")
    new_beta <- coeff_beta + solve(t(x)%*% W_Q_mean %*% x)%*% x[i ,] * a[i]
    new_alpha <- coeff_alpha + solve(t(w)%*% W_Q_precision %*% w)%*% w[i ,] * b[i]
    c(new_beta, new_alpha)
  })

  influence_mpreg$coefficients.mean <- matrix(sapply(new_coeff, function(coeffs){coeffs[1:n_beta]}), nrow = nrow(x))
  colnames(influence_mpreg$coefficients.mean) <- names(model$coefficients$mean)
  rownames(influence_mpreg$coefficients.mean) <- 1:nrow(influence_mpreg$coefficients.mean)
  influence_mpreg$coefficients.precision <- matrix(sapply(new_coeff, function(coeffs){coeffs[(n_beta+1):(n_beta+n_alpha)]}), nrow = nrow(w))
  colnames(influence_mpreg$coefficients.precision) <- names(model$coefficients$precision)
  rownames(influence_mpreg$coefficients.precision) <- 1:nrow(influence_mpreg$coefficients.precision)
}

influence_mpreg$pear.res <- stats::residuals(model, type = "pearson")

influence_mpreg$score.res <- stats::residuals(model, type = "score")

influence_mpreg
}
