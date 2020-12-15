#############################################################################################
#' @name global.influence.mixpoissonreg
#' @title Global Influence Diagnostics for Mixed Poisson Regression Models
#' @aliases cooks.distance.mixpoissonreg hatvalues.mixpoissonreg
#' @description These functions provides global influence diagnostic quantities such as Cook's distance, hat values, generalized Cook's distance (through argument
#' on \code{cooks.distance.mixpoissonreg} function), likelihood displacement (through argument
#' on \code{cooks.distance.mixpoissonreg} function) and Q-displacement (through argument
#' on \code{cooks.distance.mixpoissonreg} function).
#' @param model a \code{mixpoissonreg} object.
#' @param hat hat values \code{H[i,i]}. The default is obtained through the second-derivative of the Q-function in the spirit of Zhu et al. ();
#' @param type the type of Cook's distance to be used. The options are 'CD', the standard Cook's distance; 'GCD', the generalized Cook's distance with respect to all parameters;
#' 'GCDmean', the generalized Cook's distance with respect to the mean parameters; 'GCDprecision', the generalized Cook's distance with respect to the precision parameters;
#' 'LD', the likelihood displacement (also known as likelihood distance); 'QD', the Q-displacement. See 'details' for more informations. For 'GCD', 'GCDmean', 'LD' and 'QD', the model must be fitted with 'x' set to TRUE,
#' and for 'GCD', 'GCDprecision', 'LD' and 'QD', the model must be fitted with 'w' set to TRUE. For 'CD', if 'hat' is set to 'mean', the model must be fitted with 'x' set to TRUE, whereas
#' if 'hat' is set to 'precision', the model must be fitted with 'w' set to TRUE.
#' @param parameters the parameter to which the hat values will be computed. The options are 'mean' and 'precision'. The default is 'mean'. For hatvalues with respect to the mean
#' the model must be fitted with 'x' set to TRUE, and for hatvalues with respect to the precision the model must be fitted with 'w' set to TRUE.
#' @param do.coef logical indicating if the the approximation to the change of coefficients values after case removal are desired. The model must be fitted with x = TRUE. See details for further explanations.
#' @details Preencher
#'
#' @references Zhu et al.; Barreto-Souza and Simas; Cook and Preigibon; etc..

#' @rdname global.influence.mixpoissonreg
#' @export

hatvalues.mixpoissonreg <- function(model, parameters = c("mean", "precision")){
  if(length(parameters)>1){
    parameters = parameters[1]
  }
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

#' @rdname global.influence.mixpoissonreg
#' @export

cooks.distance.mixpoissonreg <- function(model, type = c("CD", "GCD", "GCDmean", "GCDprecision", "LD", "QD"), hat = c("mean", "precision")){
if(length(type)>1){
  type = type[1]
}
if(!(type%in%c("CD", "GCD", "GCDmean", "GCDprecision", "LD", "QD"))){
  stop("type must be one of 'CD', 'GCD', 'GCDmean'', 'GCDprecision', 'LD', 'QD'")
}

if(length(hat)>1){
  hat = hat[1]
}

if(!(hat%in%c("mean","precision"))){
  stop("hat must be one of 'mean' or 'precision'")
}

h = switch(hat,
           "mean" = {
             if(is.null(model$x)){
               stop("x component not found. fit the model again with argument x = TRUE")
             }
             hatvalues(model, parameters = hat)
             },
           "precision" = {
             if(is.null(model$w)){
               stop("w component not found. fit the model again with argument w = TRUE")
             }
             hatvalues(model, parameters = hat)
           })
dist <- switch(type,
               "CD" = {
                 p <- length(model$coefficients$mean)
                 pearson_residuals <- residuals(model, type = "Pearson")
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
                 coeff <- c(model$coefficients$mean, model$coefficients$precision)
                 new_coeff <- lapply(1:nrow(x), function(i){
                   W_Q <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "all")
                   coeff + solve(t(covariates_matrix)%*% W_Q %*% covariates_matrix)%*%covariates_matrix[i ,] * c(rep(a[i], n_beta), rep(0, n_alpha)) +
                     solve(t(covariates_matrix)%*% W_Q %*% covariates_matrix)%*%covariates_matrix[nrow(x)+i ,] * c(rep(0, n_beta), rep(b[i], n_alpha))
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
                 coeff <- c(model$coefficients$mean, model$coefficients$precision)
                 new_coeff <- lapply(1:nrow(x), function(i){
                   W_Q <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "all")
                   coeff + solve(t(covariates_matrix)%*% W_Q %*% covariates_matrix)%*%covariates_matrix[i ,] * c(rep(a[i], n_beta), rep(0, n_alpha)) +
                     solve(t(covariates_matrix)%*% W_Q %*% covariates_matrix)%*%covariates_matrix[nrow(x)+i ,] * c(rep(0, n_beta), rep(b[i], n_alpha))
                 })
                 Q_disp <- lapply(new_coeff, function(i){2*(Q_function_mixpoisson(coeff, mu, phi, y, x, w, model$link.mean, model$link.precision, model$modeltype) - Q_function_mixpoisson(i, mu, phi, y, x, w, model$link.mean, model$link.precision, model$modeltype))})
                 Q_disp <- c(unlist(Q_disp))
                 names(Q_disp) <- 1:length(Q_disp)
                 Q_disp
               }
               )

dist
}

#' @rdname global.influence.mixpoissonreg
#' @export

influence.mixpoissonreg <- function(model, do.coef = TRUE){
influence_mpreg <- list()
hat_mean <- hatvalues(model, parameters = "mean")
hat_precision <- hatvalues(model, parameters = "precision")
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

  covariates_matrix <- cbind(x, matrix(0,nrow=nrow(x), ncol = ncol(w)))
  covariates_matrix <- rbind(covariates_matrix, cbind(matrix(0, nrow = nrow(x), ncol = ncol(x)), w))
  coeff <- c(model$coefficients$mean, model$coefficients$precision)
  new_coeff <- lapply(1:nrow(x), function(i){
    W_Q <- obs_fisher_weight_matrix_mixpoisson(model, parameters = "all")
    coeff + solve(t(covariates_matrix)%*% W_Q %*% covariates_matrix)%*%covariates_matrix[i ,] * c(rep(a[i], n_beta), rep(0, n_alpha)) +
      solve(t(covariates_matrix)%*% W_Q %*% covariates_matrix)%*%covariates_matrix[nrow(x)+i ,] * c(rep(0, n_beta), rep(b[i], n_alpha))
  })
  influence_mpreg$coefficients.mean <- matrix(sapply(new_coeff, function(coeffs){coeffs[1:n_beta]}), nrow = nrow(x))
  colnames(influence_mpreg$coefficients.mean) <- names(model$coefficients$mean)
  rownames(influence_mpreg$coefficients.mean) <- 1:nrow(influence_mpreg$coefficients.mean)
  influence_mpreg$coefficients.precision <- matrix(sapply(new_coeff, function(coeffs){coeffs[(n_beta+1):(n_beta+n_alpha)]}), nrow = nrow(w))
  colnames(influence_mpreg$coefficients.precision) <- names(model$coefficients$precision)
  rownames(influence_mpreg$coefficients.precision) <- 1:nrow(influence_mpreg$coefficients.precision)


}
influence_mpreg
}