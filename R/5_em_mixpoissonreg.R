#' @import gamlss.dist

#############################################################################################
#' @title EM_mixpoisson
#' @description Function to run the Expectation-Maximization algorithm for mixed Poisson regression models.
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
#' @param em_controls only used with the 'EM' method. A list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm;
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm. \code{em_tolgrad} that defines the tolerance value
#' of the maximum-norm of the the gradient of the Q-function, the default is set to 10^(-2).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function.
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return a list containing
#' \itemize{
#'     \item niter - number of EM iterations;
#'     \item coefficients - a list containing estimated vectors 'mean' (mean-related coefficients) and 'precision' (precision-related coefficients);
#'     \item fitted.values - the estimated means.
#' }
EM_mixpoisson <- function(beta, alpha, y, x, w,
                          link.mean, link.precision, model, em_controls, optim_method, optim_controls) {
  nbeta <- length(beta)

  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  theta <- c(beta, alpha)
  count <- 0

  maxit <- em_controls$maxit
  epsilon <- em_controls$em_tol
  epsilon_grad <- em_controls$em_tolgrad

  repeat{
    theta_old <- theta
    beta <- theta[1:nbeta]
    alpha <- theta[-c(1:nbeta)]
    mu <- link_mean$linkinv(x %*% beta)
    phi <- link_precision$linkinv(w %*% alpha)
    ### E step ------------------------------
    mu_old <- mu
    phi_old <- phi

    ### M step ------------------------------
    M <- tryCatch(stats::optim(
      par = theta,
      fn = Q_function_mixpoisson,
      gr = gradient_Q_function_mixpoisson,
      muold = mu_old,
      phiold = phi_old,
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
      M <- tryCatch(stats::optim(
        par = theta,
        fn = Q_function_mixpoisson,
        gr = NULL,
        muold = mu_old,
        phiold = phi_old,
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
      M <- tryCatch(stats::optim(
        par = theta,
        fn = Q_function_mixpoisson,
        gr = NULL,
        muold = mu_old,
        phiold = phi_old,
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
      warning("The EM algorithm did not converge")
      break
    }

    theta <- M$par
    # Compute Q -----------------------------
    Q_old <- Q_function_mixpoisson(theta_old, mu_old, phi_old, y, x, w, link.mean, link.precision,model)
    Q <- M$value
    ### Convergence criterion ---------------
    term1 <- sqrt(sum((theta - theta_old)^2))
    term2 <- abs(Q - Q_old)
    term_grad <- max(abs(gradient_Q_function_mixpoisson(theta_old, mu_old, phi_old, y, x, w, link.mean, link.precision,model)))

    ### -------------------------------------
    count <- count + 1
    if (max(term1, term2) < epsilon & term_grad < epsilon_grad) {
      break
    }
    if (count >= maxit) {
      warning("The EM algorithm did not converge")
      break
    }
    ### -------------------------------------
  }

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
  out$niter <- count
  return(out)
}


#############################################################################################
#' @title Q_function_mixpoisson
#' @description Q-function for mixed Poisson regression models. This function is required in the EM-algorithm.
#' @param theta vector of parameters (all coefficients).
#' @param muold previous value of the mean parameter (mu).
#' @param phiold previous value of the precision parameter (phi).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return scalar representing the value of the Q-function at 'theta'.
Q_function_mixpoisson <- function(theta, muold, phiold, y, x, w,
                                  link.mean, link.precision,
                                  model) {
  nbeta <- ncol(x)
  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  B<- switch(model,
             "PIG" = {function(x) {-sqrt((-2*x))}},
             "NB" = {function(x) {-log(-x)}}
  )

  D<-switch(model,
            "PIG" = {function(x) {log(x)/2}},
            "NB" = {function(x) {x*log(x)-lgamma(x)}}
  )

  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  sum(y*log(mu)-mu*lambda_r(y,muold,phiold,model)+D(phi)+phi*
          (lambda_r(y,muold,phiold, model)*qsi0-
             B(qsi0)+kappa_r(y,muold,phiold, model)))+
      sum((y-1)*kappa_r(y,muold,phiold,model))
}

#############################################################################################
#' @title gradient_Q_function_mixpoisson
#' @description Function to calculate the gradient of the Q-function of mixed Poisson regression models,
#' which is required for optimization via \code{optim}.
#' @param theta vector of parameters (all coefficients).
#' @param phiold previous value of the precision parameter (phi).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return vector containing the gradient of the Q function at theta.
gradient_Q_function_mixpoisson <- function(theta, muold, phiold, y, x, w,
                                           link.mean, link.precision, model) {
  nbeta <- ncol(x)
  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]
  nalpha <- length(alpha)

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  B<- switch(model,
             "PIG" = {function(x) {-sqrt((-2*x))}},
             "NB" = {function(x) {-log(-x)}}
  )

  dD<-switch(model,
            "PIG" = {function(x) {1/(2*x)}},
            "NB" = {function(x) {log(x)+1-digamma(x)}}
  )

  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  dmudeta <- link_mean$mu.eta(x %*% beta)
  dphideta <- link_precision$mu.eta(w %*% alpha)

  grad1<-t(x)%*%( ( (y-mu*lambda_r(y,muold,phiold,model))/mu ) * dmudeta )
  grad2<-t(w)%*%(dphideta*(qsi0*lambda_r(y,muold,phiold,model)-
                        B(qsi0)+kappa_r(y,muold,phiold,model)+dD(phi)))
  c(grad1,grad2)
}



#############################################################################################
#' @title obs_fisher_information_mixpoisson
#' @description Function to compute the Fisher's information matrix for mixed Poisson regression models.
#' @param theta vector of parameters (all coefficients).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return Fisher's information matrix.

obs_fisher_information_mixpoisson <- function(theta, y, x, w,
                                              link.mean, link.precision, model) {
  d2.Q <- D2Q_Obs_Fisher_mixpoisson(theta, y, x, w, link.mean, link.precision, model)
  dd.Q <- DQ2_Obs_Fisher_mixpoisson(theta, y, x, w, link.mean, link.precision, model)
  inf_fisher <- d2.Q - dd.Q # Fisher Information Matrix.
  inf_fisher
}


#############################################################################################
#' @title D2Q_Obs_Fisher_mixpoisson
#' @description Auxiliary function to compute the observed Fisher information matrix for mixed Poisson regression models.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return Hessian of the Q-function.

D2Q_Obs_Fisher_mixpoisson <- function(theta, y, x, w, link.mean, link.precision, model) {

  nbeta <- ncol(x)
  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]
  nalpha <- length(alpha)

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  B<- switch(model,
             "PIG" = {function(x) {-sqrt((-2*x))}},
             "NB" = {function(x) {-log(-x)}}
  )

  dD<-switch(model,
             "PIG" = {function(x) {1/(2*x)}},
               "NB" = {function(x) {log(x)+1-digamma(x)}}
  )

  ddD<-switch(model,
             "PIG" = {function(x){-1/(2*x^2)}},
               "NB" = {function(x){1/x-trigamma(x)}}
  )

  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  dmudeta <- link_mean$mu.eta(x %*% beta)
  dphideta <- link_precision$mu.eta(w %*% alpha)

  d2mu <- d2mudeta2(link.mean, mu)
  d2phi <- d2phideta2(link.precision, phi)

  lambda = lambda_r(y,mu,phi,model)
  kappa = kappa_r(y,mu,phi, model)

  G1 = diag(c( d2mu*(mu*lambda - y)/mu + dmudeta^2*y/(mu^2)))
  G2 = diag(c(d2phi*(-qsi0*lambda+B(qsi0)-kappa-dD(phi))-ddD(phi)*
                        dphideta^2))

    D2QBeta = t(x)%*%G1%*%x
    D2QBeta = cbind(D2QBeta,matrix(0,nbeta,nalpha))
    D2QAlpha = t(w)%*%G2%*%w
    D2QAlpha = cbind(matrix(0,nalpha,nbeta),D2QAlpha)
    D2Q = rbind(D2QBeta,D2QAlpha)
    D2Q
}

#############################################################################################
#' @title DQ2_Obs_Fisher_mixpoisson
#' @description Auxiliary function to compute the observed Fisher information matrix for mixed Poisson regression models.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param y response vector with y_i>=0 and integer.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param w matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log" and "inverse.sqrt".
#' @param model the mixed Poisson model, "NB" or "PIG".
#' @return matrix given by the conditional expectation of the gradient of the Q-function and its tranpose.

DQ2_Obs_Fisher_mixpoisson <- function(theta, y, x, w, link.mean, link.precision, model) {
  nbeta <- ncol(x)
  beta <- theta[1:nbeta]
  alpha <- theta[-c(1:nbeta)]
  nalpha <- length(alpha)

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  B<- switch(model,
             "PIG" = {function(x) {-sqrt((-2*x))}},
             "NB" = {function(x) {-log(-x)}}
  )

  dD<-switch(model,
             "PIG" = {function(x) {1/(2*x)}},
             "NB" = {function(x) {log(x)+1-digamma(x)}}
  )

  ddD<-switch(model,
              "PIG" = {function(x){-1/(2*x^2)}},
              "NB" = {function(x){1/x-trigamma(x)}}
  )

  link_mean <- build_links_mpreg(link.mean)
  link_precision <- build_links_mpreg(link.precision)

  mu <- link_mean$linkinv(x %*% beta)
  phi <- link_precision$linkinv(w %*% alpha)

  dmudeta <- link_mean$mu.eta(x %*% beta)
  dphideta <- link_precision$mu.eta(w %*% alpha)

  d2mu <- d2mudeta2(link.mean, mu)
  d2phi <- d2phideta2(link.precision, phi)

  lambda = lambda_r(y,mu,phi,model)
  kappa = kappa_r(y,mu,phi, model)

  b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)


  gama = switch(model,
                "PIG" = {(y+1)*(y+2)*dPIG(y+2,mu,1/phi)/(mu^2*dPIG(y,mu,1/phi))},
                "NB" = {(phi+y+1)*(phi+y)/(mu+phi)^2}
  )

  rho = switch(model,
               "PIG" = {-1/2},
               "NB" = {(y+phi)/(mu+phi)*(digamma(y+phi+1)-log(mu+phi))}
  )

  nu = switch(model,
              "PIG" = {(y==0)*(1/(4*phi^2))*(phi*(2*mu+phi)+3*
                                  (1+sqrt(phi*(2*mu+phi))))+(y==1)*((1+sqrt(phi*(2*mu+phi))))*sqrt(2*mu+phi)/(4*phi^(3/2))+(y>=2)*
    mu^2*dPIG((y>=2)*(y-2),mu,1/phi)/
    (4*(y+(y==0))*(y-1+2*(y<=1))*dPIG(y,mu,1/phi))},
    "NB" = { trigamma(y+phi)+digamma(y+phi)^2-2*digamma(y+phi)*
        log(mu+phi)+log(mu+phi)^2}
  )

  G3_diag = c( dmudeta^2 * (y^2-2*y*mu*lambda+mu^2*gama)/mu^2 )
  G3 = (dmudeta * (y-mu*lambda)/mu)%*%t(dmudeta * (y-mu*lambda)/mu)
  diag(G3) = G3_diag

  G4_diag = c( (dmudeta*dphideta/mu)*(-y*b/dphideta - mu*(qsi0*gama+rho+lambda*
                                             (dD(phi)-B(qsi0)))))
  G4 = (dmudeta * (y-mu*lambda)/mu)%*%t(-b)
  diag(G4) = G4_diag

  G5_diag = c(dphideta^2*((dD(phi)-B(qsi0))^2+2*(dD(phi)-B(qsi0))*
                          (qsi0*lambda+kappa)+qsi0^2*gama+2*qsi0*rho+nu))
  G5 = b%*%t(b)
  diag(G5) = G5_diag

    DQ2Beta = t(x)%*%G3%*%x
    DQBetaAlfa = t(x)%*%G4%*%w
    DQ2Alfa = t(w)%*%G5%*%w

    DQ21 = cbind(DQ2Beta,DQBetaAlfa)
    DQ22 = cbind(t(DQBetaAlfa),DQ2Alfa)
    DQ2 = rbind(DQ21,DQ22)
    DQ2
}


#############################################################################################
#' @title obs_fisher_weight_matrix_mixpoisson
#' @description Function to compute the weight matrix based on Fisher's information matrix for mixed Poisson regression models.
#' @param object a \code{mixpoissonreg} object.
#' @param parameters the parameter to which the matrix will be returned. The options are 'all', 'mean' and 'precision'. The default is 'all'.
#' @return Fisher's information matrix.

obs_fisher_weight_matrix_mixpoisson <- function(object, parameters = c("all", "mean", "precision")) {
  if(length(parameters)>1){
    parameters = parameters[1]
  }
  link_mean <- build_links_mpreg(object$link.mean)
  link_precision <- build_links_mpreg(object$link.precision)
  mu <- object$fitted.values
  phi <- object$fitted.precisions
  eta_mean <- link_mean$linkfun(mu)
  eta_prec <- link_precision$linkfun(phi)
  dmudeta <- link_mean$mu.eta(eta_mean)
  dphideta <- link_precision$mu.eta(eta_prec)

  d2mu <- d2mudeta2(object$link.mean, mu)
  d2phi <- d2phideta2(object$link.precision, phi)

  if(is.null(object$y)){
    y = object$residuals + object$fitted.values
  } else{
    y = object$y
  }

  model <- object$modeltype

  lambda = lambda_r(y,mu,phi,model)
  kappa = kappa_r(y,mu,phi, model)

  qsi0 <- switch(model,
                 "NB" = {-1},
                 "PIG" = {-1/2}
  )

  B<- switch(model,
             "PIG" = {function(x) {-sqrt((-2*x))}},
             "NB" = {function(x) {-log(-x)}}
  )

  dD<-switch(model,
             "PIG" = {function(x) {1/(2*x)}},
             "NB" = {function(x) {log(x)+1-digamma(x)}}
  )

  ddD<-switch(model,
              "PIG" = {function(x){-1/(2*x^2)}},
              "NB" = {function(x){1/x-trigamma(x)}}
  )

  b = dphideta*(B(qsi0)-kappa - dD(phi)-qsi0*lambda)


  gama = switch(model,
                "PIG" = {(y+1)*(y+2)*dPIG(y+2,mu,1/phi)/(mu^2*dPIG(y,mu,1/phi))},
                "NB" = {(phi+y+1)*(phi+y)/(mu+phi)^2}
  )

  rho = switch(model,
               "PIG" = {-1/2},
               "NB" = {(y+phi)/(mu+phi)*(digamma(y+phi+1)-log(mu+phi))}
  )

  nu = switch(model,
              "PIG" = {(y==0)*(1/(4*phi^2))*(phi*(2*mu+phi)+3*
                                               (1+sqrt(phi*(2*mu+phi))))+(y==1)*((1+sqrt(phi*(2*mu+phi))))*sqrt(2*mu+phi)/(4*phi^(3/2))+(y>=2)*
                  mu^2*dPIG((y>=2)*(y-2),mu,1/phi)/
                  (4*(y+(y==0))*(y-1+2*(y<=1))*dPIG(y,mu,1/phi))},
              "NB" = { trigamma(y+phi)+digamma(y+phi)^2-2*digamma(y+phi)*
                  log(mu+phi)+log(mu+phi)^2}
  )

  n <- length(y)

  if(parameters %in% c("all", "mean")){
    G1 = diag(c( d2mu*(mu*lambda - y)/mu + dmudeta^2*y/(mu^2)))
    D2QBeta = G1
  }
  if(parameters == "all"){
    D2QBeta = cbind(D2QBeta,matrix(0,n,n))
  }
  if(parameters %in% c("all", "precision")){
    G2 = diag(c(d2phi*(-qsi0*lambda+B(qsi0)-kappa-dD(phi))-ddD(phi)*
                  dphideta^2))
    D2QAlpha = G2
  }

  if(parameters == "all"){
    D2QAlpha = cbind(matrix(0,n,n),D2QAlpha)
    D2Q = rbind(D2QBeta,D2QAlpha)
  }


  if(parameters == "all"){
    return(D2Q)
  } else if(parameters == "mean"){
    return(D2QBeta)
  } else if(parameters == "precision"){
    return(D2QAlpha)
  } else{
    stop("the parameters argument must be 'all', 'mean' or 'precision'.")
  }
}
