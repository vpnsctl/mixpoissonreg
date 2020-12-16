#' @import Rfast

#############################################################################################
#' @name local_influence.mixpoissonreg
#' @title Local Influence Diagnostics for Mixed Poisson Regression Models
#' @aliases local_influence.mixpoissonreg local_influence_plot.mixpoissonreg local_influence_ggplot.mixpoissonreg
#' @description These functions provides local influence diagnostic quantities. Currently the conformal normal and normal curvatures are available
#' under several perturbation schemes. The default is the conformal normal curvature since
#' it takes values on [0,1] and other nice properties (see Zhu and Lee, 2001 and Poon and Poon, 1999 for further details).
#' @param model a \code{mixpoissonreg} object.
#' @param perturbation a list or vector of perturbation schemes to be returned. The currently available schemes are
#' "case_weights", "hidden_variable", "mean_explanatory", "precision_explanatory", "simultaneous_explanatory". See Barreto-Souza and Simas (2015) for further details.
#' @param curvature the curvature to be returned, 'conformal' for the conformal normal curvature (see Zhu and Lee, 2001 and Poon and Poon, 1999) or
#' 'normal' (see Zhu and Lee, 2001 and Cook, 1986).
#' @param direction the 'max.eigen' returns the eigenvector associated to the largest eigenvalue of the perturbation matrix. The 'canonical' considers
#' the curvatures under the canonical directions, which is known as "total local curvature" (see Lesaffre and Verbeke, 1998). For conformal
#' normal curvatures both of them coincide. The default is 'canonical'.
#' @param parameters the parameter to which the local influence will be computed. The options are 'all', 'mean' and 'precision'.
#' This argument affects the 'case_weights' and 'hidden_variable' perturbation schemes. The default is 'all'.
#' @param mean.covariates a list of characters containing the mean-explanatory variables to be used in the 'mean-explanatory' and 'simultaneous-explanatory'
#' perturbation schemes. If NULL, the 'mean-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
#' mean-related covariates. The default is NULL.
#' @param precision.covariates a list of characters containing the precision-explanatory variables to be used in the 'precision-explanatory'
#' and 'simultaneous-explanatory'
#' perturbation schemes. If NULL, the 'precision-explanatory' and 'simultaneous-explanatory' perturbation schemes will be computed by perturbing all
#' precision-related covariates. The default is NULL.
#' @return a list containing the resulting perturbation schemes as elements.
#' @details Preencher
#'
#' @references Zhu and Lee; Barreto-Souza and Simas; Cook;  etc..

#' @rdname local_influence.mixpoissonreg
#' @export
local_influence.mixpoissonreg <- function(model, perturbation = c("case_weights", "hidden_variable",
                                                                  "mean_explanatory", "precision_explanatory",
                                                                  "simultaneous_explanatory"), curvature = c("conformal", "normal"),
                                          direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                          mean.covariates = NULL, precision.covariates = NULL){
loc_infl <- list()

if(length(parameters)>1){
  parameters = parameters[1]
}

if(length(direction)>1){
  direction = direction[1]
}

if(length(curvature)>1){
  curvature = curvature[1]
}

modeltype <- model$modeltype

if(is.null(model$x)){
  stop("x component not found. fit the model again with argument x = TRUE")
}
if(is.null(model$w)){
  stop("w component not found. fit the model again with argument w = TRUE")
}

n = model$nobs
x = model$x
w = model$w

beta = model$coefficients$mean
alpha = model$coefficients$precision

n_beta = length(model$coefficients$mean)
n_alpha = length(model$coefficients$precision)

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

A = diag(c(a))
B1 = diag(c(b))



for(pert in perturbation){
  switch(pert,
         "case_weights" = {
           Bmatrix = (1 - (parameters == "precision")) * (A%*%x) %*% solve(t(x)%*%Wbeta%*%x) %*% t(A%*%x) +
             (1 - (parameters == "mean")) * (B1%*%w)%*%solve(t(w)%*%Walpha%*%w)%*% t(B1%*%w)
         },
         "hidden_variable" = {
           P1 = diag(-c(mu*lambda))
           P2 = diag(c(phi*(lambda*qsi0-kappa)) )

           Bmatrix = (1 - (parameters == "precision")) * (P1%*%x) %*% solve(t(x)%*%Wbeta%*%x) %*% t(P1%*%x) +
             (1 - (parameters == "mean")) * (P2%*%w)%*%solve(t(w)%*%Walpha%*%w)%*% t(P2%*%w)

         },
         "mean_explanatory" = {
           S_x = sqrt(colVars(x))
           if(!is.null(mean.covariates)){
             expl_mean <- names(model$coefficients$mean) %in% mean.covariates
             S_x = S_x * expl_mean
           }
           S_x = diag(S_x)

           ones_nbeta = rep(1,n_beta)
           ones_n = rep(1,n)
           T1 = A%*%ones_n%*%ones_nbeta%*%S_x+Wbeta%*%x%*%beta%*%ones_nbeta%*%S_x

           Bmatrix =  T1 %*% solve(t(x)%*%Wbeta%*%x) %*% t(T1)
         },
         "precision_explanatory" = {
           S_w = sqrt(colVars(w))
           if(!is.null(precision.covariates)){
             expl_precision <- names(model$coefficients$precision) %in% precision.covariates
             S_w = S_w * expl_precision
           }

           S_w = diag(S_w)

           ones_nalpha = rep(1,n_alpha)
           ones_n = rep(1,n)
           T2 = B1%*%ones_n%*%ones_nalpha%*%S_w+Walpha%*%w%*%alpha%*%ones_nalpha%*%S_w

           Bmatrix =  T2 %*% solve(t(w)%*%Walpha%*%w) %*% t(T2)

         },
         "simultaneous_explanatory" = {
           S_x = sqrt(colVars(x))
           if(!is.null(mean.covariates)){
             expl_mean <- names(model$coefficients$mean) %in% mean.covariates
             S_x = S_x * expl_mean
           }

           S_w = sqrt(colVars(w))
           if(!is.null(precision.covariates)){
             expl_precision <- names(model$coefficients$precision) %in% precision.covariates
             S_w = S_w * expl_precision
           }

           S_x = diag(S_x)
           S_w = diag(S_w)

           ones_nbeta = rep(1,n_beta)
           ones_nalpha = rep(1,n_alpha)
           ones_n = rep(1,n)

           T1 = A%*%ones_n%*%ones_nbeta%*%S_x+Wbeta%*%x%*%beta%*%ones_nbeta%*%S_x
           T2 = B1%*%ones_n%*%ones_nalpha%*%S_w+Walpha%*%w%*%alpha%*%ones_nalpha%*%S_w

           Bmatrix = T1 %*% solve(t(x)%*%Wbeta%*%x) %*% t(T1) + T2 %*% solve(t(w)%*%Walpha%*%w) %*% t(T2)
         }
         )

  Bmatrix <- switch(curvature,
                    "conformal" = {Bmatrix/sum(diag(Bmatrix))},
                    "normal" = Bmatrix)


  loc_infl[[pert]] = switch(direction,
                                 "max.eigen" = {eigen(Bmatrix, symmetric = TRUE)$vec[,1]},
                                 "canonical" = {diag(Bmatrix)})
  names(loc_infl[[pert]]) = 1:n
  benchmark = switch(curvature,
                     "conformal" = {
                       ifelse(direction == "canonical", 1/n + 2*sd(loc_infl[[pert]]), NA) #Zhu and Lee (2001)
                     },
                     "normal" = {
                       ifelse(direction == "canonical", 2 * mean(loc_infl[[pert]]), NA) # Verbeke and Molenberghs (2000, sect. 11.3)
                     }
  )

  attr(loc_infl[[pert]], "benchmark") = benchmark
}


loc_infl
}
