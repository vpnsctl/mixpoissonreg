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
#' @param which a list or vector indicating which plots should be displayed. 	If a subset of the plots is required, specify a subset of the numbers 1:5, see caption below (and the 'Details') for the different kinds.
#' @param caption captions to appear above the plots; character vector or list of valid graphics annotations. Can be set to "" or NA to suppress all captions.
#' @param sub.caption	common title-above the figures if there are more than one. If NULL, as by default, a possible abbreviated version of deparse(x$call) is used.
#' @param kind character indicating if the plots should be done with base R or 'ggplot'.
#' @param detect.influential logical. Indicates whether the benchmark should be used to detect influential observations and identify them on the plot. If there is no benchmark available,
#' the top 'n.influential' observations will be identified in the plot by their indexes.
#' @param n.influential interger. The maximum number of influential observations to be identified on the plot.
#' @param draw.benchmark logical. Indicates whether a horizontal line identifying the benchmark should be drawn.
#' @param lty.benchmark the line type of the benchmark if drawn.
#' @param type_plot what type of plot should be drawn. The default is 'h'.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param main character; title to be placed at each plot additionally (and above) all captions.
#' @param cex.id magnification of point labels.
#' @param cex.caption	controls the size of caption.
#' @param cex.oma.main controls the size of the sub.caption only if that is above the figures when there is more than one.
#' @param include.modeltype logical. Indicates whether the model type ('NB' or 'PIG') should be displayed on the captions.
#' @param ... other graphical arguments to be passed.
#' @return a list containing the resulting perturbation schemes as elements. Each returned element has an attribute 'benchmark', which for
#' the conformal normal curvature, it is computed following Zhu and Lee (2001), and for normal curvature it is computed following Verbeke and ...
#' If the 'direction' is 'max.eigen' the 'benchmark' attribute is NA.
#'
#' The 'mean_explanatory', 'precision_explanatory' and 'simultaneous_explanatory' elements of the list contain an attribute 'covariates' indicating
#' which covariates were used in the perturbation schemes.
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

  if(pert == "mean_explanatory"){
    attr(loc_infl[[pert]], "covariates") = mean.covariates
    if(is.null(attr(loc_infl[[pert]], "covariates"))){
      attr(loc_infl[[pert]], "covariates") = "all"
    }
  }
  if(pert == "precision_explanatory"){
    attr(loc_infl[[pert]], "covariates") = precision.covariates
    if(is.null(attr(loc_infl[[pert]], "covariates"))){
      attr(loc_infl[[pert]], "covariates") = "all"
    }
  }
  if(pert == "simultaneous_explanatory"){
    attr(loc_infl[[pert]], "covariates") = list("mean" = mean.covariates, "precision" = precision.covariates)
    if(is.null(attr(loc_infl[[pert]], "covariates")$mean)){
      attr(loc_infl[[pert]], "covariates")$mean = "all"
    }
    if(is.null(attr(loc_infl[[pert]], "covariates")$precision)){
      attr(loc_infl[[pert]], "covariates")$precision = "all"
    }
  }
}


loc_infl
}

#' @rdname local_influence.mixpoissonreg
#' @export

local_influence_plot.mixpoissonreg <- function(model, which = c(1,2,3,5),
                                               caption = list("Case Weights Perturbation",
                                                               "Hidden Variable Perturbation",
                                                               "Mean Explanatory Perturbation",
                                                               "Precision Explanatory Perturbation",
                                                               "Simultaneous Explanatory Perturbation"),
                                               sub.caption = NULL,
                                               kind = c("base", "ggplot"),
                                               detect.influential = TRUE, n.influential = 5,
                                               draw.benchmark = FALSE, lty.benchmark = 2,
                                               type_plot = "h",
                                               curvature = c("conformal", "normal"),
                                               direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                               mean.covariates = NULL, precision.covariates = NULL, main = "",
                                               ask = prod(par("mfcol")) <
                                                 length(which) && dev.interactive(),
                                               cex.id = 0.75,
                                               cex.oma.main = 1.25,
                                               cex.caption = 1,
                                               include.modeltype = TRUE,
                                               ...){
if(length(kind)>1){
  kind = kind[1]
}
  if(length(direction)>1){
    direction = direction[1]
  }

  if(length(curvature)>1){
    curvature = curvature[1]
  }

  getCaption <- function(k) if (length(caption) < k)
    NA_character_
  else {
    if(include.modeltype){
      as.graphicsAnnot(paste0(caption[[k]], " - ", model$modeltype, " Regression"))
    } else {
      as.graphicsAnnot(caption[[k]])
    }
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

  one.fig <- prod(par("mfcol")) == 1

loc_infl <- local_influence(model, curvature = curvature,
                            direction = direction, parameters = parameters,
                            mean.covariates = mean.covariates,
                            precision.covariates = precision.covariates)

if (ask) {
  oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
}

if(kind == "base"){
pert <- c("case_weights", "hidden_variable",
          "mean_explanatory", "precision_explanatory",
          "simultaneous_explanatory")
for(i in 1:length(pert)){
  if(i %in% which){

    xlab <- paste0("Index\n ", sub.caption)

    graphics::plot(loc_infl[[pert[i]]], type = type_plot, main = main, ylab = ylab_infl, ...)
    graphics::mtext(getCaption(i), side = 3, cex = cex.caption)

    if (one.fig)
      title(sub = sub.caption, ...)

    if(detect.influential){
      bm <- attr(loc_infl[[pert[i]]], "benchmark")
      infl_points <- as.vector(nth(abs(loc_infl[[pert[i]]]), k = n.influential,
                                   num.of.nths = n.influential,
                                   index.return = TRUE, descending = TRUE))

      if(!is.na(bm)){
        infl_points <- (loc_infl[[pert[i]]][infl_points] > bm)
        idx_x <- as.integer(names(which(infl_points == TRUE)))
        idx_y <- loc_infl[[pert[i]]][idx_x]
        text(idx_x, idx_y, labels = idx_x, cex = cex.id, xpd = TRUE, pos = 3, offset = 0.1)
      } else{
        idx_x_pos <- infl_points[which(loc_infl[[pert[i]]][infl_points] >= 0)]
        idx_x_neg <- setdiff(infl_points, idx_x_pos)
        idx_y_pos <- loc_infl[[pert[i]]][idx_x_pos]
        idx_y_neg <- loc_infl[[pert[i]]][idx_x_neg]
        if(length(idx_x_pos)>0){
          text(idx_x_pos, idx_y_pos, labels = idx_x_pos, cex = cex.id, xpd = TRUE, pos = 3, offset = 0.1)
        }
        if(length(idx_x_neg)>0){
          text(idx_x_neg, idx_y_neg, labels = idx_x_neg, cex = cex.id, xpd = TRUE, pos = 1, offset = 0.1)
        }
      }


    }
    if(draw.benchmark & !is.na(bm)){
      bm <- attr(loc_infl[[pert[i]]], "benchmark")
      abline(a = bm, b = 0, lty = lty.benchmark)
    }

    dev.flush()
  }
}



} else if(kind == "ggplot"){

  }

if (!one.fig && par("oma")[3L] >= 1)
  mtext(sub.caption, outer = TRUE, cex = 1.25)

invisible()
}


#############################################################################################
#' @name local_influence
#' @aliases local_influence local_influence_plot local_influence_ggplot
#' @title Local Influence Diagnostics
#' @usage local_influence(model, ...)
#' @param model an object for which the local influence is desired
#' @param ... further arguments passed to or from other methods.
#' @details Preencher
#' @references Cook, Zhu and Lee, Poon and Poon, etc..
#' @seealso Preencher
#'
#' @rdname local_influence
#' @export

local_influence <- function(model, ...){
UseMethod("local_influence", model)
}

#' @rdname local_influence
#' @export

local_influence_plot <- function(model, ...){
  UseMethod("local_influence_plot", model)
}

#' @rdname local_influence
#' @export

local_influence_ggplot <- function(model, ...){
  UseMethod("local_influence_ggplot", model)
}
