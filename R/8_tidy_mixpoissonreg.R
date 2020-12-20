#' @import generics
#' @import ggplot2
#' @import tibble
#' @import magrittr


augment.mixpoissonreg <- function(x, data = model.frame(x), newdata = NULL, type.predict = c("response", "link", "precision", "variance"), 
                                  type.residuals = c("pearson", "score"), se_fit = FALSE, conf_int = TRUE, pred_int = FALSE, ...) {
  type.predict <- rlang::arg_match(type.predict)
  type.residuals <- rlang::arg_match(type.residuals)
  df <- augment_columns(x, data = data, newdata = newdata, type.predict = type.predict,type.residuals = type.residuals, se_fit)
  df$.hat <- hatvalues(x, parameters = "mean")
  df$.gencooksd <- cooks.distance(x, type = "GCD", hat = "mean") %>% unname()
  env <- x$envelope
  
  if(is.null(newdata)){
    newdata = model.frame(x)
  }
  
  if(conf_int){
    conf_int <- predict(x, newdata = newdata, interval = "confidence", ...)
    df$.fittedlwrconf <- conf_int[, "lwr"]
    df$.fitteduprconf <- conf_int[, "upr"]
  }
  
  if(pred_int){
    pred_int <- predict(x, newdata = newdata,  interval = "prediction", ...)
    df$.fittedlwrpred <- pred_int[, "lwr"]
    df$.fitteduprpred <-pred_int[, "upr"]
  }
  
  if (!is.null(env)) {
    df$.resenv <- residuals(x, type = x$residualname)
    df$.lwrenv <- env[3,]
    df$.mdnenv <- env[2,]
    df$.uprenv <- env[1,]
  } 
  df
}
  
glance.mixpoissonreg <- function()
  
tidy.mixpoissonreg <- function()
  
autoplot.mixpoissonreg <- function()

local_influence_ggplot.mixpoissonreg <- function()

  

local_influence_ggplot <- function(model, ...){
  UseMethod("local_influence_ggplot", model)
}

tidy_local_influence <- function(model, ...){
  UseMethod("tidy_local_influence", model)
}


