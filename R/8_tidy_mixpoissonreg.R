#' @import broom
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr


augment.mixpoissonreg <- function(x, data = model.frame(x), newdata = NULL, type.predict = c("response", "link", "precision", "variance"), 
                                  type.residuals = c("pearson", "score"), se_fit = FALSE, conf_int = TRUE, pred_int = FALSE, ...) {
  type.predict <- rlang::arg_match(type.predict)
  type.residuals <- rlang::arg_match(type.residuals)
  
  if(is.null(newdata)){
    newdata = model.frame(x)
    res <- residuals(x, type = type.residuals)
  } else{
    res <- NULL
  }
  
  df <- as_tibble(newdata)
  
  pred <- predict(x, newdata = newdata, type = type.predict, se.fit = se_fit)
  

  
  if(se_fit){
    df$.fitted <- pred$fit %>% unname()
    df$.se.fit <- pred$se.fit %>% unname()
  } else{
    df$.fitted <- pred %>% unname()
  }

  
  if(conf_int){
    conf_int <- predict(x, newdata = newdata, interval = "confidence", ...)
    df$.fittedlwrconf <- conf_int[, "lwr"] %>% unname()
    df$.fitteduprconf <- conf_int[, "upr"] %>% unname()
  }
  
  if(pred_int){
    pred_int <- predict(x, newdata = newdata,  interval = "prediction", ...)
    df$.fittedlwrpred <- pred_int[, "lwr"] %>% unname()
    df$.fitteduprpred <-pred_int[, "upr"] %>% unname()
  }
  
  if(!is.null(res)){
    df$.resid <- res 
    df$.hat <- hatvalues(x, parameters = "mean") %>% unname()
    df$.cooksd <- cooks.distance(x, type = "CD", hat = "mean") %>% unname()
    df$.gencooksd <- cooks.distance(x, type = "GCD", hat = "mean") %>% unname()
    
    env <- x$envelope
    
    if (!is.null(env)) {
      df$.resenv <- residuals(x, type = x$residualname)
      df$.lwrenv <- env[3,]
      df$.mdnenv <- env[2,]
      df$.uprenv <- env[1,]
    } 
  }
  
  df
}
  
glance.mixpoissonreg <- function(x, ...){
  tibble(efron.pseudo.r2 = as.numeric(x$efron.pseudo.r2), df.null = x$df.null, 
                   logLik = as.numeric(stats::logLik(x)), AIC = stats::AIC(x), 
                   BIC = stats::BIC(x), df.residual = stats::df.residual(x), 
                   nobs = stats::nobs(x), model.type = x$modeltype, est.method = x$estimation_method)
}
  
tidy.mixpoissonreg <- function(x, conf.int = FALSE, conf.level = 0.95){
  retmean <- as_tibble(summary(x)$coefficients$mean, rownames = "term")
  colnames(retmean) <- c("term", "estimate", "std.error", "statistic", 
                     "p.value")
  retmean <- retmean %>% add_column(component = "mean", .before = "term")
  retprec <- as_tibble(summary(x)$coefficients$precision, rownames = "term")
  colnames(retprec) <- c("term", "estimate", "std.error", "statistic", 
                     "p.value")
  retprec <- retprec %>% add_column(component = "precision", .before = "term")
  
  ret <- bind_rows(retmean,retprec)
  
  if (conf.int) {
    ret$join_term <- names(coef(x, parameters = "all"))
    ci <- as_tibble(confint(x, level = conf.level), rownames = "term")
    names(ci) <- c("join_term", "conf.low", "conf.high")
    ret <- left_join(ret, ci, by = "join_term") %>% select(-join_term)
  }
  ret
}
  
tidy_local_influence <- function(){
    
}
  
autoplot.mixpoissonreg <- function()

local_influence_ggplot.mixpoissonreg <- function()

  

local_influence_ggplot <- function(model, ...){
  UseMethod("local_influence_ggplot", model)
}

tidy_local_influence <- function(model, ...){
  UseMethod("tidy_local_influence", model)
}


