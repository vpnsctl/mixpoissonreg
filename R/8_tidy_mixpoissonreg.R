#' @import generics
#' @import ggplot2

local_influence_ggplot.mixpoissonreg

local_influence_ggplot <- function(model, ...){
  UseMethod("local_influence_ggplot", model)
}