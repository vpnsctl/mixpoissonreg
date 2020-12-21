#' @import broom
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import ggfortify
#' @import gridExtra
#' @import grid

#############################################################################################
#' @title augment.mixpoissonreg
#' @description 
#' @param 
#' @return 

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
    df$.resfit <- residuals(x, type = x$residualname)
    df$.hat <- hatvalues(x, parameters = "mean") %>% unname()
    df$.cooksd <- cooks.distance(x, type = "CD", hat = "mean") %>% unname()
    df$.gencooksd <- cooks.distance(x, type = "GCD", hat = "mean") %>% unname()
    
    env <- x$envelope
    
    if (!is.null(env)) {
      df$.index <- 1:nrow(df)
      temp_tbl <- tibble(.resfit = df$.resfit, .index = df$.index)
      temp_tbl <- temp_tbl %>% arrange(.resfit)
      temp_tbl$.lwrenv <- env[3,]
      temp_tbl$.mdnenv <- env[2,]
      temp_tbl$.uprenv = env[1,]
      df <- left_join(df, temp_tbl, by = c(".index", ".resfit")) %>% select(-.index)
    } 
  }
  
  df
}
  
#############################################################################################
#' @title glance.mixpoissonreg
#' @description 
#' @param 
#' @return 

glance.mixpoissonreg <- function(x, ...){
  tibble(efron.pseudo.r2 = as.numeric(x$efron.pseudo.r2), df.null = x$df.null, 
                   logLik = as.numeric(stats::logLik(x)), AIC = stats::AIC(x), 
                   BIC = stats::BIC(x), df.residual = stats::df.residual(x), 
                   nobs = stats::nobs(x), model.type = x$modeltype, est.method = x$estimation_method)
}

#############################################################################################
#' @title tidy.mixpoissonreg
#' @description 
#' @param 
#' @return   

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

#############################################################################################
#' @title tidy_local_influence.mixpoissonreg
#' @description 
#' @param 
#' @return 
  
tidy_local_influence.mixpoissonreg <- function(x, perturbation = c("case_weights", "hidden_variable",
                                                     "mean_explanatory", "precision_explanatory",
                                                     "simultaneous_explanatory"), curvature = c("conformal", "normal"),
                                 direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                 mean.covariates = NULL, precision.covariates = NULL){
    loc_infl <- local_influence(x, perturbation = perturbation, curvature = curvature, direction = direction, 
                                parameters = parameters, mean.covariates = mean.covariates, precision.covariates = precision.covariates)
  
    tidy_loc_infl <- tibble(.rows = nobs(x))
    
    for(pert in perturbation){
      tidy_loc_infl = tidy_loc_infl %>% add_column(!!pert := loc_infl[[pert]])
    }
    tidy_loc_infl
}

#############################################################################################
#' @title local_influence_benchmarks.mixpoissonreg
#' @description 
#' @param 
#' @return 

local_influence_benchmarks.mixpoissonreg <- function(x, perturbation = c("case_weights", "hidden_variable",
                                                           "mean_explanatory", "precision_explanatory",
                                                           "simultaneous_explanatory"), curvature = c("conformal", "normal"),
                                       direction = c("canonical", "max.eigen"), parameters = c("all", "mean", "precision"),
                                       mean.covariates = NULL, precision.covariates = NULL){
  loc_infl <- local_influence(x, perturbation = perturbation, curvature = curvature, direction = direction, 
                              parameters = parameters, mean.covariates = mean.covariates, precision.covariates = precision.covariates)
  benchmarks <- c()
  for(pert in perturbation){
    benchmarks <- c(benchmarks, attr(loc_infl[[pert]], "benchmark"))
  }
  benchmarks <- matrix(benchmarks, nrow = 1)
  colnames(benchmarks) <- perturbation
  benchmarks_tbl <- as_tibble(benchmarks)
  benchmarks_tbl
}

#############################################################################################
#' @title autoplot.mixpoissonreg
#' @description 
#' @param label.label vector of labels.
#' @param env_alpha alpha of the bands of the envelope (when the fitted model has envelopes)
#' @param env_fill the colour of the filling in the envelopes.
#' @param gpar_sub.caption list of gpar parameters to be used as common title in the case of multiple plots. The title will be given in sub.caption argument. See
#' the help of \code{gpar} function from the \pkg{grid} package for all the available options.
#' @return 
#' @details Based on \code{autoplot.glm} from the excellent \code{ggfortify} package, \href{https://github.com/sinhrks/ggfortify}.
  
autoplot.mixpoissonreg <- function(object, which = c(1,2,5,6), title = list("Residuals vs Obs. number",
                                                                               "Normal Q-Q",
                                                                               "Cook's distance",
                                                                               "Generalized Cook's distance",
                                                                               "Cook's dist vs Generalized Cook's dist",
                                                                               "Response vs Fitted means"),
                                   label.repel = TRUE,
                                   nrow = NULL, ncol = NULL,
                                    qqline = TRUE, ask = prod(par("mfcol")) <
                                      length(which) && dev.interactive(), include.modeltype = TRUE,
                                    include.residualtype = FALSE, sub.caption = NULL,
                                   env_alpha = 0.5, env_fill = "grey70", gpar_sub.caption = list(fontface = "bold"),
                                    colour = "#444444", size = NULL, linetype = NULL, alpha = NULL, fill = NULL, 
                                    shape = NULL, label = TRUE, label.label = NULL, label.colour = "#000000", 
                                    label.alpha = NULL, label.size = NULL, label.angle = NULL, 
                                    label.family = NULL, label.fontface = NULL, label.lineheight = NULL, 
                                    label.hjust = NULL, label.vjust = NULL, 
                                    label.n = 3, ad.colour = "#888888", ad.linetype = "dashed", ad.size = 0.2, ...){
  p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- NULL
  if (is.null(sub.caption)) {
    cal <- object$call
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
  
  plot.data <- augment(object, type.residuals = object$residualname, type.predict = "response")
  n <- nobs(object)
  plot.data$.index <- 1:n
  if(is.null(label.label)){
    plot.data$.label <- rownames(plot.data)
  } else{
    plot.data$.label <- as.vector(label.label)
  }
  
  if(is.null(object$y)){
    y = object$residuals + object$fitted.values
  } else{
    y = object$y
  }
  
  if(2 %in% which){
    ylim <- range(plot.data$.resid, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    qx <- stats::qqnorm(plot.data$.resid, ylim = ylim, 
                        plot.it = FALSE)$x
    plot.data$.qqx <- qx
    
    qprobs <- c(0.25, 0.75)
    qy <- stats::quantile(plot.data$.resid, probs = qprobs, 
                          names = FALSE, type = 7, na.rm = TRUE)
    qx <- stats::qnorm(qprobs)
    slope <- diff(qy)/diff(qx)
    int <- qy[1L] - slope * qx[1L]
  }
  
  plot.data$.obs <- as.vector(y)
  
  plot_label <- function (p, data, x = NULL, y = NULL, label = TRUE, label.label = "rownames", 
            label.colour = NULL, label.alpha = NULL, label.size = NULL, 
            label.angle = NULL, label.family = NULL, label.fontface = NULL, 
            label.lineheight = NULL, label.hjust = NULL, label.vjust = NULL, 
            label.repel = FALSE, label.show.legend = NA) 
  {
    if (!is.data.frame(data)) {
      stop(paste0("Unsupported class: ", class(data)))
    }
    if (!missing(label.colour) && !is.null(label.colour) && missing(label)) {
      label <- TRUE
    }
    if (label || label.repel) {
      if (is.null(label.colour)) {
        label.colour <- "#000000"
      }
      if (label.repel && "ggrepel" %in% rownames(installed.packages())) {
        textfunc <- ggrepel::geom_text_repel
      }
      else {
        textfunc <- ggplot2::geom_text
      }
      p <- p + geom_factory(textfunc, data, x = x, y = y, label = ".label", 
                            colour = label.colour, alpha = label.alpha, size = label.size, 
                            angle = label.angle, family = label.family, fontface = label.fontface, 
                            lineheight = label.lineheight, hjust = label.hjust, 
                            vjust = label.vjust, show.legend = label.show.legend)
    }
    p
  }
  
  flatten <- function (df) 
  {
    ismatrix <- vapply(df, is.matrix, logical(1))
    if (any(ismatrix)) {
      return(data.frame(c(df[!ismatrix], do.call(data.frame, 
                                                 df[ismatrix])), stringsAsFactors = FALSE))
    }
    else {
      return(df)
    }
  }
  
  geom_factory <- function (geomfunc, data = NULL, ...) 
  {
    mapping <- list()
    option <- list()
    columns <- colnames(data)
    for (key in names(list(...))) {
      value <- list(...)[[key]]
      if (is.null(value)) {
      }
      else if (value %in% columns) {
        mapping[[key]] <- value
      }
      else {
        option[[key]] <- value
      }
    }
    if (!is.null(data)) {
      option[["data"]] <- data
    }
    option[["mapping"]] <- do.call(ggplot2::aes_string, mapping)
    return(do.call(geomfunc, option))
  }
  
  if (is.logical(shape) && !shape) {
    if (missing(label)) {
      label <- TRUE
    }
    if (missing(label.n)) {
      label.n <- nrow(plot.data)
    }
  }
  
  dev_ask <- is.null(nrow) & is.null(ncol)
  
  if(dev_ask){
    one.fig <- prod(par("mfcol")) == 1
    
    if (ask) {
      oask <- grDevices::devAskNewPage(TRUE)
      on.exit(grDevices::devAskNewPage(oask))
    }
  }
  
  plot.data <- flatten(plot.data)
  if(any(c(1,2) %in% which)){
    residualname <- paste0(toupper(substring(object$residualname, 1, 1)), substring(object$residualname, 2))
  }
  
  if (label.n > 0L) {
    if (any(c(1,2,6) %in% which)) {
      r.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.resid)))
      r.data <- utils::head(r.data, label.n)
    }
    if (3 %in% which) {
      cd.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.cooksd)))
      cd.data <- utils::head(cd.data, label.n)
    }
    if (4 %in% which){
      gcd.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.gencooksd)))
      gcd.data <- utils::head(gcd.data, label.n)
    }
    if(5 %in% which){
      cdgcd.data <- dplyr::arrange(plot.data, dplyr::desc((.gencooksd/sum(.gencooksd))^2 + (.cooksd/sum(.cooksd))^2))
      cdgcd.data <- utils::head(cdgcd.data, label.n)
    }
    if(6 %in% which){
      respobs.data <- dplyr::arrange(plot.data, dplyr::desc(abs(.obs - .fitted)))
      respobs.data <- utils::head(respobs.data, label.n)
    }
  }
  
  .decorate.label <- function(p, data) {
    if (label & label.n > 0) {
      p <- plot_label(p = p, data = data, label = label, 
                      label.label = ".label", label.colour = label.colour, 
                      label.alpha = label.alpha, label.size = label.size, 
                      label.angle = label.angle, label.family = label.family, 
                      label.fontface = label.fontface, label.lineheight = label.lineheight, 
                      label.hjust = label.hjust, label.vjust = label.vjust, 
                      label.repel = label.repel)
    }
    p
  }
  
  .decorate.plot <- function(p, xlab = NULL, ylab = NULL, title = NULL) {
    p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
  }
  
  if (1 %in% which) {
    if(include.residualtype){
      title[[1]] = paste(residualname, title[[1]])
    }
    
    ylab <- paste0(residualname, " residuals")
    t1 <- title[[1]]
    mapping <- ggplot2::aes_string(x = ".index", y = ".resid")
    p1 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p1 <- p1 + geom_factory(geom_point, plot.data, colour = colour, 
                              size = size, linetype = linetype, alpha = alpha, 
                              fill = fill, shape = shape)
    }
    p1 <- p1 + ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype, 
                          size = ad.size, colour = ad.colour)
    p1 <- .decorate.label(p1, r.data)
    
    if(sub.caption == "" | !dev_ask){
      xlab = "Obs. number"
    } else {
      xlab = paste0("Obs. number\n",sub.caption)
    }
    p1 <- .decorate.plot(p1, xlab = xlab, ylab = ylab, 
                         title = t1)
    if(dev_ask){
      print(p1)
    }
  }
  
  if(2 %in% which){
    env <- object$envelope
    if (!is.null(env)) {
      title[[2]] <- paste0(title[[2]]," with simulated envelopes")
    } 
    t2 <- title[[2]]

    mapping <- ggplot2::aes_string(x = ".qqx", y = ".resid")
    p2 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    
    if(qqline & is.null(env)){
      p2 <- p2 + ggplot2::geom_abline(intercept = int, slope = slope, 
                                      linetype = ad.linetype, size = ad.size, colour = ad.colour)
    }
    
    if(!is.null(env)){
      p2 <- p2 + ggplot2::geom_ribbon(aes(ymin=.lwrenv, ymax=.uprenv), alpha=env_alpha, fill = env_fill)
      p2 <- p2 + ggplot2::geom_line(aes(y = .mdnenv), linetype = ad.linetype, size = ad.size, colour = ad.colour)
    }
    
    if (!is.logical(shape) || shape) {
      p2 <- p2 + geom_factory(geom_point, plot.data, colour = colour, 
                              size = size, linetype = linetype, alpha = alpha, 
                              fill = fill, shape = shape)
    }
    
    ylab <- paste0(residualname, " residuals")
    
    if(sub.caption == "" | !dev_ask){
      xlab = "Theoretical Quantiles"
    } else {
      xlab = paste0("Theoretical Quantiles\n", sub.caption)
    }
    
    p2 <- .decorate.label(p2, r.data)
    p2 <- .decorate.plot(p2, xlab = xlab, 
                         ylab = ylab, title = t2)
    if(dev_ask){
      print(p2)
    }
  }
  
  if(3 %in% which){
    t3 <- title[[3]]
    mapping <- ggplot2::aes_string(x = ".index", y = ".cooksd", 
                                   ymin = 0, ymax = ".cooksd")
    p3 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p3 <- p3 + geom_factory(geom_linerange, plot.data, 
                              colour = colour, size = size, linetype = linetype, 
                              alpha = alpha, fill = fill, shape = shape)
    }
    p3 <- .decorate.label(p3, cd.data)
    
    if(sub.caption == "" | !dev_ask){
      xlab = "Obs. Number"
    } else {
      xlab = paste0("Obs. Number\n", sub.caption)
    }
    
    p3 <- .decorate.plot(p3, xlab = xlab, ylab = "Cook's distance", 
                         title = t3)
    if(dev_ask){
      print(p3)
    }
  }
  
  if(4 %in% which){
    t4 <- title[[4]]
    mapping <- ggplot2::aes_string(x = ".index", y = ".gencooksd", 
                                   ymin = 0, ymax = ".gencooksd")
    p4 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p4 <- p4 + geom_factory(geom_linerange, plot.data, 
                              colour = colour, size = size, linetype = linetype, 
                              alpha = alpha, fill = fill, shape = shape)
    }
    p4 <- .decorate.label(p4, gcd.data)
    
    if(sub.caption == "" | !dev_ask){
      xlab = "Obs. Number"
    } else {
      xlab = paste0("Obs. Number\n", sub.caption)
    }
    
    p4 <- .decorate.plot(p4, xlab = xlab, ylab = "Generalized Cook's distance", 
                         title = t4)
    if(dev_ask){
      print(p4)
    }
  }
  
  if (5 %in% which) {
    t5 <- title[[5]]
    mapping <- ggplot2::aes_string(x = ".gencooksd", y = ".cooksd")
    p5 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p5 <- p5 + geom_factory(geom_point, plot.data, colour = colour, 
                              size = size, linetype = linetype, alpha = alpha, 
                              fill = fill, shape = shape)
    }
    p5 <- p5 + ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype, 
                                   size = ad.size, colour = ad.colour)
    p5 <- .decorate.label(p5, cdgcd.data)
    
    if(sub.caption == "" | !dev_ask){
      xlab = "Generalized Cook's distance"
    } else {
      xlab = paste0("Generalized Cook's distance\n", sub.caption)
    }
    
    p5 <- .decorate.plot(p5, xlab = xlab, ylab = "Cook's distance", 
                         title = t5)
    if(dev_ask){
      print(p5)
    }
  }
  
  if (6 %in% which) {
    t6 <- title[[6]]
    mapping <- ggplot2::aes_string(x = ".fitted", y = ".obs")
    p6 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p6 <- p6 + geom_factory(geom_point, plot.data, colour = colour, 
                              size = size, linetype = linetype, alpha = alpha, 
                              fill = fill, shape = shape)
    }
    p6 <- p6 + ggplot2::geom_abline(intercept = 0L, slope = 1L, linetype = ad.linetype, 
                                   size = ad.size, colour = ad.colour)
    p6 <- .decorate.label(p6, respobs.data)
    
    if(sub.caption == "" | !dev_ask){
      xlab = "Predicted values"
    } else {
      xlab = paste0("Predicted values\n", sub.caption)
    }
    
    p6 <- .decorate.plot(p6, xlab = xlab, ylab = "Response values", 
                         title = t6)
    if(dev_ask){
      print(p6)
    }
  }
  

  if(!dev_ask){
    plot.list <- list(p1, p2, p3, p4, p5, p6)[which]
    if(is.null(nrow))
      nrow <- 0
    if(is.null(ncol))
      ncol <- 0
    p <- new("ggmultiplot", plots = plot.list, nrow = nrow, ncol = ncol)
    
    gpar <- do.call(grid::gpar, gpar_sub.caption)
    title_multi <- grid::textGrob(sub.caption, gp=gpar)
    gridExtra::grid.arrange(grobs = p@plots, top = title_multi, gp=gpar(fontface="bold"))
  } else{
    invisible()
  }
}

#############################################################################################
#' @title local_influence_autoplot.mixpoissonreg
#' @description 
#' @param 
#' @return 

local_influence_autoplot.mixpoissonreg <- function(){
  
}

#############################################################################################
#' @title local_influence_autoplot
#' @description 
#' @param 
#' @return   

local_influence_autoplot <- function(model, ...){
  UseMethod("local_influence_ggplot", model)
}

#############################################################################################
#' @title tidy_local_influence
#' @description 
#' @param 
#' @return  

tidy_local_influence <- function(model, ...){
  UseMethod("tidy_local_influence", model)
}

#############################################################################################
#' @title local_influence_benchmarks
#' @description 
#' @param 
#' @return  

local_influence_benchmarks <- function(model, ...){
  UseMethod("local_influence_benchmarks", model)
}

