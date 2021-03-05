#' @import Formula
#############################################################################################
#' @name mixpoissonreg
#' @title Mixed Poisson Regression for Overdispersed Count Data
#' @aliases mixpoissonreg mixpoissonreg.fit
#' @description Fits mixed Poisson regression models (Poisson-Inverse Gaussian or Negative-Binomial) on data sets with response variables being count data.
#' The models can have varying precision parameter, where a linear regression structure (through a link function) is assumed to hold on the precision parameter.
#' The Expectation-Maximization algorithm for both these models (Poisson Inverse Gaussian and Negative Binomial) is an important contribution of this package.
#' Another important feature of this package is the set of functions to perform global and local influence analysis.
#' @param formula symbolic description of the model (examples: \code{y ~ x1 + ... + xnbeta} and \code{y ~ x1 + ... + xnbeta | w1 + ... + wnalpha}); see details below.
#' @param data elements expressed in formula. This is usually a data frame composed by:
#' (i) the observations formed by count data \code{z}, with z_i being non-negative integers,
#' (ii) covariates for the mean submodel (columns \code{x1, ..., xnbeta}) and
#' (iii) covariates for the precision submodel (columns \code{w1, ..., wnalphla}).
#' @param model character ("NB" or "PIG") indicating the type of model to be fitted, with
#' "NB" standing for Negative-Binomial and "PIG" standing for Poisson Inverse Gaussian. The default is "NB".
#' @param y For \code{mixpoissonreg}: logical values indicating if the response vector should be returned as component.
#'
#' For \code{mixpoissonreg.fit}: a numerical vector of response variables with length \code{n}. Each coordinate must be a nonnegative-integer.
#' @param x For \code{mixpoissonreg}: logical values indicating if the model matrix \code{x} should be returned as component.
#'
#' For \code{mixpoissonreg.fit}: a matrix of covariates with respect to the mean with dimension \code{(n,nbeta)}.
#' @param w For \code{mixpoissonreg}: logical values indicating if the model matrix \code{w} should be returned as component.
#'
#' For \code{mixpoissonreg.fit} a matrix of covariates with respect to the precision parameter. The default is \code{NULL}.
#' If not \code{NULL} must be of dimension \code{(n,nalpha)}.
#' @param method estimation method to be chosen between "EM" (Expectation-Maximization) and "ML" (Maximum-Likelihood). The default method is "EM".
#' @param residual character indicating the type of residual to be evaluated ("pearson" or "score"). The default is "pearson". Notice that they coincide for Negative-Binomial models.
#' @param envelope number of simulations (synthetic data sets) to build envelopes for residuals (with \code{100*prob\%} confidence level).
#' The default \code{envelope = 0} dismisses the envelope analysis.
#' @param prob probability indicating the confidence level for the envelopes (default: \code{prob} = 0.95).
#' If \code{envelope} = 0, \code{prob} is ignored.
#' @param model.frame logical indicating whether the model frame should be returned as component of the returned value.
#' @param link.mean optionally, a string containing the link function for the mean. If omitted, the 'log' link function will be used.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision optionally, a string containing the link function the precision parameter. If omitted and the only precision
#' covariate is the intercept, the 'identity' link function will be used, if omitted and there is a precision covariate other than the
#' intercept, the 'log' link function will be used. The possible link functions for the precision parameter are "identity" and "inverse.sqrt" (which is \eqn{\phi^{-1/2} = w_i^T alpha}).
#' @param em_controls only used with the 'EM' method. A list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm,
#' the default is set to 5000;
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5);
#' \code{em_tolgrad} that defines the tolerance value
#' of the maximum-norm of the the gradient of the Q-function, the default is set to 10^(-2).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return \code{mixpoissonreg} returns an object of class "mixpoissonreg" whereas \code{mixpoissonreg.fit}
#' returns an object of class "mixpoissonreg_fit". Both objects are given by lists containing the outputs from the model fit (Negative-Binomial or Poisson Inverse Gaussian regression).
#'
#' An object of the class "mixpoissonreg" is a list containing the following elements:
#' \itemize{
#'   \item \code{coefficients} - a list with elements "mean" and "precision" containing the estimated coefficients of the model;
#'   \item \code{call} - the formula used by the model. If using \code{mixpoissonreg.fit}, this returns \code{NULL}.
#'   \item \code{modelname} - the fitted model, NB or PIG;
#'   \item \code{modeltype} - the abbreviated model name
#'   \item \code{residualname} - the name of the chosen residual in the call, 'pearson' or 'score';
#'   \item \code{niter} - number of iterations of the EM algorithm if method = "EM" and number of iterations
#'   of the \code{optim} function, if method = "ML";
#'   \item \code{start} - the initial guesses of the parameters
#'   \item \code{intercept} - vector indicating if the intercept is present in the mean and/or in the precision regressions;
#'   \item \code{link.mean} - link function of the mean;
#'   \item \code{link.precision} - link function of the precision parameter;
#'   \item \code{fitted.values} - a vector of fitted values in the response scale;
#'   \item \code{fitted.precisions} - a vector of fitted precisions;
#'   \item \code{efron.pseudo.r2} - Efron's pseudo R^2: the squared correlation between the response variables and the predicted values;
#'   \item \code{vcov} - covariance matrix of the parameters of the fitted model;
#'   \item \code{logLik} - log-likelihood at the estimated parameters;
#'   \item \code{Qfunction} - Q-function at the estimated parameters;
#'   \item \code{x} - the covariates related to the mean (if x = TRUE);
#'   \item \code{w} - the covariates related to the precision parameter (if w = TRUE);
#'   \item \code{y} - the response variables (if y = TRUE);
#'   \item \code{model} - if requested (the default), the model frame;
#'   \item \code{formula} - the formula supplied;
#'   \item \code{nobs} - number of observations
#'   \item \code{df.null} - the residual degrees of freedom for the model with constant mean and constant precision;
#'   \item \code{df.residual} - the residual degrees of freedom of the fitted model;
#'   \item \code{estimation_method} - the estimation method, "EM" or "ML"
#'   \item \code{residuals} - vector of raw residuals, that is, the response variable minus the fitted means;
#'   \item \code{std_errors} - the standard errors of the estimated parameters;
#'   \item \code{envelope} - the numerical envelopes used to build the Q-Q plot with simulated envelopes;
#'   \item \code{terms} - (only for \code{mixpoissonreg})the \code{terms} object used;
#'   \item \code{levels} - (where relevant, only for \code{mixpoissonreg}) the levels of the factors used;
#'   \item \code{contrasts} - (where relevant, only for \code{mixpoissonreg}) the contrasts used.
#' }
#'
#' @details Among the regression models with discrete response variables, Poisson regression is the most popular
#' for modeling count data. See, for instance Sellers and Shmueli (2010).
#' It is well-known that this model is equidispersed (that is, the mean is equal to the variance),
#' which in practice may be an unrealistic
#' assumption. Several models have been introduced in the literature to overcome this problem such as
#' negative binomial (NB) and Poisson inverse gaussian (PIG) distributions (see Lawless, 1987).
#' The most common way to do this is to consider a mixed Poisson distribution, which is defined as follows.
#' Let \eqn{Z} be a positive random variable (generally being continuous) with distribution
#' function \eqn{G_{\tau}(\cdot)},
#' where \eqn{\tau} denotes the parameter vector associated to the \eqn{G} distribution. Let
#' \eqn{Y|Z=z\sim}Poisson\eqn{(\mu z)}, for
#' some constant \eqn{\mu>0}. Therefore \eqn{Y} follows a mixed Poisson (MP) distribution with probability
#' function given by
#' \deqn{P(Y=y)=\int_0^\infty\frac{e^{-\mu z}(\mu z)^y}{y!}dG_{\tau}(z),}
#' for \eqn{y=0,1,\ldots}. With this,
#' \eqn{Y} has an overdispersed distribution and hence it is a natural alternative to the Poisson distribution.
#' The most common choices for \eqn{Z} are gamma and inverse-gaussian distributions,
#' which yields \eqn{Y} following, respectively, NB and PIG distributions.
#' General properties of the MP distributions can be found in Karlis and Xekalaki (2005) and in the references therein.
#'
#' In \code{mixpoissonreg} two regression models are implemented, namely, the NB and PIG regression models.
#' We follow the definitions and notations given in Barreto-Souza and Simas (2016). The mixed Poisson regression model
#' is defined by assuming \eqn{Y_1,\ldots,Y_n} is a random sample where
#' \eqn{Y_i\sim NB(\mu_i,\phi_i)} or \eqn{Y_i\sim PIG(\mu_i,\phi_i)} for \eqn{i = 1,\ldots,n}.
#' Under this parameterization we have \eqn{E(Y_i) = \mu_i} and \eqn{Var(Y_i) = \mu_i(1+\mu_i\phi_i^{-1}b''(\xi_0))}, where
#' \eqn{b(\theta) = -\log(-\theta)} and \eqn{\xi_0 = -1} for the NB case, and \eqn{b(\theta) = -(-2\theta)^{1/2}} and \eqn{\xi_0 = -1/2} for
#' the PIG case, with \eqn{b''(\cdot)} being the second derivative of the function \eqn{b(\cdot)}.
#' The following linear relations are assumed
#' \deqn{\Lambda_1(\mu_i) = x_i^T \beta}
#'  and
#' \deqn{\Lambda_2(\phi_i) = w_i^T \alpha,}
#'  where \eqn{\beta = (\beta_1,...,\beta_p)} and \eqn{\alpha = (\alpha_1,...,\alpha_q)} are real valued vectors.
#' The terms \eqn{x_i^T} and \eqn{v_i^T} represent, respectively, the i-th row of the matrices "x" (\eqn{n\times p})
#'  and "w" (\eqn{n\times q}) containing covariates in their columns
#' (\eqn{x_{i,1}} and \eqn{v_{i,1}} may be 1 to handle intercepts).
#'
#' Therefore, the \code{mixpoissonreg} package handles up to two regression structures
#' at the same time: one for the mean parameter, one for the precision parameter. The regression structure for
#' the mean is determined through a formula \code{y ~ x1 + ... + xn}, whereas the regression structure for
#' the precision parameter is determined through the right-hand side of the formula using the separator "\code{|}". So,
#' for example, a regression with \code{x1,...,xn} as covariates for the mean and \code{z1,...,zm} as covariates for the precision
#' parameter corresponds to the formula \code{y ~ x1 + ... + xn | z1 + ... + zm}. If only there is only formula for
#' the regression structure for the mean, the regression structure for the precision parameter will only have the intercept,
#' that is, \code{y ~ x1 + ... + xn} is the same as \code{y ~ x1 + ... + xn | 1}.
#'
#' In general, in this package, the EM-algorithm estimation method obtains estimates closer to the maximum likelihood estimate than the maximum likelihood estimation method,
#' in the sense that the likelihood function evaluated at the EM-algorithm estimate is greater or equal (usually strictly greater) than the likelihood function evaluated
#' at the maximum likelihood estimate. So, unless the processing time is an issue, we strongly recommend the EM-algorithm as the estimation method.
#'
#' In Barreto-Souza and Simas (2016) two residuals were studied: the pearson residuals
#' and the score residuals. Both these residuals are implemented in the \code{mixpoissonreg}
#' package. They coincide for NB regression models. They can be accessed via
#' the \link[mixpoissonreg:residuals.mixpoissonreg]{residuals} method.
#'
#' It is also noteworthy that all the global and local influence analysis tools developed
#' in Barreto-Souza and Simas (2016) are implemented in this package. See \code{\link{influence.mixpoissonreg}},
#' \code{\link{local_influence.mixpoissonreg}}, \code{\link{local_influence_plot.mixpoissonreg}}
#' and \code{\link{local_influence_autoplot.mixpoissonreg}}.
#'
#' @references
#' DOI:10.1007/s11222-015-9601-6 \doi{10.1007/s11222-015-9601-6}(Barreto-Souza and Simas; 2016)
#'
#' URL:https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1751-5823.2005.tb00250.x (\href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1751-5823.2005.tb00250.x}{Karlis and Xekalaki; 2005})
#'
#' DOI:10.2307/3314912 \doi{10.2307/3314912}(Lawless; 1987)
#'
#' Sellers, K.F. and Shmueli, G. (2010) *A flexible regression model for count data.* Ann. Appl. Stat., 4, 943-961
#'
#' @seealso
#' \code{\link{summary.mixpoissonreg}}, \code{\link{plot.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}},
#' \code{\link{residuals.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}},\code{\link{influence.mixpoissonreg}},
#' \code{\link{cooks.distance.mixpoissonreg}},
#' \code{\link{local_influence.mixpoissonreg}}, \code{\link{local_influence_plot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
#'
#' @examples
#' # Examples using the Attendance dataset:
#' \donttest{
#' daysabs_fit <- mixpoissonreg(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' summary(daysabs_fit)
#' # Base R plot of the fit
#' plot(daysabs_fit)
#' # ggplot2 plot of the fit
#' autoplot(daysabs_fit)
#' # plot of local influence measures
#' local_influence_plot(daysabs_fit)
#' # ggplot2 plot of local influence measures
#' local_influence_autoplot(daysabs_fit)
#' # Fitting a reduced model of the sabe type as the previous one
#' daysabs_fit_red <- mixpoissonreg(daysabs ~ gender + math +
#' prog | prog, data = Attendance, model = daysabs_fit$modeltype)
#' # Likelihood ratio test:
#' lmtest::lrtest(daysabs_fit, daysabs_fit_red)
#' # Wald test:
#' lmtest::waldtest(daysabs_fit, daysabs_fit_red)
#' }
#'
#' @rdname mixpoissonreg
#' @export
mixpoissonreg <- function(formula, data, link.mean = c("log", "sqrt"),
                  link.precision = c("identity", "log", "inverse.sqrt"),
                  model = c("NB", "PIG"), method = c("EM", "ML"),
                  residual = c("pearson", "score"), y = TRUE, x = TRUE, w = TRUE,
                  envelope = 0, prob = 0.95, model.frame = TRUE, em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                  optim_method = "L-BFGS-B", optim_controls = list()) {
  call_mixpoissonreg <- match.call()
  ## Processing call
  # If data is not provided, verify the current R workspace
  if (missing(data)) {
    data <- environment(formula)
  }
  MF <- match.call(expand.dots = FALSE)
  arg_names <- match(c("formula", "data"), names(MF), nomatch = 0)
  MF <- MF[c(1, arg_names)]
  MF$drop.unused.levels <- TRUE

  #Processing em_controls:
  if(is.null(em_controls$maxit)){
    em_controls$maxit = 5000
  }

  if(is.null(em_controls$em_tol)){
    em_controls$em_tol = 10^(-5)
  }

  if(is.null(em_controls$em_tolgrad)){
    em_controls$em_tolgrad = 10^(-2)
  }

  ## Processing formula
  Formula_mpreg <- as.Formula(formula)
  if (length(Formula_mpreg)[2] < 2) {
    Formula_mpreg <- as.Formula(stats::formula(Formula_mpreg), ~1)
  } else {
    if (length(Formula_mpreg)[2] > 2) {
      Formula_mpreg <- Formula(stats::formula(Formula_mpreg, rhs = 1:2))
      warning("right hand side of formula should not have > 2 parts (ignoring 3rd and higher cases)")
    }
  }
  MF$formula <- Formula_mpreg

  ## Model.frame: Convert MF into a data matrix.
  MF[[1]] <- as.name("model.frame")
  MF <- eval(MF, parent.frame())

  model_matrix <- stats::model.frame(Formula_mpreg, data)

  ## Extract terms (covariate matrices and response vector)
  MTerms_x <- stats::terms(Formula_mpreg, data = data, rhs = 1)
  MTerms_w <- stats::delete.response(stats::terms(Formula_mpreg,
                                                  data = data, rhs = 2))
  y_mpreg <- stats::model.response(MF, "numeric")
  x_mpreg <- stats::model.matrix(MTerms_x, MF)
  w_mpreg <- stats::model.matrix(MTerms_w, MF)

  object <- mixpoissonreg.fit(y=y_mpreg,x=x_mpreg,w=w_mpreg,link.mean = link.mean, link.precision = link.precision,
                      model = model, method = method, residual = residual, envelope = envelope,
                      prob = prob, em_controls = em_controls, optim_method = optim_method, optim_controls = optim_controls)

  if(x){
    object$x <- x_mpreg
  }
  if(w){
    object$w <- w_mpreg
  }
  if(y){
    object$y <- y_mpreg
  }

  object$call <- call_mixpoissonreg
  object$terms <- list(mean = MTerms_x, precision = MTerms_w)
  object$levels <- list(mean = stats::.getXlevels(MTerms_x, MF),
                        precision = stats::.getXlevels(MTerms_w, MF))
  object$contrasts <- list(mean = attr(x, "contrasts"),
                           precision = attr(w, "contrasts"))
  object$formula <- Formula_mpreg

  if(model.frame){
    object$model <- model_matrix
  }

  class(object) <- "mixpoissonreg"
  return(object)

}


#' @rdname mixpoissonreg
#' @export

mixpoissonreg.fit <- function(x, y, w = NULL, link.mean = c("log", "sqrt"),
                      link.precision = c("identity", "log", "inverse.sqrt"),
                      model = c("NB", "PIG"), method = c("EM", "ML"),
                      residual = c("pearson", "score"), envelope = 0,
                      prob = 0.95, em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)), optim_method = "L-BFGS-B", optim_controls = list()) {

  #Processing em_controls:
  if(is.null(em_controls$maxit)){
    em_controls$maxit = 5000
  }

  if(is.null(em_controls$em_tol)){
    em_controls$em_tol = 10^(-5)
  }

  if(is.null(em_controls$em_tolgrad)){
    em_controls$em_tolgrad = 10^(-2)
  }

  n <- length(y)
  x = as.matrix(x)
  if(is.null(w)){
    w = matrix(rep(1,n), nrow = n)
  } else{
    w = as.matrix(w)
  }
  nbeta <- ncol(x)
  nalpha <- ncol(w)

  ## Check for intercepts:
  intercept_x <- sum(colSums(x==1)==n) > 0
  intercept_w <- sum(colSums(w==1)==n) > 0
  intercept <- c(intercept_x, intercept_w)


  if (nbeta== 0) {
    stop("empty matrix x is not allowed")
  }
  if (nalpha == 0) {
    stop("empty matrix w is not allowed")
  }

  ## Validation of input variables and arguments.
  if (n < 1) {
    stop("response variable is not provided")
  }
  if (min(y) < 0) {
    stop("requirement 0 <= y_i is not true for all responses")
  }
  if (sum(sapply(y, function(x){(x-round(x))!=0})) > 0){
    stop("the response variables must be integers")
  }
  #
  if (length(residual) > 1) {
    residual <- residual[1]
  }
  if (is.character(residual) == FALSE) {
    stop("residual must be a character (pearson, score)")
  }
  if (is.character(residual) == TRUE) {
    residual <- match.arg(residual, c("pearson", "score"))
  }
  #
  if (is.numeric(envelope) == FALSE) {
    stop("argument envelope must be numeric (0 or positive interger)")
  }
  if (envelope < 0 | envelope %% 1 != 0) {
    stop("consider 0 or positive integers for envelope")
  }
  if (envelope > 0 & envelope < 10) {
    stop("number of simulations to build envelopes is too small (try at least envelope = 10)")
  }
  if (is.numeric(prob) == FALSE) {
    stop("prob must be numeric in the interval (0,1)")
  }
  if (prob <= 0 | prob >= 1) {
    stop("prob must be in the interval (0,1)")
  }
  if (envelope == 0 & prob != 0.95) {
    warning(paste0("prob = ", prob, " is ignored since envelope = 0"))
  }
  #
  model <- rlang::arg_match(model)
  method <- rlang::arg_match(method)


  #

  ## Checking and processing link functions
  link.mean <- rlang::arg_match(link.mean)

  possible_link_mean <- c("log", "sqrt")

  if (!(link.mean %in% possible_link_mean)) {
    stop(paste0("link function for the mean must be one of ", possible_link_mean))
  }

  if (length(link.precision) > 1) {
    if (nalpha == 1 & intercept_w) {
      link.precision <- link.precision[1]
    } else {
      link.precision <- link.precision[2]
    }
  }
  possible_link_precision <- c("identity", "log", "inverse.sqrt")

  if (!(link.precision %in% possible_link_precision)) {
    stop(paste0("link function for the precision parameter must be one of ", possible_link_precision))
  }

  ## Set initial values.



  if (is.character(model) == TRUE) {
    aux_model <- match.arg(model, c("NB", "PIG"))
      start <- startvalues_mpreg(y, x, w, link.mean, link.precision, aux_model)
      beta <- start$beta
      alpha <- start$alpha
      start <- c(beta, alpha)
      names(start) = c(colnames(x), paste(colnames(w),".precision", sep = ""))

      modelname <- switch(aux_model,
                          "NB" = {"Negative Binomial Regression"},
                          "PIG" = {"Poisson Inverse Gaussian Regression"}
      )
      inits <- list(beta, alpha)
      if(method == "EM"){
        new_start <- tryCatch(suppressWarnings(ML_mixpoisson(beta, alpha, y, x, w,
                                   link.mean, link.precision, aux_model, optim_method, optim_controls)$coefficients),
                              error = function(e) {
                                "Error"
                              })
        if(length(new_start) > 1){
          beta <- new_start$mean
          alpha <- new_start$precision
        }

        fitted_mpreg <- EM_mixpoisson(beta, alpha, y, x, w,
                                  link.mean, link.precision, aux_model, em_controls, optim_method, optim_controls)
      } else if(method == "ML"){
        fitted_mpreg <- ML_mixpoisson(beta, alpha, y, x, w,
                                      link.mean, link.precision, aux_model, optim_method, optim_controls)
      } else{
        stop("method must be 'EM' or 'ML'")
      }


      niter <- fitted_mpreg$niter # number of iterations of the EM algorithm
      coefficients_mpreg <- fitted_mpreg$coefficients # coefficients
      fitted_values <- fitted_mpreg$fitted.values # mu
      residuals_mpreg <- switch(residual,
                    "pearson" = {pearson_residual_mixpoisson(coefficients_mpreg, y, x, w,
                                         link.mean, link.precision, aux_model)},
                    "score" = {score_residual_mixpoisson(coefficients_mpreg, y, x, w,
                                                           link.mean, link.precision, aux_model)}
      )
      SEr <- std_error_mixpoisson(coefficients_mpreg, y, x, w,
                             link.mean, link.precision, aux_model) # standard errors
      if(sum(is.nan(SEr)) > 0){
        warning("Please, try another optimization method.")
      }

      Env <- NULL
      if (envelope > 0) {
        Env <- envelope_mixpoisson(residual, method,
                                   coefficients_mpreg, x, w, envelope,
                                   prob, n, link.mean,
                                   link.precision,aux_model,
                                   em_controls, optim_method, optim_controls)
        # % of residuals inside the envelope
        Env_prop <- 100 * sum(sort(residuals_mpreg) < Env[1, ] & sort(residuals_mpreg) > Env[3, ]) / n
      }
  }

  #vcov

  fisher_obs <- obs_fisher_information_mixpoisson(c(coefficients_mpreg$mean,coefficients_mpreg$precision), y, x,
                                                  w, link.mean,
                                                  link.precision,
                                                  model)

  vcov_mpreg <- tryCatch(solve(fisher_obs), error = function(e) matrix(NA, nrow(fisher_obs), ncol(fisher_obs)))

  #Efron pseudo R^2: the squared linear correlation between the response variable and the estimated mean
  efron.pseudo.r2 <- stats::cor(y, fitted_values)^2

  #Naming variables accordingly
  names(coefficients_mpreg$mean) <- colnames(x)
  names(coefficients_mpreg$precision) <- colnames(w)
  names(fitted_values) <- 1:length(fitted_values)

  #residuals
  raw_residuals <- y - fitted_values
  raw_residuals <- c(raw_residuals)
  names(raw_residuals) <- 1:length(raw_residuals)

  #precisions
  fitted_precisions <- c(fitted_mpreg$fitted.precisions)
  names(fitted_precisions) <- 1:length(fitted_precisions)

  #log-likelihood
  logLik_mpreg <- loglik_mixpoisson(c(coefficients_mpreg$mean,coefficients_mpreg$precision), y, x,
                                    w, link.mean,
                                    link.precision,
                                    model)

  #Q-Function

  Qfunction_mpreg <- Q_function_mixpoisson(c(coefficients_mpreg$mean,coefficients_mpreg$precision), fitted_values, fitted_precisions, y, x,
                                           w, link.mean,
                                           link.precision,
                                           model)

  ###############
  object <- list()
  object$call <- NULL
  object$modelname <- modelname
  object$modeltype <- model
  object$residualname <- residual
  object$niter <- niter
  object$start <- start
  object$intercept <- intercept
  object$link.mean <- link.mean
  object$link.precision <- link.precision
  object$fitted.values <- fitted_values
  object$fitted.precisions <- fitted_precisions
  object$efron.pseudo.r2 <- efron.pseudo.r2
  object$vcov <- vcov_mpreg
  object$nobs <- length(y)
  object$df.null <- length(y) - 2
  object$df.residual <- length(y) - length(coefficients_mpreg$mean) - length(coefficients_mpreg$precision)
  object$logLik <- logLik_mpreg
  object$Qfunction <- Qfunction_mpreg
  object$coefficients <- coefficients_mpreg
  object$start <- start
  object$estimation_method <- method
  object$residuals <- raw_residuals
  object$std_errors <- SEr
  object$envelope <- Env
  if (is.null(Env) == FALSE) {
    object$envelope_prop <- Env_prop
  }
  ###############
  class(object) <- "mixpoissonreg_fit"
  return(object)
}


#############################################################################################
#' @name mixpoissonregML
#' @aliases mixpoissonregML mixpoissonregML.fit
#' @title Maximum Likelihood Mixed Poisson Regression Models for Overdispersed Count Data
#' @description Uses maximum likelihood estimators to fit mixed Poisson regression models (Poisson-Inverse Gaussian or Negative-Binomial) on data sets with response variables being count data. The models can have varying precision parameter, where a linear regression structure (through a link function) is assumed to hold on the precision parameter.
#' @param formula symbolic description of the model (examples: \code{y ~ x1 + ... + xnbeta} and \code{y ~ x1 + ... + xnbeta | w1 + ... + wnalpha}); see details below.
#' @param data elements expressed in formula. This is usually a data frame composed by:
#' (i) the observations formed by count data \code{z}, with z_i being non-negative integers,
#' (ii) covariates for the mean submodel (columns \code{x1, ..., xnbeta}) and
#' (iii) covariates for the precision submodel (columns \code{w1, ..., wnalphla}).
#' @param model character ("NB" or "PIG") indicating the type of model to be fitted, with
#' "NB" standing for Negative-Binomial and "PIG" standing for Poisson Inverse Gaussian. The default is "NB".
#' @param y For \code{mixpoissonregML}: logical values indicating if the response vector should be returned as component.
#'
#' For \code{mixpoissonregML.fit}: a numerical vector of response variables with length \code{n}. Each coordinate must be a nonnegative-integer.
#' @param x For \code{mixpoissonregML}: logical values indicating if the model matrix \code{x} should be returned as component.
#'
#' For \code{mixpoissonregML.fit}: a matrix of covariates with respect to the mean with dimension \code{(n,nbeta)}.
#' @param w For \code{mixpoissonregML}: logical values indicating if the model matrix \code{w} should be returned as component.
#'
#' For \code{mixpoissonregML.fit} a matrix of covariates with respect to the precision parameter. The default is \code{NULL}. If not \code{NULL} must be of dimension \code{(n,nalpha)}.
#' @param residual character indicating the type of residual to be evaluated ("pearson" or "score"). The default is "pearson". Notice that they coincide for Negative-Binomial models.
#' @param envelope number of simulations (synthetic data sets) to build envelopes for residuals (with \code{100*prob\%} confidence level).
#' The default \code{envelope = 0} dismisses the envelope analysis.
#' @param prob probability indicating the confidence level for the envelopes (default: \code{prob} = 0.95).
#' If \code{envelope} = 0, \code{prob} is ignored.
#' @param model.frame logical indicating whether the model frame should be returned as component of the returned value.
#' @param link.mean optionally, a string containing the link function for the mean. If omitted, the 'log' link function will be used.
#' The possible link functions for the mean are "log" and "sqrt".
#' @param link.precision optionally, a string containing the link function the precision parameter. If omitted and the only precision
#' covariate is the intercept, the 'identity' link function will be used, if omitted and there is a precision covariate other than the
#' intercept, the 'log' link function will be used. The possible link functions for the precision parameter are "identity" and "inverse.sqrt" (which is \eqn{\phi^{-1/2} = w_i^T alpha}).
#' @param em_controls only used with the 'EM' method. A list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000;
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5). \code{em_tolgrad} that defines the tolerance value
#' of the maximum-norm of the the gradient of the Q-function, the default is set to 10^(-2).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return \code{mixpoissonregML} returns an object of class "mixpoissonreg" whereas \code{mixpoissonregML.fit}
#' returns an object of class "mixpoissonreg_fit". Both objects are given by lists containing the outputs from the model fit (Negative-Binomial or Poisson Inverse Gaussian regression).
#'
#' An object of the class "mixpoissonreg" is a list containing the following elements:
#' \itemize{
#'   \item \code{coefficients} - a list with elements "mean" and "precision" containing the estimated coefficients of the model;
#'   \item \code{call} - the formula used by the model. If using \code{mixpoissonreg.fit}, this returns \code{NULL}.
#'   \item \code{modelname} - the fitted model, NB or PIG;
#'   \item \code{modeltype} - the abbreviated model name
#'   \item \code{residualname} - the name of the chosen residual in the call, 'pearson' or 'score';
#'   \item \code{niter} - number of iterations of the EM algorithm if method = "EM" and number of iterations
#'   of the \code{optim} function, if method = "ML";
#'   \item \code{start} - the initial guesses of the parameters
#'   \item \code{intercept} - vector indicating if the intercept is present in the mean and/or in the precision regressions;
#'   \item \code{link.mean} - link function of the mean;
#'   \item \code{link.precision} - link function of the precision parameter;
#'   \item \code{fitted.values} - a vector of fitted values in the response scale;
#'   \item \code{fitted.precisions} - a vector of fitted precisions;
#'   \item \code{efron.pseudo.r2} - Efron's pseudo R^2: the squared correlation between the response variables and the predicted values;
#'   \item \code{vcov} - covariance matrix of the parameters of the fitted model;
#'   \item \code{logLik} - log-likelihood at the estimated parameters;
#'   \item \code{Qfunction} - Q-function at the estimated parameters;
#'   \item \code{x} - the covariates related to the mean (if x = TRUE);
#'   \item \code{w} - the covariates related to the precision parameter (if w = TRUE);
#'   \item \code{y} - the response variables (if y = TRUE);
#'   \item \code{model} - if requested (the default), the model frame;
#'   \item \code{formula} - the formula supplied;
#'   \item \code{nobs} - number of observations
#'   \item \code{df.null} - the residual degrees of freedom for the model with constant mean and constant precision;
#'   \item \code{df.residual} - the residual degrees of freedom of the fitted model;
#'   \item \code{estimation_method} - the estimation method, "EM" or "ML"
#'   \item \code{residuals} - vector of raw residuals, that is, the response variable minus the fitted means;
#'   \item \code{std_errors} - the standard errors of the estimated parameters;
#'   \item \code{envelope} - the numerical envelopes used to build the Q-Q plot with simulated envelopes;
#'   \item \code{terms} - (only for \code{mixpoissonreg})the \code{terms} object used;
#'   \item \code{levels} - (where relevant, only for \code{mixpoissonreg}) the levels of the factors used;
#'   \item \code{contrasts} - (where relevant, only for \code{mixpoissonreg}) the contrasts used.

#' }
#'
#' @details Among the regression models with discrete response variables, Poisson regression is the most popular
#' for modeling count data. See, for instance Sellers and Shmueli (2010).
#' It is well-known that this model is equidispersed (that is, the mean is equal to the variance),
#' which in practice may be an unrealistic
#' assumption. Several models have been introduced in the literature to overcome this problem such as
#' negative binomial (NB) and Poisson inverse gaussian (PIG) distributions (see Lawless, 1987).
#' The most common way to do this is to consider a mixed Poisson distribution, which is defined as follows.
#' Let \eqn{Z} be a positive random variable (generally being continuous) with distribution
#' function \eqn{G_{\tau}(\cdot)},
#' where \eqn{\tau} denotes the parameter vector associated to the \eqn{G} distribution. Let
#' \eqn{Y|Z=z\sim}Poisson\eqn{(\mu z)}, for
#' some constant \eqn{\mu>0}. Therefore \eqn{Y} follows a mixed Poisson (MP) distribution with probability
#' function given by
#' \deqn{P(Y=y)=\int_0^\infty\frac{e^{-\mu z}(\mu z)^y}{y!}dG_{\tau}(z),}
#' for \eqn{y=0,1,\ldots}. With this,
#' \eqn{Y} has an overdispersed distribution and hence it is a natural alternative to the Poisson distribution.
#' The most common choices for \eqn{Z} are gamma and inverse-gaussian distributions,
#' which yields \eqn{Y} following, respectively, NB and PIG distributions.
#' General properties of the MP distributions can be found in Karlis and Xekalaki (2005) and in the references therein.
#'
#' In \code{mixpoissonreg} two regression models are implemented, namely, the NB and PIG regression models.
#' We follow the definitions and notations given in Barreto-Souza and Simas (2016). The mixed Poisson regression model
#' is defined by assuming \eqn{Y_1,\ldots,Y_n} is a random sample where
#' \eqn{Y_i\sim NB(\mu_i,\phi_i)} or \eqn{Y_i\sim PIG(\mu_i,\phi_i)} for \eqn{i = 1,\ldots,n}.
#' Under this parameterization we have \eqn{E(Y_i) = \mu_i} and \eqn{Var(Y_i) = \mu_i(1+\mu_i\phi_i^{-1}b''(\xi_0))}, where
#' \eqn{b(\theta) = -\log(-\theta)} and \eqn{\xi_0 = -1} for the NB case, and \eqn{b(\theta) = -(-2\theta)^{1/2}} and \eqn{\xi_0 = -1/2} for
#' the PIG case, with \eqn{b''(\cdot)} being the second derivative of the function \eqn{b(\cdot)}.
#' The following linear relations are assumed
#' \deqn{\Lambda_1(\mu_i) = x_i^T \beta}
#'  and
#' \deqn{\Lambda_2(\phi_i) = w_i^T \alpha,}
#'  where \eqn{\beta = (\beta_1,...,\beta_p)} and \eqn{\alpha = (\alpha_1,...,\alpha_q)} are real valued vectors.
#' The terms \eqn{x_i^T} and \eqn{v_i^T} represent, respectively, the i-th row of the matrices "x" (\eqn{n\times p})
#'  and "w" (\eqn{n\times q}) containing covariates in their columns
#' (\eqn{x_{i,1}} and \eqn{v_{i,1}} may be 1 to handle intercepts).
#'
#' Therefore, the \code{mixpoissonreg} package handles up to two regression structures
#' at the same time: one for the mean parameter, one for the precision parameter. The regression structure for
#' the mean is determined through a formula \code{y ~ x1 + ... + xn}, whereas the regression structure for
#' the precision parameter is determined through the right-hand side of the formula using the separator "\code{|}". So,
#' for example, a regression with \code{x1,...,xn} as covariates for the mean and \code{z1,...,zm} as covariates for the precision
#' parameter corresponds to the formula \code{y ~ x1 + ... + xn | z1 + ... + zm}. If only there is only formula for
#' the regression structure for the mean, the regression structure for the precision parameter will only have the intercept,
#' that is, \code{y ~ x1 + ... + xn} is the same as \code{y ~ x1 + ... + xn | 1}.
#'
#' In general, in this package, the EM-algorithm estimation method obtains estimates closer to the maximum likelihood estimate than the maximum likelihood estimation method,
#' in the sense that the likelihood function evaluated at the EM-algorithm estimate is greater or equal (usually strictly greater) than the likelihood function evaluated
#' at the maximum likelihood estimate. So, unless the processing time is an issue, we strongly recommend the EM-algorithm as the estimation method.
#'
#' In Barreto-Souza and Simas (2016) two residuals were studied: the pearson residuals
#' and the score residuals. Both these residuals are implemented in the \code{mixpoissonreg}
#' package. They coincide for NB regression models. They can be accessed via
#' the \link[mixpoissonreg:residuals.mixpoissonreg]{residuals} method.
#'
#' It is also noteworthy that all the global and local influence analysis tools developed
#' in Barreto-Souza and Simas (2016) are implemented in this package. See \code{\link{influence.mixpoissonreg}},
#' \code{\link{local_influence.mixpoissonreg}}, \code{\link{local_influence_plot.mixpoissonreg}}
#' and \code{\link{local_influence_autoplot.mixpoissonreg}}.
#'
#' @references
#' DOI:10.1007/s11222-015-9601-6 \doi{10.1007/s11222-015-9601-6}(Barreto-Souza and Simas; 2016)
#'
#' URL:https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1751-5823.2005.tb00250.x (\href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1751-5823.2005.tb00250.x}{Karlis and Xekalaki; 2005})
#'
#' DOI:10.2307/3314912 \doi{10.2307/3314912}(Lawless; 1987)
#'
#' Sellers, K.F. and Shmueli, G. (2010) *A flexible regression model for count data.* Ann. Appl. Stat., 4, 943-961
#'
#' @seealso
#' \code{\link{summary.mixpoissonreg}}, \code{\link{plot.mixpoissonreg}}, \code{\link{autoplot.mixpoissonreg}},
#' \code{\link{residuals.mixpoissonreg}}, \code{\link{predict.mixpoissonreg}},\code{\link{influence.mixpoissonreg}},
#' \code{\link{cooks.distance.mixpoissonreg}},
#' \code{\link{local_influence.mixpoissonreg}}, \code{\link{local_influence_plot.mixpoissonreg}}, \code{\link{local_influence_autoplot.mixpoissonreg}}
#'
#' @examples
#' # Examples using the Attendance dataset:
#' \donttest{
#' daysabs_fit_ml <- mixpoissonregML(daysabs ~ gender + math +
#' prog | gender + math + prog, data = Attendance)
#' summary(daysabs_fit_ml)
#' # Base R plot of the fit
#' plot(daysabs_fit_ml)
#' # ggplot2 plot of the fit
#' autoplot(daysabs_fit_ml)
#' # plot of local influence measures
#' local_influence_plot(daysabs_fit_ml)
#' # ggplot2 plot of local influence measures
#' local_influence_autoplot(daysabs_fit_ml)
#' # Fitting a reduced model of the sabe type as the previous one
#' daysabs_fit_ml_red <- mixpoissonregML(daysabs ~ gender + math +
#' prog | prog, data = Attendance, model = daysabs_fit_ml$modeltype)
#' # Likelihood ratio test:
#' lmtest::lrtest(daysabs_fit_ml, daysabs_fit_ml_red)
#' # Wald test:
#' lmtest::waldtest(daysabs_fit_ml, daysabs_fit_ml_red)
#' }
#'
#' @rdname mixpoissonregML
#' @export
mixpoissonregML <- function(formula, data, link.mean = c("log", "sqrt"),
                          link.precision = c("identity", "log", "inverse.sqrt"),
                          model = c("NB", "PIG"),
                          residual = c("pearson", "score"), y = TRUE, x = TRUE, w = TRUE,
                          envelope = 0, prob = 0.95, model.frame = TRUE, em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                          optim_method = "L-BFGS-B", optim_controls = list()) {
  call_mixpoissonreg <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }

  object <- mixpoissonreg(formula = formula, data = data, link.mean = link.mean,
                          link.precision = link.precision,
                          model = model, method = "ML",
                          residual = residual, y = y, x = x, w = w,
                          envelope = envelope, prob = prob, model.frame = model.frame, em_controls = em_controls,
                          optim_method = optim_method, optim_controls = optim_controls)

  object$call <- call_mixpoissonreg
  return(object)

}


#' @rdname mixpoissonregML
#' @export

mixpoissonregML.fit <- function(x, y, w = NULL, link.mean = c("log", "sqrt"),
                              link.precision = c("identity", "log", "inverse.sqrt"),
                              model = c("NB", "PIG"),
                              residual = c("pearson", "score"), envelope = 0,
                              prob = 0.95, em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)), optim_method = "L-BFGS-B", optim_controls = list()) {

  object <- mixpoissonreg.fit(y=y, x=x, w = w, link.mean = link.mean,
             link.precision = link.precision,
             model = model, method = "ML",
             residual = residual, envelope = envelope,
             prob = prob, em_controls = em_controls, optim_method = optim_method, optim_controls = optim_controls)

  return(object)
}
