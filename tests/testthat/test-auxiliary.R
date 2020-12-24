set.seed(3333)

startvalues_mpreg(rpois(30, 10), as.matrix(rexp(30)), as.matrix(rep(1,30)), "log", "log", "NB")

build_links_mpreg("log")

build_links_mpreg("inverse.sqrt")

expect_error(build_links_mpreg("logit"))

d2mudeta2("log", 1)

d2mudeta2("sqrt", 1)

d2phideta2("identity", 1)

d2phideta2("log", 1)

d2phideta2("inverse.sqrt", 1)

generate_data_mixpoisson(list(mean = 1, precision = 1), as.matrix(rep(1,30)), as.matrix(rep(1,30)), 1, "log", "log", "NB")

generate_data_mixpoisson(list(mean = 1, precision = 1), as.matrix(rep(1,30)), as.matrix(rep(1,30)), 1, "log", "log", "PIG")

expect_warning(envelope_mixpoisson("pearson", "EM", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "NB", em_controls = list(maxit = 50, em_tol = 10^(-1), em_tolgrad = 10^(-1)),
                                    optim_method = "L-BFGS-B", optim_controls = list()))

expect_warning(envelope_mixpoisson("score", "EM", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "NB", em_controls = list(maxit = 50, em_tol = 10^(-1), em_tolgrad = 10^(-1)),
                    optim_method = "L-BFGS-B", optim_controls = list()))

envelope_mixpoisson("pearson", "ML", list(mean = c(exp(-1)), precision = 1),  as.matrix(rep(1,30)), as.matrix(rep(1,30)), 2, 0.95, 30, "sqrt", "identity", "NB", em_controls = list(maxit = 50, em_tol = 10^(-1), em_tolgrad = 10^(-1)),
                    optim_method = "Nelder-Mead", optim_controls = list())

envelope_mixpoisson("score", "ML", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "NB", em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                    optim_method = "L-BFGS-B", optim_controls = list())

envelope_mixpoisson("pearson", "EM", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "PIG", em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                    optim_method = "L-BFGS-B", optim_controls = list())

envelope_mixpoisson("score", "EM", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "PIG", em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                    optim_method = "L-BFGS-B", optim_controls = list())

envelope_mixpoisson("pearson", "ML", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "PIG", em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                    optim_method = "L-BFGS-B", optim_controls = list())

expect_warning(envelope_mixpoisson("score", "ML", list(mean = c(1), precision = 1),  as.matrix(rexp(30)), as.matrix(rep(1,30)), 2, 0.95, 30, "log", "log", "PIG", em_controls = list(maxit = 5000, em_tol = 10^(-5), em_tolgrad = 10^(-2)),
                    optim_method = "L-BFGS-B", optim_controls = list()))

pearson_residual_mixpoisson(list(mean = 1, precision = 1), rexp(30), as.matrix(rep(1,30)), as.matrix(rep(1,30)), "log", "log", "NB")

pearson_residual_mixpoisson(list(mean = 1, precision = 1), rexp(30), as.matrix(rep(1,30)), as.matrix(rep(1,30)), "log", "log", "PIG")

score_residual_mixpoisson(list(mean = 1, precision = 1), rexp(30), as.matrix(rep(1,30)), as.matrix(rep(1,30)), "log", "log", "NB")

score_residual_mixpoisson(list(mean = 1, precision = 1), rexp(30), as.matrix(rep(1,30)), as.matrix(rep(1,30)), "log", "log", "PIG")

expect_warning(std_error_mixpoisson(list(mean = 1, precision = 1), rexp(30), as.matrix(rep(1,30)), as.matrix(rep(1,30)), "log", "log", "NB"))

lambda_r(1,1,1,"NB")

lambda_r(1,1,1,"PIG")

kappa_r(1,1,1,"NB")

lambda_r(1,1,1,"PIG")

