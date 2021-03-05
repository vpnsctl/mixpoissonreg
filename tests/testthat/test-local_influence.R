set.seed(33333)

fit_ml1 <- mixpoissonregML(daysabs ~ prog + math, data = Attendance, 
                           optim_controls = list(maxit=1))

expect_warning(local_influence(fit_ml1))

expect_warning(local_influence(fit_ml1, mean.covariates = "prog"))

expect_warning(local_influence(fit_ml1, curvature = "normal"))

expect_warning(local_influence(fit_ml1, curvature = "normal", mean.covariates = "prog"))

expect_warning(local_influence(fit_ml1, direction = "max.eigen"))

expect_warning(local_influence(fit_ml1, direction = "max.eigen", mean.covariates = "prog"))

expect_warning(local_influence(fit_ml1, direction = "max.eigen", curvature = "normal"))

expect_warning(local_influence(fit_ml1, direction = "max.eigen", curvature = "normal", mean.covariates = "prog"))

fit_ml2 <- mixpoissonregML(daysabs ~ prog + math | prog, data = Attendance)

local_influence(fit_ml2)

local_influence(fit_ml2, precision.covariates = "prog")

local_influence(fit_ml2, precision.covariates = "prog", mean.covariates = "prog")

local_influence(fit_ml2, curvature = "normal")

local_influence(fit_ml2, curvature = "normal", precision.covariates = "prog")

local_influence(fit_ml2, direction = "max.eigen")

local_influence(fit_ml2, direction = "max.eigen", precision.covariates = "prog")

local_influence(fit_ml2, direction = "max.eigen", precision.covariates = "prog", mean.covariates = "prog")

local_influence(fit_ml2, direction = "max.eigen", curvature = "normal")

local_influence(fit_ml2, direction = "max.eigen", curvature = "normal", precision.covariates = "prog")

local_influence(fit_ml2, direction = "max.eigen", curvature = "normal", precision.covariates = "prog", mean.covariates = "prog")

local_influence(fit_ml2, curvature = "normal", precision.covariates = "prog", mean.covariates = "prog")

