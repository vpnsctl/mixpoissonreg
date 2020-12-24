set.seed(33333)

fit_ml1 <- mixpoissonregML(daysabs ~ prog + math, data = Attendance)

local_influence(fit_ml1)

local_influence(fit_ml1, curvature = "normal")

local_influence(fit_ml1, direction = "max.eigen")

local_influence(fit_ml1, direction = "max.eigen", curvature = "normal")

fit_ml2 <- mixpoissonregML(daysabs ~ prog + math | prog, data = Attendance)

local_influence(fit_ml2)

local_influence(fit_ml2, curvature = "normal")

local_influence(fit_ml2, direction = "max.eigen")

local_influence(fit_ml2, direction = "max.eigen", curvature = "normal")