set.seed(33333)

fit_ml1env <- expect_warning(mixpoissonregML(daysabs ~ prog + math, data = Attendance, envelope = 10))

fit_ml1 <- mixpoissonregML(daysabs ~ math, data = Attendance)

plot(fit_ml1, which = 1:6)

autoplot(fit_ml1, nrow = 2)


autoplot(fit_ml1, which = 1:6)

expect_warning(local_influence_plot(fit_ml1, which = 1:5))

expect_warning(local_influence_plot(fit_ml1, which = 1:5, direction = "max.eigen"))

local_influence_plot(fit_ml1, which = 2)

expect_warning(local_influence_autoplot(fit_ml1, which = 1:5))

expect_warning(local_influence_autoplot(fit_ml1, which = 1:5, direction = "max.eigen"))

expect_warning(local_influence_autoplot(fit_ml1, which = 2))

expect_warning(local_influence_autoplot(fit_ml1, nrow = 2))



plot(fit_ml1env)

plot(fit_ml1env, which = 2)

autoplot(fit_ml1env, which = 2)
