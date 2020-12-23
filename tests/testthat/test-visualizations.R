fit_ml1 <- mixpoissonregML(daysabs ~ math, data = Attendance)

plot(fit_ml1)

plot(fit_ml1, which = 2)

autoplot(fit_ml1)

autoplot(fit_ml1, nrow = 2)

autoplot(fit_ml1, which = 2)

expect_warning(local_influence_plot(fit_ml1))

local_influence_plot(fit_ml1, which = 2)

expect_warning(local_influence_autoplot(fit_ml1))

expect_warning(local_influence_autoplot(fit_ml1, which = 2))

expect_warning(local_influence_autoplot(fit_ml1, nrow = 2))
