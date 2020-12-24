set.seed(33333)

fit1 <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance,
                      em_controls = list(maxit = 5)))
fit2 <- expect_warning(mixpoissonreg(daysabs ~ prog + math | math, data = Attendance,
                      em_controls = list(maxit = 5)))

fit1env <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance, envelope = 10,
                         em_controls = list(maxit = 5)))
fit2env <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance, envelope = 10,
                         em_controls = list(maxit = 5)))

fit1pig <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG",
                         em_controls = list(maxit = 5)))
fit2pig <- expect_warning(mixpoissonreg(daysabs ~ prog + math | math, data = Attendance, model = "PIG",
                         em_controls = list(maxit = 5)))

fit1pigenv <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG", envelope = 10,
                            em_controls = list(maxit = 5)))
fit2pigenv <- expect_warning(mixpoissonreg(daysabs ~ prog + math | math, data = Attendance, model = "PIG", envelope = 10,
                            em_controls = list(maxit = 10)))

fit_ml1 <- mixpoissonregML(daysabs ~ prog + math, data = Attendance)

fit_ml2 <- mixpoissonregML(daysabs ~ prog + math | math, data = Attendance)

fit_ml1env <- expect_warning(mixpoissonregML(daysabs ~ prog + math, data = Attendance, envelope = 10))
fit_ml2env <- mixpoissonregML(daysabs ~ prog + math | math, data = Attendance, envelope = 10)

fit_ml1pig <- expect_warning(mixpoissonregML(daysabs ~ prog + math, data = Attendance, model = "PIG"))
fit_ml2pig <- expect_warning(mixpoissonregML(daysabs ~ prog + math | math, data = Attendance, model = "PIG"))

fit_ml1pigenv <- expect_warning(mixpoissonregML(daysabs ~ prog + math, data = Attendance, model = "PIG", envelope = 10))
fit_ml2pigenv <- expect_warning(mixpoissonregML(daysabs ~ prog + math | math, data = Attendance, model = "PIG", envelope = 10))

fit1sq <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance, link.mean = "sqrt",
                                       em_controls = list(maxit = 10)))
fit2sq <- expect_warning(mixpoissonreg(daysabs ~ math | math, data = Attendance, link.precision = "inverse.sqrt",
                                       em_controls = list(maxit = 10)))

fit_ml1sq <- mixpoissonregML(daysabs ~ prog + math, data = Attendance, link.mean = "sqrt")
fit_ml2sq <- expect_warning(mixpoissonregML(daysabs ~ math | math, data = Attendance, link.precision = "inverse.sqrt"))

fit1pigsq <- mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG", link.mean = "sqrt")
fit2pigsq <- expect_warning(mixpoissonreg(daysabs ~ math | math, data = Attendance,
                                          model = "PIG", link.precision = "inverse.sqrt",
                                          em_controls = list(maxit = 5)))

fit_ml1pigsq <- mixpoissonregML(daysabs ~ prog + math, data = Attendance, model = "PIG", link.mean = "sqrt")
fit_ml2pigsq <- expect_warning(mixpoissonregML(daysabs ~ math | math, data = Attendance, model = "PIG", link.precision = "inverse.sqrt"))


print(fit1)
summary(fit1)
coef(fit1)
coef(fit1, parameters = "mean")
coef(fit1, parameters = "precision")
vcov(fit1)
vcov(fit1, parameters = "mean")
vcov(fit1, parameters = "precision")
predict(fit1)
predict(fit1, interval = "confidence")
expect_warning(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2))
predict(fit1, type = "link", se.fit = TRUE)
predict(fit1, type = "precision")
predict(fit1, type = "variance")
residuals(fit1)
residuals(fit1, type = "score")
residuals(fit1pig)
residuals(fit1pig, type = "score")

terms(fit1)

coeftest.mixpoissonreg(fit1)

coefci(fit1)

predict(fit1, newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))
expect_warning(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2), newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))
predict(fit1, type = "link", se.fit = TRUE, newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))
predict(fit1, type = "precision", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))
predict(fit1, type = "variance", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))

logLik(fit1)

obs_fisher_weight_matrix_mixpoisson(fit2)

obs_fisher_weight_matrix_mixpoisson(fit2, parameters = "mean")

obs_fisher_weight_matrix_mixpoisson(fit2, parameters = "precision")

print(fit_ml1)
summary(fit_ml1)

print(summary(fit_ml1))

summary(fit_ml1env)

print(summary(fit_ml1env))

lmtest::lrtest(fit2, fit1)

lmtest::waldtest(fit2, fit1)

build_links_mpreg("inverse.sqrt")

d2mudeta2("log", 1)

d2mudeta2("sqrt", 1)

d2phideta2("inverse.sqrt", 1)

update(fit1, . ~ . - 1)


x <- rexp(30)

y <- rNBI(30, mu = exp(2 - x))

mixpoissonregML(y ~ x)

expect_error(mixpoissonreg( x ~ -1))

expect_error(mixpoissonreg( x ~ x | - 1))


