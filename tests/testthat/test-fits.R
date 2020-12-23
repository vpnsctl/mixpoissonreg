fit1 <- mixpoissonreg(daysabs ~ math, data = Attendance)
fit2 <- mixpoissonreg(daysabs ~ math | math, data = Attendance,
                      em_controls = list(maxit = 50))

fit1env <- mixpoissonreg(daysabs ~ math, data = Attendance, envelope = 10,
                         em_controls = list(maxit = 50))
fit2env <- mixpoissonreg(daysabs ~ math, data = Attendance, envelope = 10,
                         em_controls = list(maxit = 50))

fit1pig <- mixpoissonreg(daysabs ~ math, data = Attendance, model = "PIG")
fit2pig <- mixpoissonreg(daysabs ~ math | math, data = Attendance, model = "PIG",
                         em_controls = list(maxit = 50))

fit1pigenv <- mixpoissonreg(daysabs ~ math, data = Attendance, model = "PIG", envelope = 10,
                            em_controls = list(maxit = 50))
fit2pigenv <- mixpoissonreg(daysabs ~ math | math, data = Attendance, model = "PIG", envelope = 10,
                            em_controls = list(maxit = 50))

fit_ml1 <- mixpoissonregML(daysabs ~ math, data = Attendance)
fit_ml2 <- mixpoissonregML(daysabs ~ math | math, data = Attendance)

fit_ml1env <- mixpoissonregML(daysabs ~ math, data = Attendance, envelope = 10)
fit_ml2env <- mixpoissonregML(daysabs ~ math | math, data = Attendance, envelope = 10)

fit_ml1pig <- expect_warning(mixpoissonregML(daysabs ~ math, data = Attendance, model = "PIG"))
fit_ml2pig <- expect_warning(mixpoissonregML(daysabs ~ math | math, data = Attendance, model = "PIG"))

fit_ml1pigenv <- expect_warning(mixpoissonregML(daysabs ~ math, data = Attendance, model = "PIG", envelope = 10))
fit_ml2pigenv <- expect_warning(mixpoissonregML(daysabs ~ math | math, data = Attendance, model = "PIG", envelope = 10))

fit1sq <- mixpoissonreg(daysabs ~ math, data = Attendance, link.mean = "sqrt")
fit2sq <- expect_warning(mixpoissonreg(daysabs ~ math | math, data = Attendance, link.precision = "inverse.sqrt"))

fit_ml1sq <- mixpoissonregML(daysabs ~ math, data = Attendance, link.mean = "sqrt")
fit_ml2sq <- expect_warning(mixpoissonregML(daysabs ~ math | math, data = Attendance, link.precision = "inverse.sqrt"))

fit1pigsq <- mixpoissonreg(daysabs ~ math, data = Attendance, model = "PIG", link.mean = "sqrt")
fit2pigsq <- expect_warning(mixpoissonreg(daysabs ~ math | math, data = Attendance,
                                          model = "PIG", link.precision = "inverse.sqrt",
                                          em_controls = list(maxit = 50)))

fit_ml1pigsq <- mixpoissonregML(daysabs ~ math, data = Attendance, model = "PIG", link.mean = "sqrt")
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

print(fit_ml1)
summary(fit_ml1)
summary(fit_ml1env)

