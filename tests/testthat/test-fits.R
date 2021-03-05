set.seed(33333)

fit1 <- mixpoissonreg(daysabs ~ prog + math, data = Attendance,
                      em_controls = list(maxit = 1))
fit2 <- mixpoissonreg(daysabs ~ prog + math | math, data = Attendance,
                      em_controls = list(maxit = 1))

fit1env <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math, data = Attendance, envelope = 19,
                         em_controls = list(maxit = 1)))
fit2env <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math, data = Attendance, envelope = 10,
                         em_controls = list(maxit = 1)))

fit1pig <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG",
                         em_controls = list(maxit = 1)))
fit2pig <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math | math, data = Attendance, model = "PIG",
                         em_controls = list(maxit = 1)))

fit1pigenv <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG", envelope = 10,
                            em_controls = list(maxit = 1)))
fit2pigenv <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math | math, data = Attendance, model = "PIG", envelope = 10,
                            em_controls = list(maxit = 1)))

fit_ml1 <- mixpoissonregML(daysabs ~ prog + math, data = Attendance)

fit_ml2 <- mixpoissonregML(daysabs ~ prog + math | math, data = Attendance)

fit_ml1env <- mixpoissonregML(daysabs ~ prog + math, data = Attendance, envelope = 10,
                              optim_controls = list(maxit=1))
fit_ml2env <- mixpoissonregML(daysabs ~ prog + math | math, data = Attendance, envelope = 10,
                              optim_controls = list(maxit=1))

fit_ml1pig <- suppressWarnings(mixpoissonregML(daysabs ~ prog + math, data = Attendance, model = "PIG",
                               optim_controls = list(maxit=1)))
fit_ml2pig <- suppressWarnings(mixpoissonregML(daysabs ~ prog + math | math, data = Attendance, model = "PIG",
                               optim_controls = list(maxit=1)))

fit_ml1pigenv <- suppressWarnings(mixpoissonregML(daysabs ~ prog + math, data = Attendance, model = "PIG", envelope = 10,
                                  optim_controls = list(maxit=1)))
fit_ml2pigenv <- suppressWarnings(mixpoissonregML(daysabs ~ prog + math | math, data = Attendance, model = "PIG", envelope = 10,
                                  optim_controls = list(maxit=1)))

fit1sq <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math, data = Attendance, link.mean = "sqrt",
                                       em_controls = list(maxit = 1)))
fit2sq <- suppressWarnings(mixpoissonreg(daysabs ~ math | math, data = Attendance, link.precision = "inverse.sqrt",
                                       em_controls = list(maxit = 1)))

fit_ml1sq <- mixpoissonregML(daysabs ~ prog + math, data = Attendance, link.mean = "sqrt")
fit_ml2sq <- suppressWarnings(mixpoissonregML(daysabs ~ math | math, data = Attendance, link.precision = "inverse.sqrt"))

fit1pigsq <- suppressWarnings(mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG", link.mean = "sqrt",
                           em_controls = list(maxit = 1)))
fit2pigsq <- suppressWarnings(mixpoissonreg(daysabs ~ math | math, data = Attendance,
                                          model = "PIG", link.precision = "inverse.sqrt",
                                          em_controls = list(maxit = 1)))

fit_ml1pigsq <- mixpoissonregML(daysabs ~ prog + math, data = Attendance, model = "PIG", link.mean = "sqrt",
                                optim_controls = list(maxit=1))
fit_ml2pigsq <- suppressWarnings(mixpoissonregML(daysabs ~ math | math, data = Attendance, model = "PIG", link.precision = "inverse.sqrt",
                                                 optim_controls = list(maxit=1)))


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
predict(fit1, interval = "confidence", type = "link")
suppressWarnings(predict(fit1, interval = "confidence", type = "precision"))
suppressWarnings(predict(fit1, interval = "confidence", type = "variance"))

suppressWarnings(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2))
suppressWarnings(predict(fit1, interval = "prediction", type = "link", nsim_pred = 2, nsim_pred_y = 2))
suppressWarnings(predict(fit1, interval = "prediction", type = "precision", nsim_pred = 2, nsim_pred_y = 2))
suppressWarnings(predict(fit1, interval = "prediction", type = "variance", nsim_pred = 2, nsim_pred_y = 2))
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
suppressWarnings(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2, newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational")))))
suppressWarnings(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2, type = "link", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational")))))
suppressWarnings(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2, type = "precision", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational")))))
suppressWarnings(predict(fit1, interval = "prediction", nsim_pred = 2, nsim_pred_y = 2, type = "variance", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational")))))

predict(fit1, interval = "confidence", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))
predict(fit1, interval = "confidence", type = "link", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))
suppressWarnings(predict(fit1, interval = "confidence", type = "precision", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational")))))
suppressWarnings(predict(fit1, interval = "confidence",  type = "variance", newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational")))))


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

update(fit1, . ~ . - 1)

set.seed(3333)

x <- rexp(30)

y <- rNBI(30, mu = exp(2 - x))

mixpoissonregML(y ~ x)

expect_error(mixpoissonreg( x ~ -1))

expect_error(mixpoissonreg( y ~ x | - 1))

expect_error(mixpoissonreg( ~ 1))

expect_error(mixpoissonreg(daysabs ~ prog + math, data = Attendance, prob = -1))

expect_error(mixpoissonreg(I(daysabs-300) ~ prog + math, data = Attendance))

expect_error(mixpoissonreg(I(daysabs+0.5) ~ prog + math, data = Attendance))

expect_error(mixpoissonreg(daysabs ~ prog + math, data = Attendance, residual = 1))

expect_error(mixpoissonreg(daysabs ~ prog + math, data = Attendance, envelope = "bla"))

expect_error(mixpoissonreg(daysabs ~ prog + math, data = Attendance, envelope = -1))

expect_error(mixpoissonreg(daysabs ~ prog + math, data = Attendance, prob = "bla"))

suppressWarnings(mixpoissonreg(daysabs ~ math, data = Attendance, prob = 0.9))


expect_error(mixpoissonregML( x ~ -1))

expect_error(mixpoissonregML( y ~ x | - 1))

expect_error(mixpoissonregML( ~ 1))

expect_error(mixpoissonregML(daysabs ~ prog + math, data = Attendance, prob = -1))

expect_error(mixpoissonregML(I(daysabs-300) ~ prog + math, data = Attendance))

expect_error(mixpoissonregML(I(daysabs+0.5) ~ prog + math, data = Attendance))

expect_error(mixpoissonregML(daysabs ~ prog + math, data = Attendance, residual = 1))

expect_error(mixpoissonregML(daysabs ~ prog + math, data = Attendance, envelope = "bla"))

expect_error(mixpoissonregML(daysabs ~ prog + math, data = Attendance, envelope = -1))

expect_error(mixpoissonregML(daysabs ~ prog + math, data = Attendance, prob = "bla"))

suppressWarnings(mixpoissonregML(daysabs ~ math, data = Attendance, prob = 0.9))

x1 <- rexp(30)

x2 <- rnorm(30)

y <- rpois(30, exp(2+2*x1 - 2*x2))

suppressWarnings(mixpoissonreg(y ~ x1, envelope = 10, model = "PIG", em_controls = list(maxit = 20)))

suppressWarnings(mixpoissonreg(y ~ x1, envelope = 10, model = "NB", em_controls = list(maxit = 20)))

set.seed(1216)

x1 <- rexp(30)

x2 <- rnorm(30)

y <- rpois(30, exp(2))

expect_error(mixpoissonregML(y ~ x1+x2 | x2 + x1 -1, envelope = 500, model = "NB"))

expect_error(mixpoissonregML(y ~ x1 + x2 | x1 + x2, envelope = 500, model = "PIG"))

