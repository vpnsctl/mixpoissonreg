set.seed(33333)

fit1 <- mixpoissonreg(daysabs ~ prog + math, data = Attendance,
                                     em_controls = list(maxit = 5))

augment(fit1)

augment(fit1, newdata = data.frame(math = 1, prog = factor(c("Academic", "Vocational"), levels = c("General", "Academic", "Vocational"))))

expect_warning(augment(fit1, pred_int = TRUE, nsim_pred = 5, nsim_y = 5))

augment(fit1, type.predict = "link", se_fit = TRUE)

glance(fit1)

tidy(fit1)

tidy(fit1, conf.int = TRUE)