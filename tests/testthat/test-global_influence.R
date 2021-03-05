set.seed(33333)

fit1 <- mixpoissonreg(daysabs ~ prog + math, data = Attendance,
                                     em_controls = list(maxit = 1))
fit2 <- mixpoissonreg(daysabs ~ prog + math | math, data = Attendance,
                                     em_controls = list(maxit = 1))

hatvalues(fit1)

hatvalues(fit1, parameters = "precision")

hatvalues(fit2)

hatvalues(fit2, parameters = "precision")

cooks.distance(fit1)

cooks.distance(fit1, hat = "precision")

cooks.distance(fit1, type = "GCD")

cooks.distance(fit1, type = "GCDmean")

cooks.distance(fit1, type = "GCDprecision")

cooks.distance(fit1, type = "LD")

cooks.distance(fit1, type = "QD")

cooks.distance(fit2)

cooks.distance(fit2, hat = "precision")

cooks.distance(fit2, type = "GCD")

cooks.distance(fit2, type = "GCDmean")

cooks.distance(fit2, type = "GCDprecision")

cooks.distance(fit2, type = "LD")

cooks.distance(fit2, type = "QD")

influence(fit1, do.coef = TRUE)

influence(fit2, do.coef = TRUE)

fit1_PIG <- expect_warning(mixpoissonreg(daysabs ~ prog + math, data = Attendance, model = "PIG",
                                     em_controls = list(maxit = 1)))
fit2_PIG <- expect_warning(mixpoissonreg(daysabs ~ prog + math | math, data = Attendance, model = "PIG",
                                     em_controls = list(maxit = 1)))

hatvalues(fit1_PIG)

hatvalues(fit1_PIG, parameters = "precision")

hatvalues(fit2_PIG)

hatvalues(fit2_PIG, parameters = "precision")

cooks.distance(fit1_PIG)

cooks.distance(fit1_PIG, hat = "precision")

cooks.distance(fit1_PIG, type = "GCD")

cooks.distance(fit1_PIG, type = "GCDmean")

cooks.distance(fit1_PIG, type = "GCDprecision")

cooks.distance(fit1_PIG, type = "LD")

cooks.distance(fit1_PIG, type = "QD")

cooks.distance(fit2_PIG)

cooks.distance(fit2_PIG, hat = "precision")

cooks.distance(fit2_PIG, type = "GCD")

cooks.distance(fit2_PIG, type = "GCDmean")

cooks.distance(fit2_PIG, type = "GCDprecision")

cooks.distance(fit2_PIG, type = "LD")

cooks.distance(fit2_PIG, type = "QD")

influence(fit1_PIG, do.coef = TRUE)

influence(fit2_PIG, do.coef = TRUE)
