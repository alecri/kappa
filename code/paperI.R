library(dosresmeta)
library(tidyverse)
library(scales)
library(knitr)
library(metafor)
library(hetmeta)
library(kableExtra)

data("coffee_mort")
group_by(coffee_mort, id) %>%
  summarize(min = min(dose)) %>% c()


# single study
legrady <- subset(coffee_mort, id == 1)
covar.logrr(cases = cases, n = n, y = logrr, v = I(se^2), type = type,
            data = legrady)

lin_le <- dosresmeta(logrr ~ dose, se = se, type = type, 
                     cases = cases, n = n, data = legrady)
summary(lin_le)
predict(lin_le, delta = 3, expo = TRUE)

quadr_le <- dosresmeta(logrr ~ dose + I(dose^2), se = se, type = type,
                       cases = cases, n = n, data = legrady)
summary(quadr_le)
predict(quadr_le, newdata = data.frame(dose = 0:6), expo = TRUE)


# multiple studies

lin_i <- coffee_mort %>%
  split(.$id) %>%
  map(~ dosresmeta(logrr ~ dose, se = se, type = type,
                   cases = cases, n = n, data = .x))
lin_bi <- map_dbl(lin_i, ~ coef(.x))
lin_vi <- map_dbl(lin_i, ~ vcov(.x))
head(cbind(bi = lin_bi, vi = lin_vi))
mvmeta(lin_bi, lin_vi, method = "mm")

lin <- dosresmeta(logrr ~ dose, id = id, se = se, type = type,
           cases = cases, n = n, data = coffee_mort, method = "mm")
summary(lin)
predict(lin, delta = 3, expo = TRUE)

quadr_i <- coffee_mort %>%
  split(.$id) %>%
  map(~ dosresmeta(logrr ~ dose + I(dose^2), se = se, type = type,
                   cases = cases, n = n, data = .x))
quadr_bi <- t(map_df(quadr_i, ~ coef(.x)))
quadr_vi <- map(quadr_i, ~ vcov(.x))

mvmeta(quadr_bi ~ 1, quadr_vi, method = "mm")

quadr <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, type = type,
                    cases = cases, n = n, data = coffee_mort)
summary(quadr)

xref <- 0
pred <- data.frame(dose = c(xref, seq(0, 8, .1))) %>%
  predict(quadr, newdata = ., expo = T) %>%
  cbind(lin = predict(lin, newdata = ., expo = T))
ggplot(pred, aes(dose, pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line() + geom_ribbon(alpha = .1) +
  geom_line(aes(y = lin.pred), linetype = "dashed") +
  scale_y_continuous(trans = "log", breaks = scales::pretty_breaks())


filter(pred, dose %in% 0:9) %>%
  select(-lin.dose) %>% unique() %>% 
  knitr::kable(digits = 3, format = "latex", booktabs = T, row.names = FALSE) %>%
  kable_styling(position = "center")


#lin_rma <- rma.uni(yi = lin_bi, vi = lin_vi)
lin_rma <- rma.uni(yi = lin$bi, vi = unlist(lin$Si), method = "DL")
hetmeta(lin_rma)



newd <- data.frame(dose = c(xref, seq(0, 6, .1)))
newd %>%
  cbind(map(lin$bi, ~ exp(.x*newd$dose - 2))) %>%
  #purrr::set_names(c("dose", paste("study", 1:nrow(lin$bi)))) %>%
  gather(study, pred, -dose) %>%
  ggplot(aes(dose, pred, group = study)) + geom_line() +
  scale_y_continuous(trans = "log", breaks = scales::pretty_breaks()) +
  labs(x = "Coffee consumption (cups/day)", y = "Relative Risk")
newd %>%
  cbind(map(array_branch(quadr$bi, 1), ~ exp(.x[1]*newd$dose + .x[2]*newd$dose^2))) %>%
  #purrr::set_names(c("dose", paste("study", 1:nrow(quadr$bi)))) %>%
  gather(study, pred, -dose) %>%
  ggplot(aes(dose, pred, group = study)) + geom_line() +
  scale_y_continuous(trans = "log", breaks = c(.5, 1, 2, 5, 10), limits = c(.25, 10)) +
  labs(x = "Coffee consumption (cups/day)", y = "Relative Risk")
