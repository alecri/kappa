library(dosresmeta)
library(tidyverse)
library(rms)
library(gridExtra)
library(cowplot)

data("coffee_mort")

# linear analysis (analysis A)
lin <- dosresmeta(logrr ~ dose, id = id, se = se, type = type,
                  cases = cases, n = n, data = coffee_mort, method = "fixed")
summary(lin)
gof(lin)

# decorrelated plot
pglin <- gof(lin)$tdata %>%
  mutate(dose = coffee_mort$dose[coffee_mort$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "Linear") 
pglin

# quadratic analysis (analysis B)
quadr <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, type = type,
                    cases = cases, n = n, data = coffee_mort, method = "fixed")
summary(quadr)
gof(quadr)

pgquadr <- gof(quadr)$tdata %>%
  mutate(dose = coffee_mort$dose[coffee_mort$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "Quadratic") 
pgquadr


# rcs spline (analysis C)
k <- quantile(coffee_mort$dose, c(.1, .5, .9))
spl <- dosresmeta(logrr ~ rcs(dose, k), id = id, se = se, type = type,
                  cases = cases, n = n, data = coffee_mort, method = "fixed")
summary(spl)
gof(spl)

pgspl <- gof(spl)$tdata %>%
  mutate(dose = coffee_mort$dose[coffee_mort$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "RCS") 
pgspl


# meta-regression (analysis D)
spl_reg <- dosresmeta(logrr ~ rcs(dose, k), id = id, se = se, type = type,
                  cases = cases, n = n, data = coffee_mort, method = "fixed",
                  mod = ~ gender + area)
summary(spl_reg)
gof(spl_reg)

pgspl_reg <- gof(spl_reg)$tdata %>%
  mutate(dose = coffee_mort$dose[coffee_mort$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "RCS + interaction") 
pgspl_reg


# grid of pics
plot_grid(pglin, pgquadr, pgspl, pgspl_reg,
  align = 'vh', labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)


# table
tab_gof <- data.frame(
  analysis = c("A", "B", "C", "D"),
  model = c("Linear", "Quadratic",  "RCS", "RCS + interaction")
  ) %>%
  cbind(
    map(list(lin, quadr, spl, spl_reg), ~ gof(.x)) %>%
      map(~ c(deviance = .x$deviance$D, df = .x$deviance$df, p = .x$deviance$p,
              R2 = .x$R2, R2adj = .x$R2adj)) %>%
      do.call("rbind", .)
  )
tab_gof
