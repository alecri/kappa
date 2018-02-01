library(dosresmeta)
library(tidyverse)
library(rms)
library(gridExtra)
library(cowplot)
library(knitr)
library(kableExtra)

data("coffee_cvd")

# linear analysis
lin <- dosresmeta(logrr ~ dose, id = id, se = se, type = type,
                  cases = cases, n = n, data = coffee_cvd, method = "fixed")
summary(lin)
gof(lin)

pglin <- gof(lin)$tdata %>%
  mutate(dose = coffee_cvd$dose[coffee_cvd$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "Linear") 
#pglin

# quadratic analysios
quadr <- dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, type = type,
                    cases = cases, n = n, data = coffee_cvd, method = "fixed")
summary(quadr)
gof(quadr)

pgquadr <- gof(quadr)$tdata %>%
  mutate(dose = coffee_cvd$dose[coffee_cvd$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "Quadratic") 
#pgquadr


# rcs spline
k <- quantile(coffee_cvd$dose, c(.1, .5, .9))
spl <- dosresmeta(logrr ~ rcs(dose, k), id = id, se = se, type = type,
                  cases = cases, n = n, data = coffee_cvd, method = "fixed")
summary(spl)
gof(spl)

pgspl <- gof(spl)$tdata %>%
  mutate(dose = coffee_cvd$dose[coffee_cvd$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "RCS") 
#pgspl


# meta-regression
spl_reg <- dosresmeta(logrr ~ rcs(dose, k), id = id, se = se, type = type,
                      cases = cases, n = n, data = coffee_cvd, method = "fixed",
                      mod = ~ area + smoking)
summary(spl_reg)
gof(spl_reg)

pgspl_reg <- gof(spl_reg)$tdata %>%
  mutate(dose = coffee_cvd$dose[coffee_cvd$se != 0]) %>%
  ggplot(aes(dose, tresiduals)) +
  geom_point() + geom_smooth(se = F, color = "black") +
  ylim(c(-4, 4)) + labs(x = "Coffee consumption (cups/day)",
                        y = "Decorellated residuals", title = "RCS + interaction") 
#pgspl_reg


# grid of pics
plot_grid(pglin, pgquadr, pgspl, pgspl_reg,
          align = 'vh', labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)


# table
tab_gof <- data.frame(
  analysis = c("A", "B", "C", "D"),
  model = c("Linear", "Quadratic", 
            paste("RCS", footnote_marker_alphabet(1, "latex")),
            paste("RCS + interaction", footnote_marker_alphabet(2, "latex")))
) %>%
  cbind(
    map(list(lin, quadr, spl, spl_reg), ~ gof(.x)) %>%
      map(~ c(deviance = .x$deviance$D, df = .x$deviance$df, p = .x$deviance$p,
              R2 = .x$R2, R2adj = .x$R2adj)) %>%
      do.call("rbind", .)
  )

# customized print for latex
kable(tab_gof, "latex", align = "c", booktabs = T, escape = F, digits = 3,
      caption = "prova", 
      col.names = c("Analysis", "Model", "Deviance", "df", "$p$~value", "\\textrm{$\\mathrm{R^2}$}", "\\textrm{$\\mathrm{R_{\\textrm{adj}}^2}$}")) %>%
  footnote(alphabet = c("3 knots located at the 10th, 50th, and 90th percentiles of the distribution of coffee.", 
                        "As in c) + interaction with gender of participants (only men, only women, both sexes) and geographical area (Europe, USA, Japan) included as categorical study-level covariates."),
           footnote_as_chunk = F, threeparttable = T)
