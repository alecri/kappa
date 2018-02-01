## SEER data

library(dosresmeta)
library(tidyverse)
library(xtable)
library(Epi)

# data selected from paper IV
load("data/breast_1501.RData")

# linear model on IPD
lin_IPD <- glm(Y ~ age + histology + maritalStatus + race +
                 tumor_size + radiation + grade + ER + PR + node_pos, 
               data = breast_1501, family = "poisson", offset = log(surv_time))
summary(lin_IPD)

# categorical model on IPD
breast_1501$age_cat <- cut(breast_1501$age, c(20, 49, 60, 70, 97),
                           include.lowest = T)
table(breast_1501$age_cat)
pos <- 1
cat_1501 <- glm(Y ~ relevel(age_cat, 1) + histology + maritalStatus + race +
                  tumor_size + radiation + grade + ER + PR + node_pos, 
                data = breast_1501, family = "poisson", offset = log(surv_time))
summary(cat_1501)

# creatig aggregated data
cat_res <- ci.exp(cat_1501)[1:4,] %>% 
  data.frame() %>%
  mutate(
    ci = paste0("(", round(X2.5., 2), ", ", round(X97.5., 2), ")"),
    logrr = coef(cat_1501)[1:4],
    se = diag(vcov(cat_1501))[1:4]^.5,
    refpos = c(pos, setdiff(1:4, pos))
  ) %>%
  arrange(refpos)
cat_res[pos, ] <- c(1, 1, 1, NA, 0, 0, pos)
tab_1501 <-breast_1501 %>%
  group_by(age_cat) %>%
  summarise(dose = mean(age),
            cases = sum(Y),
            PT = sum(surv_time)) %>%
  cbind(cat_res)

# linear model on AD
lin_ad <- dosresmeta(logrr ~ dose, type = "ir", se = se, cases = cases, n = PT,
                     data = tab_1501)
predict(lin_ad, delta = 1)
summary(lin_ad)
# lin_ad_indep <- dosresmeta(logrr ~ dose, type = "ir", se = se, cases = cases, n = PT,
#                            data = tab_1501, covariance = "indep")
# summary(lin_ad_indep)


# # print for latex
# tab_1501 %>%
#   select(-`X2.5.`, -`X97.5.`) %>%
#   `colnames<-`(c("Age category", "Mean age", "Cases",
#                  "Person-years", "IRR", "95% CI", "log IRR", "SE")) %>%
#   xtable(digits = 2, align = "ccccccccc", label = "tab:breast_ad",
#          caption = "Aggregated dose-response data on the adjusted association between age and breast cancer mortality based on one registry from the SEER program.
#          ") %>%
#   print(include.rownames = FALSE, caption.placement = "top")


ggplot(tab_1501, aes(x = dose, y = exp.Est..)) + 
  geom_point() +
  geom_errorbar(aes(ymin = X2.5., ymax = X97.5.)) +
  geom_line(aes(y = predict(lin_ad, xref_pos = pos, expo = T)$pred))
