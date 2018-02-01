## motivating example for pointwise

library(dosresmeta)
library(tidyverse)
library(gridExtra)


data("milk_mort")

milk_mort %>%
  group_by(id) %>%
  summarise(min = min(dose), max = max(dose))

id_milk <- c(14, 11, 9, 3) 
milk_sub <- filter(milk_mort, id %in% id_milk)

milk_sub %>%
  group_by(id) %>%
  summarise(min = min(dose), max = max(dose))

fps <- fpgrid()
shift <- 5
scale <- 10

modi_milk <- lapply(split(milk_sub, milk_sub$id), function(d){
  lapply(split(fps, 1:nrow(fps)), function(p){
    dosresmeta(logrr ~ fracpol(dose, p = p, shift, scale),
               se = se, type = type, cases = cases, n = n, data = d)
  })
})

# fit_fp <- milk_mort %>%
#   group_by(id) %>%
#   nest()
# for (i in 1:nrow(fps)){
#   fit_fp[[paste0("fit", i)]] <- map(fit_fp$data, ~ dosresmeta(logrr ~ fracpol(dose, p = fps[i, ], shift, scale), 
#                                         se = se, type = type, 
#                                         cases = cases, n = n, data = .x))
# }


AIC_i <- lapply(modi_milk, function(m) sapply(m, AIC))
best_i <- sapply(AIC_i, which.min)
fps[best_i, ]
modi_milk_sel <- Map(function(m, i) m[[i]], modi_milk, best_i)


p_milk <- Map(function(d, m){
  newdata <- data.frame(dose = seq(min(d$dose), max(d$dose), length.out = 100)) %>%
    mutate(pred = predict(m, ., expo = T)$pred)
  ggplot(d, aes(dose, rr)) +
    geom_errorbar(aes(ymin = lb, ymax = ub)) +
    scale_y_continuous(trans = "log", breaks = c(.5, .6, .7, .8, .9, 1, 1.2, 1.5, 2, 2.5)) +
    labs(y = "Relative risk", x = "Milk consumption (ml/day)", title = paste("Study ID", d$id[1])) +
    geom_line(data = newdata, aes(x = dose, y = pred))
}, d = split(milk_sub, milk_sub$id), m = modi_milk_sel)

arrangeGrob(grobs = p_milk, ncol = 2, nrow = 2)
