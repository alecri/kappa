# Number of downloads of the dosresmeta and hetmeta packages
# https://cranlogs.r-pkg.org/

library(cranlogs)
library(tidyverse)
library(lubridate)

ndown <- cran_downloads(packages = c("dosresmeta", "hetmeta"), 
                        from = "2013-09-09", to = "2017-09-01") %>%
  mutate(time = floor_date(date, "month")) %>%
  group_by(package, time) %>%
  summarise(n = sum(count))

ndown %>%
  filter(n > 0) %>%
  ggplot(aes(time, n, group = package, col = package)) +
  geom_line() + geom_point(aes(shape = package)) +
  labs(x = "", y = "Monthly downloads from Rstudio CRAN mirror") +
  scale_x_date(date_breaks = "6 month", date_minor_breaks = "3 month",
               date_labels = "%b %Y") + ylim(c(0, 500)) + 
  theme_light() + 
  theme(panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1))
