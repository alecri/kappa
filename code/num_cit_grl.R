## number of citations of Greenland and Longnecker (1992)
# https://scholar.google.com/citations?view_op=view_citation&hl=en&citation_for_view=HcvPl18AAAAJ:IjCSPb-OGe4C

library(tidyverse)

source("code/functions.R")
citations <- get_article_cite_history("HcvPl18AAAAJ", "IjCSPb-OGe4C")

citations %>%
  subset(year < 2017) %>%
  ggplot(aes(year, cites)) + 
  geom_point(size = 2) + geom_line() +
  #geom_smooth(se = F) +
  scale_x_continuous(breaks = seq(1992, 2017, 4)) +
  labs(y = "Number of citations", x = "Year") +
  theme_light() + 
  theme(panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1))