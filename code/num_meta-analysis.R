## number of published meta-analysis in pubmed
## last updated: 171219

library(RISmed)
library(tidyverse)

year <- 1980:2017
count <- vector()
for (i in seq_along(year)){
  count[i] <- QueryCount(EUtilsSummary("meta-analysis", type = "esearch", db = "pubmed", 
                                       mindate = year[i], maxdate = year[i], retmax = 90000))
}
count_meta <- data.frame(year = year, count = count)
save(count_meta, file = "data/count_meta.RData")

count_meta %>%
  subset(year <= 2017) %>%
  ggplot(aes(year, count)) + 
  geom_point(size = 2) + geom_line() +
  scale_x_continuous(breaks = seq(1980, 2017, 5)) +
  labs(y = "Number of meta-analysis", x = "Year") +
  theme_light() + 
  theme(panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1))
