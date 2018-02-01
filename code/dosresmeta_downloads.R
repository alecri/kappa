# Number of downloads of the dosresmeta with info of country

library(tidyverse)
library(lubridate)
library(rworldmap)

num_dwnld <- read.table("data/pkgs_dwnld.txt", stringsAsFactors = FALSE) %>%
  mutate(
    date = as.Date(date),
    month = floor_date(date, "month"),
    week = floor_date(date, "week")
  ) %>%
  filter(package == "dosresmeta", month < "2017-11-01")

dosresmeta_count <- num_dwnld %>%
  group_by(month) %>%
  count()

dosresmeta_map <- num_dwnld %>%
  group_by(country, package) %>%
  count() %>% group_by(country) %>%
  summarise(n = sum(n)) %>%
  arrange(desc(n))
save(dosresmeta_count, dosresmeta_map, file = "data/count_dosresmeta.RData")

# time serie
ggplot(dosresmeta_count, aes(month, n)) + 
  geom_point(size = 2) + geom_line() +
  scale_x_date(date_breaks = "6 month", date_minor_breaks = "6 month", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Monthly downloads of dosresmeta", x = "")

# world map
map <- joinCountryData2Map(dosresmeta_map, joinCode = "ISO2", nameJoinColumn = "country")
old_par <- par()
par(mai = c(0, 0, 0.2, 0), xaxs = "i", yaxs = "i")
op <- palette(heat.colors(5)[5:1])
cutVector <- c(1, 50, 100, 500, 1500, 4100)
#classify the data to a factor
map@data[["n_cat"]] <- cut(map@data[["n"]], cutVector, include.lowest = TRUE, right = F)
#rename the categories
levels(map@data[["n_cat"]]) <- c("1-50", "50-100", "100-500", "500-1500", ">1500")
mapParams <- mapCountryData(map, nameColumnToPlot = "n_cat", catMethod = "categorical",
                            mapTitle = "", addLegend = FALSE, 
                            colourPalette = "palette", missingCountryCol = "white")
do.call(addMapLegendBoxes, c(mapParams, x = 'left', horiz = F, bg = NA, bty = "n",
                             title = ""))
suppressWarnings(par(old_par))

# alternative world map
map <- joinCountryData2Map(dosresmeta_map, joinCode = "ISO2", nameJoinColumn = "country")
old_par <- par()
par(mai = c(0, 0, 0.2, 0), xaxs = "i", yaxs = "i")
mapParams <- mapCountryData(map, nameColumnToPlot = "n", catMethod = "pretty",
                            mapTitle = "", numCats = 9, addLegend = FALSE, colourPalette = "heat")
do.call(addMapLegend, c(mapParams, legendWidth = 0.5, legendMar = 3, 
                        legendShrink = .5, legendLabels = "all"))
suppressWarnings(par(old_par))



