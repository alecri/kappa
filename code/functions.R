## Collection of usefull functions


## Number of citation of an article by year
# taken and mdified (fixed bugs) by scholar::get_article_cite_history
get_article_cite_history <- function(id, article){
  library(magrittr)
  id <- id[1]
  url_base <- paste0("http://scholar.google.com/citations?", 
                     "view_op=view_citation&hl=en&citation_for_view=")
  url_tail <- paste(id, article, sep = ":")
  url <- paste0(url_base, url_tail)
  doc <- httr:::GET(url, handle = getOption("scholar_handle")) %>% 
    xml2::read_html()
  years <- doc %>% rvest::html_nodes(xpath = "//*/span[@class='gsc_vcd_g_t']") %>% 
    rvest::html_text() %>% as.numeric()
  vals <- doc %>% rvest::html_nodes(xpath = "//*/span[@class='gsc_vcd_g_al']") %>% 
    rvest::html_text() %>% as.numeric()
  df <- data.frame(year = years, cites = vals)
  tmp <- merge(data.frame(year = min(years):max(years)), df, 
               all.x = TRUE)
  tmp[is.na(tmp)] <- 0
  tmp$pubid <- article
  return(tmp)
}
