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


# easily get prediction for rma.uni models
pred.rma <- function(object, transf, string = 0, digits = 2){
  pred <- if (missing(transf)){
    predict(object)
  } else {
    predict(object, transf = transf)
  }
  pred <- suppressWarnings(as.double(unlist(pred)))[c(1, 3, 4)]
  pred <- round(pred, digits)
  names(pred) <- c("pred", "ci.lb", "ci.ub")
  if (string == 1){
    pred <- paste0(pred[1], " (", pred[2], ", ", pred[3], ")")
  }
  pred
}


# label for forest plot
Qlab <- function(mod, total = TRUE){
  str0 <- if (total){
    #      "Overall (I-squared = "
    "Overall (Rb = "
  } else {
    #      "Subtotal (I-squared = "
    "Subtotal (Rb = "
  }
  #   I2str <- round(mod$I2, 0)
  Rbstr <- round(hetmeta(mod)$Rb, 0)
  str1 <- "%, p"
  pQstr <- if (mod$QEp <0.01){
    " < 0.01)"
  } else {
    paste0(" = ", round(mod$QEp, 2), ")")
  }
  #   paste0(str0, I2str, str1, pQstr)
  paste0(str0, Rbstr, str1, pQstr)
}
