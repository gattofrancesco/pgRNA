getannotatedgenes <- function(geneinfodir,noputative=F,onlymetabolic=F){
  suppressMessages(library(dplyr))
  geneinfo <- read.delim(geneinfodir,header = T)
  df       <- geneinfo
  if (noputative)     df <- filter(df,GENENAME!="")
  if (onlymetabolic)  df <- filter(df,!is.na(METABOLIC))
  tokeep  <- as.character(df$GENEID)
  tokeep
}