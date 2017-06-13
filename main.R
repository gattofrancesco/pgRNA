#!/usr/bin/env Rscript
rm(list=ls())

#Install Rmarkdown
install.packages("rmarkdown",dependencies = T)
library(rmarkdown)

# Define default args
nmismatches = 0
chunk       = 1
nchunks     = 1
region      = "coding"
toprank     = 10
maxgRNA     = 200000
interimdir  = "interim/" 
minGC       = 0.1 
selgenes    = 'YOR317W,YER015W,YIL009W,YMR246W,YMR008C,YMR006C,YOL011W'
genome      = "sacCer3"

# Read args from command line
args=(commandArgs(TRUE)) #
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#Launch RMarkdown script
library(rmarkdown)
library(knitr)
render(input="main.Rmd",
       params = 
       list(nmismatches=nmismatches,
            chunk=chunk,
            nchunks=nchunks,
            region=region,
            toprank=toprank,
            maxgRNA=maxgRNA,
            interimdir=interimdir, 
            minGC=minGC, 
            selgenes=selgenes,
            genome=genome),
            output_format = 'html_document')