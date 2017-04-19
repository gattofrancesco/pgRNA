# Define args
resdir      <- NULL
args=(commandArgs(TRUE)) #
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#Issue command
source("lumphitsindexes.R")
if (!grepl("hitlist",dir(resdir))){
  message("No hitlist found in directory: ", resdir,". Skipping...")
} else {
  lumphitsindexes(resdir)
}