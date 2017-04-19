lumphitsindexes <- function(resdir,region,nmismatches){
  # resdir    <- "res/"
  indexes   <- grepl("hitlist",dir(resdir)) & 
               grepl(paste0("nmismatch_",nmismatches),dir(resdir)) &
               grepl(region,dir(resdir))
  files     <- dir(resdir)
  files.hit <- paste0(resdir,files[indexes])
  for (f in files.hit){
    message("Loading ",f)
    if (f == files.hit[1]){
      load(f)
      p        <- promoters
      hits_all <- hitsindexes
    } else {
      load(f)
      stopifnot(identical(p,promoters))
      hits_all <- cbind(hits_all,hitsindexes)
    }
  }
  hits_all
}