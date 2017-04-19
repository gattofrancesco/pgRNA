aligngRNAtoregions <- function(gRNAs.ext=NULL,regions=NULL,gRNAfile=NULL,regfile=NULL,
                                 nmismatches,nchunks=1,chunk=1,filename=NULL){
  suppressMessages(library(Biostrings))
  
  #Check gRNAs.ext and regions
  if (is.null(gRNAs.ext)&is.null(regions)){
    if (is.null(gRNAfile)|is.null(regfile)){
      stop("Either gRNAfile and regfile or gRNAs.ext and regions must be supplied")}
    else {
      load(gRNAfile)
      load(regfile)
      regions <- res[["regions"]]
      genes     <- res[["genes"]]    
    }
  }
  
  #Given n mismatches, find gRNA alignments to target regions
  message("Matching all gRNAs to target regions for chunk #",chunk," of ",nchunks,
          " with max ",nmismatches," mismatch(es)")

  #Divide gRNAs into n chunks and align only # chunk
  chunksize   <- ceiling(length(gRNAs.ext)/nchunks)
  chunkinds   <- split(1:length(gRNAs.ext), 
                       ceiling(seq_along(1:length(gRNAs.ext))/chunksize))
  chunkinds_inscript <- chunkinds[[chunk]]
  
  #Start alignment
  gRNAs.set  <- DNAStringSet(gRNAs.ext[chunkinds_inscript])
  gRNAs.comp <- reverseComplement(gRNAs.set)
  if (nmismatches==0){
    #Use exact matching with dictionary
    message("Running exact matching on cis strand...")
    gRNAs.dict      <- PDict(gRNAs.set)
    nhitsxregion.cis   <- sapply(regions,function(x){
      countPDict(gRNAs.dict,DNAString(x),
                 algorithm="auto",max.mismatch=0,
                 with.indels=FALSE, fixed=TRUE)
      }
      )
    message("Running exact matching on trans strand...")
    gRNAs.comp.dict <- PDict(gRNAs.comp)
    nhitsxregion.trans <- sapply(regions,function(x){
      countPDict(gRNAs.comp.dict,DNAString(x),
                 algorithm="auto",max.mismatch=0,
                 with.indels=FALSE, fixed=TRUE)
    }
    )
    nhitsxregion <- nhitsxregion.cis + nhitsxregion.trans
  } else if (nmismatches>0){
    #Use unprocessed dictionary
    message("Running inexact matching on cis strand...")
    nhitsxregion.cis   <- sapply(regions,function(x){
      countPDict(gRNAs.set,DNAString(x),
                 algorithm="auto",max.mismatch=nmismatches,
                 with.indels=FALSE, fixed=TRUE)
    }
    )
    message("Running inexact matching on trans strand...")
    nhitsxregion.trans <- sapply(regions,function(x){
      countPDict(gRNAs.comp,DNAString(x),
                 algorithm="auto",max.mismatch=nmismatches,
                 with.indels=FALSE, fixed=TRUE)
    }
    )
    nhitsxregion <- nhitsxregion.cis + nhitsxregion.trans
  }

  #Postprocess: transform into gRNAs x region
  hitsindexes <- t(nhitsxregion)
  colnames(hitsindexes) <- as.character(gRNAs.set)
  
  #Return
  if (!is.null(filename)){
    save(list=c("gRNAs.ext","hitsindexes","regions"), 
         file = paste0(filename,chunk,"_nmismatch_",nmismatches,".rda"))
  }
  hitsindexes
}