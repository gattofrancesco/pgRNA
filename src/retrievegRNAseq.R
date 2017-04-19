retrievegRNAseq <- function(regions=NULL,genes=NULL,regfile=NULL,
                            gRNA.size=20,PAM.size=3,minGC=0.1,filename=NULL){
  #Check if results exist
  if (!is.null(filename)){ 
    if (file.exists(filename)){
       message("gRNAfile.rda found in: ", filename,". Skipping...")
       load(filename)
       return(gRNAs.ext)
    } 
  }
  
  suppressMessages(library(Biostrings))
  suppressMessages(library(seqinr))
  
  #Check regions
  if (is.null(regions)&is.null(genes)){
    if (is.null(regfile)){
      stop("Either regfile or genes and regions must be supplied")}
    else {
      load(regfile)
      regions <- res[["regions"]]
      genes   <- res[["genes"]]    
      }
  }
  
  ## Find gRNA
  gRNA.pattern <- ""
  starts.gRNA  <- ""
  ends.gRNA    <- ""
  gRNA.tab     <- ""
  gRNAs.ext    <- c()
  extSeq.4     <- ""
  seqleft      <- length(regions)
  
  #No idea what the inner loop does - rename variable and comment
  for (region in regions){
    PAM = "GG"
    seq.len  <- nchar(region)
    #Search position of PAM
    pos.PAMs <- unlist(gregexpr(PAM, region, perl = TRUE,
                                ignore.case = TRUE, fixed = FALSE)) -1
    #Remove less than 20 nucleotides PAM pos.
    pos.PAMs <- pos.PAMs[pos.PAMs != -1 & pos.PAMs > gRNA.size]  
    #Start and end of gRNAs
    starts.gRNA <- pos.PAMs - gRNA.size 
    ends.gRNA   <- starts.gRNA + gRNA.size - 1
    seq       <- as.character(Views(region,
                              start = starts.gRNA,
                              end = ends.gRNA))
    extSeq    <- as.character(Views(region,
                                 start = starts.gRNA ,
                                 end = ends.gRNA + PAM.size))
    gRNAs.ext <- c(gRNAs.ext,extSeq)
    PAM3="CC"
    pos.PAMs.3 <- unlist(gregexpr(PAM3, region, perl = TRUE, 
                                  ignore.case = TRUE, fixed = FALSE)) +1
    pos.PAMs.3 <- pos.PAMs.3[pos.PAMs.3 != -1 & pos.PAMs.3 > gRNA.size & pos.PAMs.3 < seq.len - 20] #Hard-coded?
    starts.gRNA.3 <- pos.PAMs.3 + 2
    ends.gRNA.3 <- (starts.gRNA.3 + gRNA.size - 1)
    seq.3   <- as.character(Views(region, 
                                start = starts.gRNA.3,
                                end = ends.gRNA.3))
    extSeq.3 <- as.character(Views(region,
                                   start = starts.gRNA.3 - PAM.size ,
                                   end = ends.gRNA.3 ))
    if (length(extSeq.3 > 0)){
      
      for (k in 1:length(extSeq.3)){ 
        mattt = ""
        tmpp =""
        tmpp <- DNAString(extSeq.3[k])
        tmpp <- reverseComplement(tmpp)
        tmmp <- as.character(tmpp)
        mattt[1] <- tmmp
        extSeq.4 <- c(extSeq.4,mattt[1])
        
      }
      extSeq.4 <- extSeq.4[-1]
      gRNAs.ext <- c(gRNAs.ext,extSeq.4)
    }
    extSeq.4=""
    seqleft <- seqleft - 1
    message(seqleft," sequences remaining...")
  }
  
  ## Filter bad quality gRNA 
  GCcont <- sapply(gRNAs.ext,function(x){
    myseq <- s2c(x)
    gcseq <- GC(myseq)
    }
    )
  message(sum(GCcont<=minGC)," gRNAs have GC content <= ",minGC,". Discarded.")
  gRNAs.ext <- gRNAs.ext[GCcont>minGC]
  
  #Return
  if (!is.null(filename)) save(list = c("gRNAs.ext"),file = filename)
  gRNAs.ext  
}