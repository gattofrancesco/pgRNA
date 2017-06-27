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
  r   <- DNAStringSet(regions)
  
  #Forward strand
  PAM <- DNAString("GG") #Protospacer Adjacent Motif 3 nucleotides
  
  #Search position of PAM in each region
  pos.PAMs <- vmatchPattern(PAM,r,
                algorithm="auto",max.mismatch=0,
                with.indels=FALSE, fixed=TRUE)
  
  #Filter out PAM sequence if located at <20bp from the first nucleotide
  pos.PAMs.f <- lapply(pos.PAMs,function(x){y <- x[start(x) != -1 & start(x) > gRNA.size]})
  
  #Define gRNAs region
  starts.gRNA <- lapply(pos.PAMs.f,function(x){start(x) - gRNA.size - 1})    #position of all the starting nucleotide of each gRNA
  ends.gRNA   <- lapply(pos.PAMs.f,function(x){start(x) - 2}) #position of the end of each gRNAs
  gRNAs.ext.f   <- c()
  for (i in 1:length(regions)){
    allgRNAs <- Views(regions[[i]], ##retrieve the sequence between each start and end which corresponds to the 20 gRNA sequences 
                      start = starts.gRNA[[names(regions)[i]]] ,
                      end   = ends.gRNA[[names(regions)[i]]] + PAM.size)
    region.length <- length(subject(allgRNAs))
    allgRNAs.f    <- allgRNAs[!(start(allgRNAs)==0|end(allgRNAs)>region.length)] #gRNAs is fully contained in region
    gRNAs.ext.f   <- c(gRNAs.ext.f,as.character(allgRNAs.f))
  }
  
  #Reverse strand
  PAM.c <- DNAString("CC") #Protospacer Adjacent Motif 3 nucleotides
  
  #Search position of PAM in each complement region
  pos.PAMs.c <- vmatchPattern(PAM.c,r,
                            algorithm="auto",max.mismatch=0,
                            with.indels=FALSE, fixed=TRUE)
  
  #Define gRNAs region
  starts.gRNA <- lapply(pos.PAMs.c,function(x){end(x) + 2})    #position of all the starting nucleotide of each gRNA
  ends.gRNA   <- lapply(pos.PAMs.c,function(x){end(x) + gRNA.size + 1}) #position of the end of each gRNAs
  gRNAs.ext.c   <- c()
  for (i in 1:length(regions)){
    allgRNAs <- Views(regions[[i]], ##retrieve the sequence between each start and end which corresponds to the 20 gRNA sequences 
                      start = starts.gRNA[[names(regions)[i]]] - PAM.size,
                      end   = ends.gRNA[[names(regions)[i]]])
    region.length <- length(subject(allgRNAs))
    allgRNAs.f    <- allgRNAs[!(start(allgRNAs)==0|end(allgRNAs)>region.length)] #gRNAs is fully contained in region
    gRNAs.ext.c   <- c(gRNAs.ext.c,as.character(allgRNAs.f))
  }
  gRNAs.ext.c <- as.character(reverseComplement(DNAStringSet(gRNAs.ext.c)))
  gRNAs.ext   <- c(gRNAs.ext.f,gRNAs.ext.c)
  
  ## Filter bad quality gRNA 
  GCcont <- sapply(gRNAs.ext,function(x){
    myseq <- s2c(x)
    gcseq <- GC(myseq)
    }
    )
  message(sum(GCcont<=minGC)," gRNAs have GC content <= ",minGC,". Discarded.")
  gRNAs.ext <- gRNAs.ext[GCcont>minGC]
  
  ## Filter A/T/C/G content 
  gRNAs.ext.alpha  <- alphabetFrequency(DNAStringSet(gRNAs.ext))
  gRNAs.ext.noATCG <- gRNAs.ext.alpha[,-which(colnames(gRNAs.ext.alpha)%in%c("A","C","T","G"))]
  ATCG.filter      <- apply(gRNAs.ext.noATCG==0,1,all)
  message(sum(!ATCG.filter)," gRNAs have at least one nucleotide different from ATCG. Discarded.")
  gRNAs.ext <- gRNAs.ext[ATCG.filter]
  
  #Return
  if (!is.null(filename)) save(list = c("gRNAs.ext"),file = filename)
  gRNAs.ext  
}