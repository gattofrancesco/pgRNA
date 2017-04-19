analyzealignments <- function(hitsindexes,region,nmismatches,regions=NULL,genes=NULL,regfile=NULL,
                              toprank=100,maxgRNA=2*10^5,selgenes=NULL,filename=NULL){
  suppressMessages(library(ggplot2))
  suppressMessages(library(reshape2))
  suppressMessages(library(seqinr))
  
  #Check target sequences
  if (is.null(regions)&is.null(genes)){
    if (is.null(regfile)){
      stop("Either regfile or genes and regions must be supplied")}
    else {
      load(regfile)
      regions <- res[["regions"]]
      genes   <- res[["genes"]]    
    }
  }
  
  #Preprocess results
  hits         <- hitsindexes
  if (!is.null(selgenes)){
    ind2keep <- match(selgenes,genes)[!is.na(match(selgenes,genes))]
    message(length(ind2keep)," genes in target sequences were retained")
    hits  <- hits[ind2keep,]
    genes <- genes[ind2keep]
    regions <- regions[genes]
  }
  nhitspergRNA <- colSums(hits!=0)
  
  #Plot 1: n hits per gRNA
  m <- melt(nhitspergRNA)
  p <- ggplot(m,aes(x=value)) + geom_histogram(binwidth = 1) + theme_bw() + xlab(paste("# hits per gRNA")) + ylab("Counts") +
    theme(axis.text = element_text(size=14),axis.title = element_text(size=14)) + #ylim(0,10)
    scale_y_continuous(trans="log10")
  print(p)
  filename <- paste0("plots/nhitsxgRNA_",region,"_nmismatches",nmismatches,ifelse(!is.null(selgenes),"_annotatedonly.pdf",".pdf"))
  ggsave(filename = filename,p)
  
  #Plot 2: density hits per gRNA
  p <- ggplot(m,aes(x=value)) + geom_density() + theme_bw() + xlab(paste("N hits per gRNA")) + ylab("Frequency")
  print(p)
  filename <- paste0("plots/densityhitsxgRNA_",region,"_nmismatches",nmismatches,ifelse(!is.null(selgenes),"_annotatedonly.pdf",".pdf"))
  ggsave(filename = filename,p)

  #Plot 3: top gRNAs according to n of hits
  nhitspergRNA.sort <- sort(nhitspergRNA,decreasing = T)
  nhitspergRNA.sort <- nhitspergRNA.sort[1:maxgRNA] #Discard too promiscous
  nhitspergRNA.sort <- nhitspergRNA.sort[!duplicated(names(nhitspergRNA.sort))][1:toprank]
  m2 <- data.frame(gRNA=factor(names(nhitspergRNA.sort),levels=names(nhitspergRNA.sort)),
                   melt(nhitspergRNA.sort))
  p <- ggplot(m2,aes(x=gRNA,y=value)) + geom_bar(stat = "identity",width = 0.5) + 
    theme_bw() + xlab(paste("Top",toprank,"multitarget gRNAs")) + ylab(paste("N",region,"sequence hits")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text = element_text(size=14),axis.title = element_text(size=14))
  print(p)
  filename <- paste0("plots/top",toprank,"multitargetgRNAs_",region,"_nmismatches",nmismatches,ifelse(!is.null(selgenes),"_annotatedonly.pdf",".pdf"))
  ggsave(filename = filename,p)
  
  #Tab 1-2: Realign top gRNAs tp target sequences
  TOPgRNA <- names(which(nhitspergRNA.sort > 1))
  if (length(TOPgRNA > 1)){
    print("Matching top gRNAs to target sequences...")
    
    #Help function to align each gRNA and score alignment
    suppressMessages(suppressMessages(library(Biostrings)))
    scoreAlignment <- function(gRNA,regions){
      tmp <- sapply(regions,function(x){
        matchPattern(gRNA, x, algorithm="auto",max.mismatch=nmismatches, with.indels=FALSE, fixed=TRUE)
      })
    }
    
    best.tab = matrix(nrow=length(TOPgRNA),ncol=m2[1,2])
    rownames(best.tab) <- TOPgRNA
    
    for (gRNA in TOPgRNA){
      print(paste("Doing...",gRNA))
      scores         <- scoreAlignment(gRNA,regions)
      nhitsxregion <- sapply(scores,function(x) length(start(x)))
      nhitsxregion <- as.matrix(nhitsxregion)
      posiz <- which(nhitsxregion != 0)
      posiz <- as.matrix(posiz)
      best.tab[gRNA,1:length(posiz)] <- genes[posiz]
    }
  }
  
  if (is.null(filename)){
    filename   <- paste0("res/top",toprank,"gRNA_hitlistin",region,"_nmismatches",nmismatches,
                       ifelse(!is.null(selgenes),"_annotatedonly.csv",".csv"))
  }
  write.csv(best.tab, filename)
  }