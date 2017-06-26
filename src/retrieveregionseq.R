retrieveregionseq <- function(genome="sacCer3",
                              region=c("promoter","coding","chromosome"),
                              selgenes=NULL, #List of gene to filter from genome
                              filename=NULL){
  
  #Check if results exist
  if (!is.null(filename)){ 
    if (file.exists(filename)){
      message("regfile.rda found in: ", filename,". Skipping...")
      load(filename)
      return(res)
    }
  }
  
  #Params
  region   <- match.arg(region)
  
  #Load correct genome
  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(BSgenome))
  
  loadTxdb <- function(genome){
    if (genome == "sacCer3"){
      suppressMessages(library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene))
      txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
    } else {
      txdb <- makeTxDbFromUCSC(genome=genome, tablename="knownGene")  
    }
    txdb
  }
  txdb <- loadTxdb(genome)
  
  loadBS <- function(genome){
    if (genome == "sacCer3"){
      suppressMessages(library(BSgenome.Scerevisiae.UCSC.sacCer3))
      bsgenome <- BSgenome.Scerevisiae.UCSC.sacCer3
    } else {
      bsgenome <- getBSgenome(genome, masked=FALSE)
    }
    bsgenome
  }
  bsgenome <- loadBS(genome)
  
  #Retrieve SGD unique genes
  sgd.genome    <- as.data.frame(genes(txdb))
  sgd.genes     <- sgd.genome[,6]  #List of SGD genes - hard coded'
  if (!is.null(selgenes)){
    selgenes <- unlist(strsplit(selgenes,","))
    if (any(selgenes %in% sgd.genes)){
      sgd.genes   <- intersect(selgenes,sgd.genes)
    } else {
      stop("No genes in genome found within selected list of genes. Must be, e.g., '",sgd.genes[1],",'",sgd.genes[2],"'")
    }
  }
  sgd.genes.u   <- sgd.genes[!duplicated(sgd.genes)]
  
  #Define region (unless "chromosome")
  if (region!="chromosome"){
    if (region=="promoter"){
      upstream   <- 500
      downstream <- rep(0,length(sgd.genes))
        names(downstream) <- sgd.genes
    } else if (region == "coding"){
      upstream   <- 0
      downstream <- sgd.genome[,3]-sgd.genome[,2]
      names(downstream) <- sgd.genes
    }
    message("Extracting sequences from: ",region," region - Upstream: ",upstream, " Downstream (median): ",median(downstream))
  
    #Fetch target region sequence
    fetchregionseq <- function(gene,upstream,downstream,txdb,bsgenome){
      chromosomal.loc <- transcriptsBy(txdb, by="gene") [gene]
      regionseq     <- NA
      tryCatch({
        regionseq <- getPromoterSeq(chromosomal.loc, bsgenome, 
                                      upstream=upstream, 
                                      downstream=downstream[gene])
        },
        error=function(e){}
        )
      regionseq.ch  <- as.character(regionseq[[1]])[[1]]
      regionseq.ch
    }

  #Apply to each gene or chromosome
    message("Fetching region seq for ",length(sgd.genes.u)," genes...")
    regions.all <- sapply(sgd.genes.u, fetchregionseq,upstream,downstream,txdb,bsgenome)
  } else {
    chrs <- seqnames(bsgenome)
    message("Fetching region seq for ",length(chrs)," chromosomes...")
    regions.all <- sapply(chrs,function(x){as.character(bsgenome[[x]])})
  }  
  
  #Filter out duplicates or null seqs
  filter    <- intersect(which(!duplicated(regions.all)),which(!is.na(regions.all)))
  message(length(filter)," region seq are non-null or unique")
  
  #Return
  regions   <- regions.all[filter]
  if (region!="chromosome"){#Use chr names if region="chromosome"
    genes <- sgd.genes[filter]
  } else {
    genes <- names(regions)
  }
  res       <- list(regions=regions,genes=genes)
  if (!is.null(filename)) save(list = c("res"),file = filename)
  res
}