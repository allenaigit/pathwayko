.apply_parallel_pathway_analysis <- function(do.methods=NULL, DE=NULL, kpg=NULL, pwys=NULL, DataObject=NULL, targetPathways=NULL, named_FC=NULL, named_PVAL=NULL, spia.data.dir=NULL, path_to_xml=NULL){
  cat("\nINFO: pathway analysis...\n")
  joblist <- list()
  RONTOTOOLS_PE_RES <- NULL
  RONTOTOOLS_PDIS_RES <- NULL
  SPIA_RES <- NULL
  PADOG_RES <- NULL
  GSA_RES <- NULL
  SAFE_RES <- NULL
  GSEA_RES <- NULL

  if(do.methods$run_RONTOTOOLS_PE){
    RONTOTOOLS_PE_RES <- mcparallel(expr=.apply_RONTOTOOLS(fc=DE, kpg=kpg, named_FC=named_FC, named_PVAL=named_PVAL, targetPathways=targetPathways, verbose=FALSE, flag="PE"),
          name="RONTOTOOLS_PE_RES", 
          mc.set.seed = TRUE, 
          silent = TRUE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

  if(do.methods$run_RONTOTOOLS_PDIS){
    RONTOTOOLS_PDIS_RES <- mcparallel(expr=.apply_RONTOTOOLS(fc=DE, kpg=kpg, named_FC=named_FC, named_PVAL=named_PVAL, targetPathways=targetPathways, verbose=FALSE, flag="PDIS"),
          name="RONTOTOOLS_PDIS_RES", 
          mc.set.seed = TRUE, 
          silent = TRUE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

  if(do.methods$run_SPIA){
    data.dir <- spia.data.dir
    SPIA_RES <- mcparallel(expr=.apply_SPIA(de=DE, named_FC=named_FC, data.dir=data.dir, targetPathways=targetPathways, verbose=FALSE),
          name="SPIA_RES", 
          mc.set.seed = TRUE, 
          silent = TRUE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

  if(do.methods$run_PADOG){
    PADOG_RES <- mcparallel(expr=.apply_PADOG(DataObject = DataObject, targetPathways = targetPathways, pwys = pwys, verbose=FALSE),
          name="PADOG_RES", 
          mc.set.seed = TRUE, 
          silent = TRUE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

  if(do.methods$run_GSA){
    GSA_RES <- mcparallel(expr=.apply_GSA(DataObject = DataObject, targetPathways = targetPathways, pwys = pwys),
          name="GSA_RES", 
          mc.set.seed = TRUE, 
          silent = TRUE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

  if(do.methods$run_SAFE){
    SAFE_RES <- mcparallel(expr=.apply_SAFE(DataObject = DataObject, targetPathways = targetPathways, pwys = pwys),
          name="SAFE_RES", 
          mc.set.seed = TRUE, 
          silent = TRUE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

 
  if(do.methods$run_GSEA){
    GSEA_RES <- mcparallel(expr=.apply_GSEA(DataObject = DataObject, targetPathways = targetPathways, pwys = pwys),
          name="GSEA_RES", 
          mc.set.seed = TRUE, 
          silent = FALSE, 
          mc.affinity = NULL, 
          mc.interactive = FALSE, 
          detached = FALSE)
  }

  # Collect all jobs
  joblist <- list(
                SPIA_RES,
                RONTOTOOLS_PE_RES,
                RONTOTOOLS_PDIS_RES,
                PADOG_RES,                  
                GSA_RES,        
                GSEA_RES,
                SAFE_RES
              )
  collected <- mccollect(jobs=joblist[!(sapply(joblist,is.null))], wait = TRUE, timeout = 0, intermediate = FALSE)
  if(!length(collected)>0){
    cat("WARNING: no pathway analysis performed.\n")
    return(NULL)
  }
  #save(collected,file="collected_debug.RData")
  return(collected)
}



.apply_SPIA <- function(de, named_FC, data.dir="./", targetPathways, verbose){
  all <- names(named_FC)
  RES <- SPIA::spia(de = de, all = all, organism = "mmu", data.dir = data.dir, 
                    pathids = NULL, nB = 2000, plots = FALSE, verbose = verbose, 
                    beta = NULL, combine = "fisher")
  RES$Response <- rep(FALSE, length(RES$ID))
  for(i in 1:length(targetPathways)){
    # target pathway pathway present
    if(any(grepl(targetPathways[i],RES$ID))){
      RES[grep(targetPathways[i],RES$ID),"Response"] <- TRUE
    }
  }
  Rank = 1:length(RES$ID)
  RES <- cbind(RES, Predictor = RES$pG, Rank)
  return(RES)
}



.apply_RONTOTOOLS <- function(fc, kpg, named_FC, named_PVAL, targetPathways, verbose, flag="PE"){
  pv <- named_PVAL[names(fc)]
  # prepend "mmu:" to match examples in ROntoTools' manual
  ref <- unlist(lapply(names(named_FC),function(x){paste("mmu:",x,sep="")}))
  names(fc) <- sapply(names(fc),function(x){paste("mmu:",x,sep="")})
  names(pv) <- sapply(names(pv),function(x){paste("mmu:",x,sep="")})

  RES <- NULL
  # pe sub method
  if(toupper(flag)=="PE"){
    # set weight if running PE
    kpg <- ROntoTools::setNodeWeights(kpg, weights = ROntoTools::alphaMLG(pv), defaultWeight = 1)
    RES <- ROntoTools::pe(x = fc, graphs = kpg, ref = ref, nboot = 2000, verbose = verbose)

  }
  # pDis sub method
  else if(toupper(flag)=="PDIS"){
    RES <- ROntoTools::pDis(x = fc, graphs = kpg, ref = ref, nboot = 2000, verbose = verbose)
  }else{
    stop("ERROR: unknown flag in .apply_RONTOTOOLS.")
  }

  RES <- Summary(RES)

  RES$ID <- gsub("^path:mmu","",rownames(RES))
  RES$Response <- rep(FALSE, length(rownames(RES)))

  for(i in 1:length(targetPathways)){
    # target pathway pathway present
    if(any(grepl(targetPathways[i],RES$ID))){
      RES[grep(targetPathways[i],RES$ID),"Response"] <- TRUE
    }
  }
  Predictor <- RES$pComb
  Predictor[is.na(Predictor)] <- 1.0
  Rank = 1:length(RES$ID)
  RES <- cbind(RES, Predictor = Predictor, Rank)
  return(RES)
}



.apply_SAFE <- function(DE=NULL, DataObject=NULL, pwys=NULL, targetPathways=NULL){

  gslist <- lapply(pwys, function(X)
    {
      .getGenesByPwy <- function(Y)
      {
        ts <- vapply(nodes(Y), getType, character(1))
        genes <- unique(unlist(lapply(nodes(Y)[ts == "gene"], getName)))
        return(genes)
      }
      genes <- .getGenesByPwy(X)
      genes <- sub("^[a-z]{3}:", "", genes)
      genes <- sort(genes)
      return(genes)
    })
  names(gslist) <- lapply(names(gslist), function(X, pwys){
      return(pwys[[X]]@pathwayInfo@number)
    }, pwys=pwys)

  C.mat <- getCmatrix(keyword.list=gslist, as.matrix=TRUE)
  X.mat <- as.data.frame(exprs(DataObject$eset))
  geneInBoth <- intersect(rownames(X.mat),rownames(C.mat))
  # for X.mat and C.mat need to have same dimension
  C.mat <- C.mat[geneInBoth,]
  X.mat <- X.mat[geneInBoth,]
  y.vec <- ifelse(DataObject$eset$description=="case",1,0)
  
  #RES <- safe::safe(X.mat=X.mat, y.vec=y.vec, C.mat=C.mat, Pi.mat=2000, global="Fisher", args.global=list(one.sided=FALSE, genelist.length=length(DE)))
  RES <- safe(X.mat=X.mat, y.vec=y.vec, C.mat=C.mat, Pi.mat=2000)
  RES <- safe.toptable(RES,number=length(gslist),description=FALSE)
  RES <- RES[order(RES$P.value),]
  RES$ID <- RES$GenesetID

  RES$Response <- rep(FALSE, length(RES[,"ID"]))
  for(i in 1:length(targetPathways)){
    # target pathway pathway present
    if(any(grepl(targetPathways[i],RES[,"ID"]))){
      RES[grep(targetPathways[i],RES$ID),"Response"] <- TRUE
    }
  }
  Predictor <- as.numeric(RES$P.value)
  Predictor[is.na(Predictor)] <- 1.0
  Rank = 1:length(RES$ID)
  RES <- cbind(RES, Predictor = Predictor, Rank)

  return(RES)
}



.apply_GSEA <- function(DE=NULL, DataObject=NULL, pwys=NULL, targetPathways=NULL){
  ranks <- DataObject$exprTable[,c("ID","t")]
  ranks <- setNames(ranks$t, ranks$ID)


  gslist <- lapply(pwys, function(X)
    {
      .getGenesByPwy <- function(Y)
      {
        ts <- vapply(nodes(Y), getType, character(1))
        genes <- unique(unlist(lapply(nodes(Y)[ts == "gene"], getName)))
        return(genes)
      }
      genes <- .getGenesByPwy(X)
      genes <- sub("^[a-z]{3}:", "", genes)
      genes <- sort(genes)
      return(genes)
    })
  names(gslist) <- lapply(names(gslist), function(X, pwys){
      return(pwys[[X]]@pathwayInfo@number)
    }, pwys=pwys)

  
  RES <- fgsea(pathways=gslist, stats=ranks, eps=0.0)
  RES <- RES[order(RES$pval),]
  RES$leadingEdge <- NULL
  RES$ID <- RES$pathway
  RES <- as.data.frame(RES)
  RES$Response <- rep(FALSE, length(RES[,"ID"]))
  for(i in 1:length(targetPathways)){
    # target pathway pathway present
    if(any(grepl(targetPathways[i],RES[,"ID"]))){
      RES[grep(targetPathways[i],RES$ID),"Response"] <- TRUE
    }
  }
  Predictor <- as.numeric(RES$pval)
  Predictor[is.na(Predictor)] <- 1.0
  Rank = 1:length(RES$ID)
  RES <- cbind(RES, Predictor = Predictor, Rank)

  return(RES)
}



# Certain set based analysis doesn't require nor consider user supplied DE genes
# Hence their results does not depend on user supplied DE genes, and such analysis
# will be run without user supplied DE genes
.apply_PADOG<- function(DataObject = NULL, targetPathways = NULL, pwys = NULL, verbose=TRUE){ 
  # Modified from PADOG manual
  esetm <- exprs(DataObject$eset)
  # Convert expression matrix rownames back to probe ID for padog
  suppressMessages(require(DataObject$db,character.only=TRUE))
  
  # convert rownames of expr matrix from EntrezID back to probe IDs          
  rownames(esetm) <- suppressMessages(select(get(DataObject$db),rownames(esetm),columns="PROBEID",keytype="ENTREZID")[,2])

  group <- unlist(lapply(colnames(exprs(DataObject$eset)),function(X){
                  if(X%in%DataObject$control){
                    return("c")
                  }else if(X%in%DataObject$case){
                    return("d")
                  }else{
                    return("NA")
                  }
                }))

  paired <- FALSE
  block <- NULL
  targetgs <- targetPathways
  annotation <- DataObject$db
  gslist <- lapply(pwys, function(X)
    {
      .getGenesByPwy <- function(Y)
      {
        ts <- vapply(nodes(Y), getType, character(1))
        genes <- unique(unlist(lapply(nodes(Y)[ts == "gene"], getName)))
        return(genes)
      }
      genes <- .getGenesByPwy(X)
      genes <- sub("^[a-z]{3}:", "", genes)
      genes <- sort(genes)
      return(genes)
    })

  gs.names <- unlist(lapply(names(gslist), function(X, pwys)
    {
      return(pwys[[X]]@pathwayInfo@title)
    },pwys=pwys))
  names(gslist) <- lapply(names(gslist), function(X, pwys){
      return(pwys[[X]]@pathwayInfo@number)
    }, pwys=pwys)
  names(gs.names) <- names(gslist)

  RES <- padog(esetm = esetm,
                      group = group,
                      paired = paired,
                      block = block,
                      targetgs = targetgs,
                      annotation = annotation,
                      gslist = gslist,
                      gs.names = gs.names,
                      organism = "mmu",
                      verbose = verbose,
                      Nmin = 3,
                      # NI acts like an upper bound rather than a strict requirement. 
                      # e.g. PADOG may choose to run less iterations disregarding NI parameter
                      NI = 2000, 
                      plots = FALSE,
                      dseed = 1,
                      parallel = FALSE
                      )
  RES <- RES[order(RES$Ppadog),]
  RES$Response <- rep(FALSE, length(RES$ID))
  for(i in 1:length(targetPathways)){
    # target pathway pathway present
    if(any(grepl(targetPathways[i],RES$ID))){
      RES[grep(targetPathways[i],RES$ID),"Response"] <- TRUE
    }
  }
  Rank = 1:length(RES$ID)
  RES <- cbind(RES, Predictor = RES$Ppadog, Rank)

  return(RES)
}


.apply_GSA <- function(DataObject=NULL, pwys=NULL, targetPathways=NULL){
  x <- exprs(DataObject$eset)
  genenames <- rownames(x)
  y <- ifelse(DataObject$eset$description=="control",1,2)

  gslist <- lapply(pwys, function(X)
    {
      .getGenesByPwy <- function(Y)
      {
        ts <- vapply(nodes(Y), getType, character(1))
        genes <- unique(unlist(lapply(nodes(Y)[ts == "gene"], getName)))
        return(genes)
      }
      genes <- .getGenesByPwy(X)
      genes <- sub("^[a-z]{3}:", "", genes)
      genes <- sort(genes)
      return(genes)
    })
  names(gslist) <- lapply(names(gslist), function(X, pwys){
      return(pwys[[X]]@pathwayInfo@number)
    }, pwys=pwys)

  resp.type <- "Two class unpaired"

  RES <- GSA(x=x, y=y, nperms=2000, genesets=gslist, resp.type=resp.type, genenames=genenames)

  ps <- cbind(RES$pvalues.lo, RES$pvalues.hi)
  ps <- 2 * apply(ps, 1, min)
  scores <- RES$GSA.scores
  RES <- cbind(names(gslist), scores, ps)
  colnames(RES) <- c("ID","SCORE","PVAL")
  rownames(RES) <- names(gslist)
  RES <- as.data.frame(RES)
  RES <- RES[order(RES$PVAL),]



  RES$Response <- rep(FALSE, length(RES[,"ID"]))
  for(i in 1:length(targetPathways)){
    # target pathway pathway present
    if(any(grepl(targetPathways[i],RES[,"ID"]))){
      RES[grep(targetPathways[i],RES$ID),"Response"] <- TRUE
    }
  }
  Predictor <- as.numeric(RES$PVAL)
  Predictor[is.na(Predictor)] <- 1.0
  Rank = 1:length(RES$ID)
  RES <- cbind(RES, Predictor = Predictor, Rank)

  return(RES)
}

