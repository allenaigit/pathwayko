generateROC <- function(collected=NULL,result_output=NULL,filename_prefix=NULL,targetPathways=NULL,roc.partial.se.lowerbound=NULL,roc.partial.sp.lowerbound=NULL){
  cat("\nINFO: plotting ROC...\n")
  ret <- NULL

  # Pull only Response and Predictor data for ROC
  jobs <- lapply(collected, function(X,t){
      ret <- list()
      diff <- length(t)-length(which(X[,"Response"]==TRUE))
      if(diff>0){
        ret$Response <- c(X[,"Response"],rep("TRUE",diff))
        ret$Predictor <- c(X[,"Predictor"],rep(1.0,diff))
      }else{
        ret$Response <- X[,"Response"]
        ret$Predictor <- X[,"Predictor"]
      }
      return(ret)
    },t=targetPathways)

  # Generate ROC objects and extract AUCs
  rocobjs <- lapply(jobs, function(X){
      roc(response=X$Response, predictor=X$Predictor, direction=">", percent=TRUE, quiet=TRUE, auc=TRUE, ci=TRUE, ci.of="auc")
    })
  aucs <- lapply(rocobjs, function(X){
      format(as.numeric(X$auc),digits=4,nsmall=2)
    })

  # Save AUCs
  auc.out <- data.frame(lapply(rocobjs, function(X){
              format(as.numeric(X$auc),digits=5,nsmall=3)
            }))
  colnames(auc.out) <- gsub("_RES$","_AUC",colnames(auc.out))
  rownames(auc.out) <- gsub("_$","",filename_prefix)
  write.csv(auc.out,file=paste0(result_output,filename_prefix,"AUC.csv"))
  ret$AUC <- auc.out

  # Save other stats
  # stats is based on the best threshold determined by pROC, which may choose multiple best points, hence multiple rows
  stats <- lapply(rocobjs,function(X){
      coords(
        X,
        "best",
        ret=c(  
            "threshold",
            "specificity",
            "sensitivity",
            "accuracy",
            "precision",
            "recall",
            "tpr",
            "fpr",
            "tnr",
            "fnr",
            "fdr",
            "npv",
            "ppv",
            "tn",
            "tp",
            "fn",
            "fp"
            ),
        transpose=FALSE)
    })
  stats <- rbind(unlist(stats))
  rownames(stats) <- gsub("_$","",filename_prefix)
  write.csv(stats,file=paste0(result_output,filename_prefix,"stats.csv"))
  ret$stats <- stats

  # stats1 forces only one best threshold by maximizing specificity
  stats1 <- lapply(rocobjs,function(X){
      tmp <- coords(
        X,
        "best",
        ret=c(  
            "threshold",
            "specificity",
            "sensitivity",
            "accuracy",
            "precision",
            "recall",
            "tpr",
            "fpr",
            "tnr",
            "fnr",
            "fdr",
            "npv",
            "ppv",
            "tn",
            "tp",
            "fn",
            "fp"
            ),
        transpose=FALSE)
      if(nrow(tmp)>1){
        tmp <- tmp[order(tmp$specificity,decreasing=TRUE),]
        tmp <- tmp[1,]
      }
      return(tmp)
    })
  stats1 <- rbind(unlist(stats1))
  rownames(stats1) <- gsub("_$","",filename_prefix)
  write.csv(stats1,file=paste0(result_output,filename_prefix,"stats1.csv"))
  ret$stats1 <- stats1

  stats2 <- lapply(rocobjs,function(X){
      coords(
        X,
        "best",
        ret=c(  
            "threshold",
            "specificity",
            "sensitivity",
            "accuracy",
            "precision",
            "recall",
            "tpr",
            "fpr",
            "tnr",
            "fnr",
            "fdr",
            "npv",
            "ppv",
            "tn",
            "tp",
            "fn",
            "fp"
            ),
        transpose=TRUE)
    })
  stats2 <- as.data.frame(stats2)
  colnames(stats2) <- gsub("_RES$","",colnames(stats2))
  write.csv(stats2,file=paste0(result_output,filename_prefix,"stats2.csv"))
  ret$stats2 <- stats2


  # Generate partial AUCs and save them
  pauc.sp.cor <- 
    lapply(jobs, function(X){
      roc(response=X$Response, 
        predictor=X$Predictor, 
        direction=">", 
        percent=TRUE, 
        quiet=TRUE, 
        auc=TRUE,
        partial.auc=c(100, roc.partial.sp.lowerbound), 
        partial.auc.correct=TRUE, 
        allow.invalid.partial.auc.correct=FALSE,
        partial.auc.focus="sp"
        )
    })
  pauc.sp.cor.out <- data.frame(lapply(pauc.sp.cor, function(X){
      format(as.numeric(X$auc),digits=5,nsmall=3)
    }))
  colnames(pauc.sp.cor.out) <- gsub("_RES$","_pAUC_SP",colnames(pauc.sp.cor.out))
  rownames(pauc.sp.cor.out) <- gsub("_$","",filename_prefix)
  write.csv(pauc.sp.cor.out,file=paste0(result_output,filename_prefix,"pAUC_SP_",roc.partial.sp.lowerbound,"_Corrected.csv"))
  ret$pauc_sp_cor <- pauc.sp.cor.out

  pauc.sp.ori <- 
  lapply(jobs, function(X){
    roc(response=X$Response, 
      predictor=X$Predictor, 
      direction=">", 
      percent=TRUE, 
      quiet=TRUE, 
      auc=TRUE,
      partial.auc=c(100, roc.partial.sp.lowerbound), 
      partial.auc.correct=FALSE, 
      allow.invalid.partial.auc.correct=FALSE,
      partial.auc.focus="sp"
      )
  })
  pauc.sp.ori.out <- data.frame(lapply(pauc.sp.ori, function(X){
      format(as.numeric(X$auc),digits=5,nsmall=3)
    }))
  colnames(pauc.sp.ori.out) <- gsub("_RES$","_pAUC_SP",colnames(pauc.sp.ori.out))
  rownames(pauc.sp.ori.out) <- gsub("_$","",filename_prefix)
  write.csv(pauc.sp.ori.out,file=paste0(result_output,filename_prefix,"pAUC_SP_",roc.partial.sp.lowerbound,"_Original.csv"))
  ret$pauc_sp_ori <- pauc.sp.ori.out


  pauc.se.cor <- 
  lapply(jobs, function(X){
    roc(response=X$Response, 
      predictor=X$Predictor, 
      direction=">", 
      percent=TRUE, 
      quiet=TRUE, 
      auc=TRUE,
      partial.auc=c(100, roc.partial.se.lowerbound), 
      partial.auc.correct=TRUE, 
      allow.invalid.partial.auc.correct=FALSE,
      partial.auc.focus="se"
      )
  })
  pauc.se.cor.out <- data.frame(lapply(pauc.se.cor, function(X){
      format(as.numeric(X$auc),digits=5,nsmall=3)
    }))
  colnames(pauc.se.cor.out) <- gsub("_RES$","_pAUC_SE",colnames(pauc.se.cor.out))
  rownames(pauc.se.cor.out) <- gsub("_$","",filename_prefix)
  write.csv(pauc.se.cor.out,file=paste0(result_output,filename_prefix,"pAUC_SE_",roc.partial.se.lowerbound,"_Corrected.csv"))
  ret$pauc_se_cor <- pauc.se.cor.out


  pauc.se.ori <- 
  lapply(jobs, function(X){
    roc(response=X$Response, 
      predictor=X$Predictor, 
      direction=">", 
      percent=TRUE, 
      quiet=TRUE, 
      auc=TRUE,
      partial.auc=c(100, roc.partial.se.lowerbound), 
      partial.auc.correct=FALSE, 
      allow.invalid.partial.auc.correct=FALSE,
      partial.auc.focus="se"
      )
  })
  pauc.se.ori.out <- data.frame(lapply(pauc.se.ori, function(X){
      format(as.numeric(X$auc),digits=5,nsmall=3)
    }))
  colnames(pauc.se.ori.out) <- gsub("_RES$","_pAUC_SE",colnames(pauc.se.ori.out))
  rownames(pauc.se.ori.out) <- gsub("_$","",filename_prefix)
  write.csv(pauc.se.ori.out,file=paste0(result_output,filename_prefix,"pAUC_SE_",roc.partial.se.lowerbound,"_Original.csv"))
  ret$pauc_se_ori <- pauc.se.ori.out



  # ROC Combined
  colors <- brewer.pal(max(length(rocobjs),3),"Set3")
  pdf(file=paste0(result_output,filename_prefix,"ROC_combined.pdf"),height=7,width=7)
  for(i in 1:length(rocobjs)){
    if(i==1){
      plot.roc(
            x=rocobjs[[i]],
            col=colors[i],
            quiet=TRUE, 
            main=paste0(filename_prefix,"ROC_combined")
          )
    }else{
      plot.roc(
            x=rocobjs[[i]],
            add=TRUE,
            col=colors[i],
            quiet=TRUE
          )
    }
  }
  legend("bottomright",
      legend=paste0(gsub("_RES$","",names(rocobjs))," (",aucs,")"),
      col=colors,
      lwd=3,
      cex=0.8)
  dev.off()

  # Individual ROCs with:
  # 1. AUC with 95% confidence interval and the best threshold with corresponding specificity and sensitivity
  # 2. Partial AUCs while locking sensitivity percentage in [100, roc.partial.se.lowerbound] and specificity percentage in [100, roc.partial.sp.lowerbound]
  # 3. Confidence intervals of both SE and SP with bar plots
  for(i in 1:length(rocobjs)){
    # 1
    pdf(file=paste0(result_output,filename_prefix,gsub("_RES$","_ROC",names(rocobjs)[[i]]),".pdf"),height=7,width=7)
    plot.roc(
          x=rocobjs[[i]],
          col="red",
          lwd=3,
          main=paste0(filename_prefix,gsub("_RES$","_ROC",names(rocobjs)[[i]])),
          print.thres="best", 
          print.thres.col="black", 
          print.thres.cex=2.0,
          print.thres.adj = c(0.5, 1.5),
          print.thres.pattern="Best p threshold:\n %.3f (%.1f%%, %.1f%%)", 
          print.thres.pattern.cex=1.0,
          reuse.auc=TRUE,
          print.auc=TRUE, 
          print.auc.col="black",
          print.auc.pattern="AUC: %.2f%% [%.1f%% - %.1f%%]",
          print.auc.cex=1.0,
          quiet=TRUE
        ) 
    legend("bottomright",
      legend=paste0(gsub("_RES$","",names(rocobjs[i]))),
      col="red",
      lwd=3,
      cex=1.0)        
    dev.off()

    # 2
    pdf(file=paste0(result_output,filename_prefix,gsub("_RES$","_pAUC",names(rocobjs)[[i]]),".pdf"),height=7,width=7)
    plot.roc(
          x=rocobjs[[i]]$response, 
          predictor=rocobjs[[i]]$predictor,
          direction=">",
          col="black",
          lwd=3,
          percent=TRUE,
          main=paste0(filename_prefix,gsub("_RES$","_pAUC",names(rocobjs)[[i]])),
          partial.auc=c(100, roc.partial.se.lowerbound), 
          partial.auc.correct=TRUE, 
          allow.invalid.partial.auc.correct=FALSE,
          partial.auc.focus="se",
          print.auc=TRUE, 
          print.auc.pattern=paste0("Corrected pAUC (100%%-",roc.partial.se.lowerbound,"%% SE):\n%.2f%%"), 
          print.auc.col="#008600",
          auc.polygon=TRUE, 
          auc.polygon.col="#008600",
          max.auc.polygon=TRUE, 
          max.auc.polygon.col="#00860022",
          quiet=TRUE
            ) 
    plot.roc(
          x=rocobjs[[i]]$response, 
          predictor=rocobjs[[i]]$predictor,
          direction=">",
          lwd=3,
          percent=TRUE,
          add=TRUE,
          type="n",
          partial.auc=c(100, roc.partial.sp.lowerbound), 
          partial.auc.correct=TRUE, 
          allow.invalid.partial.auc.correct=FALSE,
          partial.auc.focus="sp",
          print.auc=TRUE, 
          print.auc.y=40,
          print.auc.pattern=paste0("Corrected pAUC (100%%-",roc.partial.sp.lowerbound,"%% SP):\n%.2f%%"), 
          print.auc.col="#1c61b6",
          auc.polygon=TRUE, 
          auc.polygon.col="#1c61b6",
          max.auc.polygon=TRUE, 
          max.auc.polygon.col="#1c61b622",
          quiet=TRUE
            )
    dev.off()

    # 3
    pdf(file=paste0(result_output,filename_prefix,gsub("_RES$","_CI",names(rocobjs)[[i]]),".pdf"),height=7,width=7)
    rocobj <- plot.roc(
                x=rocobjs[[i]]$response,
                predictor=rocobjs[[i]]$predictor,
                direction=">",
                percent=TRUE,
                lwd=3,
                col="red",
                main=paste0(filename_prefix,gsub("_RES$","_ConfidenceInterval",names(rocobjs)[[i]])),
                ci=TRUE, 
                of="se", 
                specificities=seq(0, 100, 4), 
                ci.type="bars",
                quiet=TRUE,
                progress="none"
              )
    plot(ci.sp(rocobj, sensitivities=seq(0, 100, 4), progress="none"), type="bars")
    dev.off()
  }

  return(ret)
}



GenerateDEGenes<- function(method="HES1", KEGGgraph=NULL, DataObject=NULL, named_FC=NULL, named_PVAL=NULL, alpha=NULL, result_output=NULL, filename_prefix=NULL, mc.cores=1){
  DE <- NULL
  # HighEdgeS
  if(method%in%c("HES1","HES2","HES3")){
    if(is.null(KEGGgraph)||is.null(DataObject)||is.null(named_FC)||is.null(named_PVAL)){
      stop(paste0("ERROR: in GenerateDEGenes missing argument.\n"))
    }

    hes <- highedges(KEGGgraph, DataObject, named_FC, named_PVAL, mc.cores)
    mechanismGraph <- NULL

    if(method=="HES1"){
      DE <- hes$testData[1:sum(hes$histObject$counts[hes$changePoints[[1]]:length(hes$histObject$counts)])]
      DE <- unique(unlist(str_split(names(DE),pattern="[|]")))
      DE <- named_FC[DE]
      mechanismGraph <- generateMechanismGraph(hes=hes,CP_index=1)
    }
    if(method=="HES2"){
      DE <- hes$testData[1:sum(hes$histObject$counts[hes$changePoints[[2]]:length(hes$histObject$counts)])]
      DE <- unique(unlist(str_split(names(DE),pattern="[|]")))
      DE <- named_FC[DE]
      mechanismGraph <- generateMechanismGraph(hes=hes,CP_index=2)
    }
    if(method=="HES3"){
      DE <- hes$testData[1:sum(hes$histObject$counts[hes$changePoints[[3]]:length(hes$histObject$counts)])]
      DE <- unique(unlist(str_split(names(DE),pattern="[|]")))
      DE <- named_FC[DE]
      mechanismGraph <- generateMechanismGraph(hes=hes,CP_index=3)
    }

    pdf(file=paste(result_output,filename_prefix,"HighEdgeS_plot.pdf",sep=""))
    suppressWarnings(plot(mechanismGraph, main = paste("KO gene ", DataObject$KOgeneSymbol, " Mechanism Graph", sep = "")))
    dev.off()

    pdf(file=paste(result_output,filename_prefix,"Histogram_plot.pdf",sep=""))
    plot(hes$histObject, xlab = 'Edges Scores',main = paste("KO gene ", DataObject$KOgeneSymbol, " Edge Score Histogram", sep = ""))
    abline(v=hes$histObject$breaks[hes$changePoints[[1]]],col="blue")
    abline(v=hes$histObject$breaks[hes$changePoints[[2]]],col="green")
    abline(v=hes$histObject$breaks[hes$changePoints[[3]]],col="red")
    legend(x="center",legend=c("HES1","HES2","HES3"),col=c("blue","green","red"),lwd=1.2)
    dev.off()

  # Classical (percentage of genes in KEGG)
  }else if(method=="CLA"){
    if(is.null(DataObject)||is.null(named_FC)||is.null(alpha)){
      stop(paste0("ERROR: in GenerateDEGenes missing argument.\n"))
    }
    all <- DataObject$exprTable[DataObject$exprTable$ID%in%names(named_FC),]
    desiredNum <- ceiling(nrow(all)*alpha)

    DE1 <- all[all$adj.P.Val < 0.1, "ID"]
    if(length(DE1)>=desiredNum){
      DE <- DE1[1:desiredNum]
    }else{
      DE2 <- all[all$P.Value < 0.01, "ID"]
      if(length(DE2)>=desiredNum){
        DE <- DE2[1:desiredNum]
      }else{
        all <- all[order(abs(all$logFC),decreasing=TRUE),"ID"]
        DE <- all[1:desiredNum]
      }
    }
    DE <- named_FC[DE]
  }else{
    stop(paste0("ERROR: in GenerateDEGenes, unknown method ", method,".\n"))
  }

  if(length(DE)>0){
    message(paste0("\nINFO: ",length(DE)," genes are considered differentially expressed."))
    return(DE)
  }else{
      stop(paste0("ERROR: in GenerateDEGenes returned zero DE gene.\n"))
  }
}



generateKEGGInfo <- function(do.parallel=TRUE, mc.cores=1, path_to_xml=NULL,DataObject=NULL){
  pwys <- parseKEGGxml(path_to_xml = path_to_xml, mc.cores=mc.cores)
  kpg <- keggPathwayGraphs(pwys = pwys, mc.cores=mc.cores)
  KEGGgraph <-mergeGraphs(kpg, edgemode="directed" )
  nodes(KEGGgraph) <- gsub(paste0("^",DataObject$Organism,":"),"",nodes(KEGGgraph))

  named_FC <- as.numeric(DataObject$exprTable[,"logFC"])
  names(named_FC) <- as.character(DataObject$exprTable[,"ID"])
  named_FC <- named_FC[names(named_FC)%in%nodes(KEGGgraph)]
  named_PVAL <- as.numeric(DataObject$exprTable[,"adj.P.Val"])
  names(named_PVAL) <- as.character(DataObject$exprTable[,"ID"])
  named_PVAL <- named_PVAL[names(named_PVAL)%in%nodes(KEGGgraph)]

  return(list(pwys=pwys,kpg=kpg,KEGGgraph=KEGGgraph,named_FC=named_FC,named_PVAL=named_PVAL))
}



# find pathways with KO genes and return their IDs
getTargetPathways <- function(DataObject = NULL,pwys = NULL){
  # rownames of pwys has to be correctly set!
  if(is.null(DataObject)||is.null(pwys)){
    stop("ERROR: DataObject and/or KEGG pathways cannot be NULL.")
  }
  targetPathways <- NULL
  index <- sort(unique(unlist(lapply(DataObject$KOgeneEntrez,grep,pwys))))
  targetPathways <- unlist(lapply(pwys[index],getName))
  names(targetPathways) <- NULL
  targetPathways <- gsub("path:mmu","",targetPathways)
  return(targetPathways)
}



getUserInput <- function(){
  cat("\nINFO: Acquiring user parameters...\n")
  input <- NULL
  run_SPIA <- FALSE
  run_RONTOTOOLS_PE <- FALSE
  run_RONTOTOOLS_PDIS <- FALSE
  run_PADOG <- FALSE
  run_GSA <- FALSE
  run_SAFE <- FALSE
  run_GSEA <- FALSE

  # Preprocessed .RData
  cat("\nINFO: files in current working directory:\n")
  rdatalist <- list.files(path=getwd(),pattern="*PREP.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
  show(rdatalist)
  input <- readline(prompt = "\nEnter preprocessed data index: ")
  dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
  if(length(ls(dataEnv))>1){
    stop("\nERROR: loaded more than one R object. Exit.\n")
  }
  DataObject <- get(ls(dataEnv),dataEnv)

  # KEGG xmls
  cat("\nINFO: directories in current working directory:\n")
  rdatalist <- list.dirs(path=getwd(),full.names=TRUE,recursive=TRUE)
  rdatalist <- append(x=rdatalist, values=grep("KEGG",list.dirs(path=system.file("extdata",package="pathwayko"),full.names=TRUE,recursive=TRUE),value=TRUE))
  show(rdatalist)
  input <- readline(prompt = "\nEnter KEGG xml directory index: ")
  path_to_xml <- rdatalist[[as.integer(input)]]

  # mmuSPIA.RData
  cat("\nINFO: SPIA data in current working directory:\n")
  rdatalist <- list.files(path=getwd(),pattern="mmuSPIA.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
  rdatalist <- append(x=rdatalist, values=list.files(path=system.file("extdata",package="pathwayko"),pattern="mmuSPIA.RData",full.names=TRUE,recursive=TRUE, include.dirs=FALSE))
  show(rdatalist)
  input <- readline(prompt = "\nEnter SPIA data index: ")
  spia.data.dir <- rdatalist[[as.integer(input)]]
  spia.data.dir <- substr(spia.data.dir,1,nchar(spia.data.dir)-13)


  # Make dir for results
  result_output <- paste(DataObject$title,"_result_output",sep="")
  if(!dir.exists(result_output)){
    dir.create(result_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
  }
  result_output <- paste(result_output,"/",sep="")

  input <- readline(prompt = "\nRun SPIA?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_SPIA <- TRUE
  }

  input <- readline(prompt = "\nRun ROntoTools PE?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_RONTOTOOLS_PE <- TRUE
  }

  input <- readline(prompt = "\nRun ROntoTools PDIS?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_RONTOTOOLS_PDIS <- TRUE
  }

  input <- readline(prompt = "\nRun PADOG?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_PADOG <- TRUE
  }

  input <- readline(prompt = "\nRun GSA?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_GSA <- TRUE
  }

  input <- readline(prompt = "\nRun SAFE?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_SAFE <- TRUE
  }

  input <- readline(prompt = "\nRun GSEA?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_GSEA <- TRUE
  }

  DE.method <- NULL
  while(TRUE){
    input <- readline(prompt = "\nHow to choose DE genes? ('HES1','HES2','HES3' for HighEdgeS lower, change point, upper bound, 'CLA' for Classical)\nNote: geneset-based methods will not use DE genes regardlessly: ")
    DE.method <- input
    if(DE.method%in%c("HES1","HES2","HES3","CLA")){
      break
    }else{
      cat(paste0("\nINFO: invalid method choice ",DE.method,", please try again.\n"))
    }
  }
  
  DE.alpha <- NULL
  if(DE.method=="CLA"){
    while(TRUE){
    input <- readline(prompt = paste0("\nFraction of genes to be considered? (e.g. 0.05): "))
    DE.alpha <- as.numeric(input)
    if(DE.alpha<=0 || DE.alpha>=1.0){
      cat(paste0("\nINFO: invalid alpha value", DE.alpha, ", please try again.\n"))
    }else{
      cat(paste0("\nINFO: using alpha ",DE.alpha," for choosing DE genes.\n"))
      break
    }
    }
  }

  while(TRUE){
    input <- readline(prompt = "\nLower bound of specificity (X-axis) in partial ROC? (e.g. 90): ")
    tmp1 <- as.integer(input)

    input <- readline(prompt = "\nLower bound of sensitivity (Y-axis) in partial ROC? (e.g. 90): ")
    tmp2 <- as.integer(input)

    if((!is.na(tmp1))&&(!is.na(tmp2))&&(tmp1>0)&&(tmp1<100)&&(tmp2>0)&&(tmp2<100)){
      roc.partial.sp.lowerbound <- tmp1
      roc.partial.se.lowerbound <- tmp2
      break
    }else{
      cat("\nINFO: please enter valid lower bounds between 0 and 100 exclusive. Try again.\n")
    }
  }

  filename_prefix <- paste0(DataObject$title,"_",DataObject$normMethod,"_",DataObject$KOgeneSymbol,"_",DE.method,"_")

  return(list(DataObject=DataObject,
              result_output=result_output,
              filename_prefix=filename_prefix,
              run_SPIA=run_SPIA,
              run_RONTOTOOLS_PE=run_RONTOTOOLS_PE,
              run_RONTOTOOLS_PDIS=run_RONTOTOOLS_PDIS,
              run_PADOG=run_PADOG,
              run_GSA=run_GSA,
              run_SAFE=run_SAFE,
              run_GSEA=run_GSEA,
              DE.method=DE.method,
              DE.alpha=DE.alpha,
              roc.partial.sp.lowerbound=roc.partial.sp.lowerbound,
              roc.partial.se.lowerbound=roc.partial.se.lowerbound,
              path_to_xml=path_to_xml,
              spia.data.dir=spia.data.dir
              ))
}



getUserInputBatch <- function(){
  cat("\nINFO: Acquiring user parameters...\n")
  input <- NULL
  run_SPIA <- FALSE
  run_RONTOTOOLS_PE <- FALSE
  run_RONTOTOOLS_PDIS <- FALSE
  run_PADOG <- FALSE
  run_GSA <- FALSE
  run_SAFE <- FALSE
  run_GSEA <- FALSE

  DataFiles <- NULL
  # Preprocessed .RData
  while(TRUE){
    cat("\nINFO: files in current working directory:\n")
    rdatalist <- list.files(path=getwd(),pattern="*PREP.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
    show(rdatalist)
    input <- readline(prompt = "\nProcess all?: (Y/N) ")
    if(input=="Y"||input=="y"){
      DataFiles <- rdatalist
      cat(paste0("\nINFO: successfully queued ",length(DataFiles)," RData files.\n"))
      break
    }else{
      cat("\nINFO: Please make adjustments, and press any key to continue.\n")
      input <- readline(prompt="Continue?")
    }
  }

  # KEGG xmls
  cat("\nINFO: directories in current working directory:\n")
  rdatalist <- list.dirs(path=getwd(),full.names=TRUE,recursive=TRUE)
  rdatalist <- append(x=rdatalist, values=grep("KEGG",list.dirs(path=system.file("extdata",package="pathwayko"),full.names=TRUE,recursive=TRUE),value=TRUE))
  show(rdatalist)
  input <- readline(prompt = "\nEnter KEGG xml directory index: ")
  path_to_xml <- rdatalist[[as.integer(input)]]

  # mmuSPIA.RData
  cat("\nINFO: SPIA data in current working directory:\n")
  rdatalist <- list.files(path=getwd(),pattern="mmuSPIA.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
  rdatalist <- append(x=rdatalist, values=list.files(path=system.file("extdata",package="pathwayko"),pattern="mmuSPIA.RData",full.names=TRUE,recursive=TRUE, include.dirs=FALSE))
  show(rdatalist)
  input <- readline(prompt = "\nEnter SPIA data index: ")
  spia.data.dir <- rdatalist[[as.integer(input)]]
  spia.data.dir <- substr(spia.data.dir,1,nchar(spia.data.dir)-13)

  input <- readline(prompt = "\nRun SPIA?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_SPIA <- TRUE
  }
  input <- readline(prompt = "\nRun ROntoTools PE?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_RONTOTOOLS_PE <- TRUE
  }
    input <- readline(prompt = "\nRun ROntoTools PDIS?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_RONTOTOOLS_PDIS <- TRUE
  }

  input <- readline(prompt = "\nRun PADOG?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_PADOG <- TRUE
  }

  input <- readline(prompt = "\nRun GSA?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_GSA <- TRUE
  }

  input <- readline(prompt = "\nRun SAFE?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_SAFE <- TRUE
  }

  input <- readline(prompt = "\nRun GSEA?: (Y/N) ")
  if(input=="Y"||input=="y"){
    run_GSEA <- TRUE
  }

  DE.method <- NULL
  while(TRUE){
    input <- readline(prompt = "\nHow to choose DE genes? ('HES1','HES2','HES3' for HighEdgeS lower, change point, upper bound, 'CLA' for Classical)\nNote: geneset-based methods will not use DE genes regardlessly: ")
    DE.method <- input
    if(DE.method%in%c("HES1","HES2","HES3","CLA")){
      break
    }else{
      cat(paste0("\nINFO: invalid method choice ",DE.method,", please try again.\n"))
    }
  }
  
  DE.alpha <- NULL
  if(DE.method=="CLA"){
    while(TRUE){
      input <- readline(prompt = paste0("\nFraction of genes to be considered? (e.g. 0.05): "))
      DE.alpha <- as.numeric(input)
      if(DE.alpha<=0 || DE.alpha>=1.0){
        cat(paste0("\nINFO: invalid alpha value", DE.alpha, ", please try again.\n"))
      }else{
        cat(paste0("\nINFO: using alpha ",DE.alpha," for choosing DE genes.\n"))
        break
      }
    }
  }

  while(TRUE){
    input <- readline(prompt = "\nLower bound of specificity (X-axis) in partial ROC? (e.g. 90): ")
    tmp1 <- as.integer(input)

    input <- readline(prompt = "\nLower bound of sensitivity (Y-axis) in partial ROC? (e.g. 90): ")
    tmp2 <- as.integer(input)

    if((!is.na(tmp1))&&(!is.na(tmp2))&&(tmp1>0)&&(tmp1<100)&&(tmp2>0)&&(tmp2<100)){
      roc.partial.sp.lowerbound <- tmp1
      roc.partial.se.lowerbound <- tmp2
      break
    }else{
      cat("\nINFO: please enter valid lower bounds between 0 and 100 exclusive. Try again.\n")
    }
  }

  return(list(
              DataFiles=DataFiles,
              run_SPIA=run_SPIA,
              run_RONTOTOOLS_PE=run_RONTOTOOLS_PE,
              run_RONTOTOOLS_PDIS=run_RONTOTOOLS_PDIS,
              run_PADOG=run_PADOG,
              run_GSA=run_GSA,
              run_SAFE=run_SAFE,
              run_GSEA=run_GSEA,
              DE.method=DE.method,
              DE.alpha=DE.alpha,
              roc.partial.sp.lowerbound=roc.partial.sp.lowerbound,
              roc.partial.se.lowerbound=roc.partial.se.lowerbound,
              path_to_xml=path_to_xml,
              spia.data.dir=spia.data.dir
              ))
}



saveDEInfo <- function(DE=NULL,file=NULL){
  if(is.null(DE)){
    return(NULL)
  }else{
    suppressMessages(require(org.Mm.eg.db))
    DEInformation <- as.data.frame(suppressMessages(select(org.Mm.eg.db,names(DE),columns=c("ENTREZID","SYMBOL","GENENAME"),keyType="ENTREZID")))
    if(!is.null(file)){
      write.csv(DEInformation, file=file)
    }
    return(DEInformation)
  }
}



savePathwayInfo <- function(targetPathways=NULL,pwys=NULL,file=NULL){
  if(is.null(targetPathways)){
    return(NULL)
  }else if(is.null(pwys)){
    stop("\nERROR: missing arguments in savePathwayInfo.\n")
  }else{
    tmp <- lapply(pwys[paste0("mmu",targetPathways)],function(X){
       pInfo <- X@pathwayInfo
       data.frame(pInfo@number,pInfo@name,pInfo@org,pInfo@title,pInfo@link)
      })
    pathwayInformation <- NULL
    for(i in 1:length(tmp)){
      pathwayInformation <- rbind(pathwayInformation,tmp[[i]])
    }
    colnames(pathwayInformation) <- gsub("pInfo.","",colnames(pathwayInformation))
    if(!is.null(file)){
      write.csv(pathwayInformation,file=file)
    }
    return(pathwayInformation)
  }
}



savePathwayAnalysis <- function(collected=NULL,filePath=NULL){
  if(is.null(collected)||is.null(filePath)){
    return(NULL)
  }
  ret <- NULL
  if(!is.null(collected$RONTOTOOLS_PE_RES)){
    ret$RONTOTOOLS_PE <- collected[["RONTOTOOLS_PE_RES"]]
    write.csv(collected[["RONTOTOOLS_PE_RES"]],file=paste0(filePath,"RONTOTOOLS_PE.csv"))
  }
  if(!is.null(collected$RONTOTOOLS_PDIS_RES)){
    ret$RONTOTOOLS_PDIS <- collected[["RONTOTOOLS_PDIS_RES"]]
    write.csv(collected[["RONTOTOOLS_PDIS_RES"]],file=paste0(filePath,"RONTOTOOLS_PDIS.csv"))
  }
  if(!is.null(collected$SPIA_RES)){
    ret$SPIA <- collected[["SPIA_RES"]]
    write.csv(collected[["SPIA_RES"]],file=paste0(filePath,"SPIA.csv"))
  }
  if(!is.null(collected$PADOG_RES)){
    ret$PADOG <- collected[["PADOG_RES"]]
    write.csv(collected[["PADOG_RES"]],file=paste0(filePath,"PADOG.csv"))
  }
  if(!is.null(collected$GSA_RES)){
    ret$GSA <- collected[["GSA_RES"]]
    write.csv(collected[["GSA_RES"]],file=paste0(filePath,"GSA.csv"))
  }
  if(!is.null(collected$SAFE_RES)){
    ret$SAFE <- collected[["SAFE_RES"]]
    write.csv(collected[["SAFE_RES"]],file=paste0(filePath,"SAFE.csv"))
  }
  if(!is.null(collected$GSEA_RES)){
    ret$GSEA <- collected[["GSEA_RES"]]
    write.csv(collected[["GSEA_RES"]],file=paste0(filePath,"GSEA.csv"))
  }
  return(ret)
}
