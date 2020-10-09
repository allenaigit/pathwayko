#' @title performing the KO pathway enrichment analysis on single dataset
#' @usage pathwayko()
#' @description It invokes utilities, highedges, functions and methods to calculate the scores of all edges (Hanoudi et al., 2017) in a known global KEGG graph (Zhang and Wiemann, 2009) to yield the distribution of edge scores, and thus automatically determines a change point on the distribution of edge scores by the change-point analysis method (Killick and Eckley, 2014). Outputs include the KO gene-associated sub-graph, the list of differential expression (DE) genes, and the list of true positive KO KEGG pathways.
#' 		References: Hanoudi S., et al. (2017) Identifying biologically relevant putative mechanisms in a given phenotype comparison. PLoS ONE 12, e0176950; Zhang J. D. and Wiemann S. (2009) KEGGgraph: a graph approach to KEGG pathway in R and bioconductor. Bioinformatics, 25, 1470–1471; Killick R. and Eckley I. (2014) changepoint: An R package for changepoint analysis. J. Statistic Software. 58,1–19.
#' @param batch boolean, should only be set to TRUE by pathwayko_batch function
#' @param batch_DataObject should only be set by pathwayko_batch function
#' @param batch_uip should only be set by pathwayko_batch function
#' @return NULL if encountered error, TRUE when batch is false and an internal 
#' 	object with results if batch is true
#' @export

pathwayko <- function(batch=FALSE,batch_DataObject=NULL,batch_uip=NULL){
	cat("\nINFO: starting pathwayko...\n")

	############################################
	# 	Internal Parameters:                   # 
	#	Do not tamper with these unless        #
	#   you know what you're doing             #
	############################################

	maxCores <- detectCores(logical=TRUE)
	do.parallel <- TRUE
	mc.cores <- 1
	if(maxCores<2){
		do.parallel <- FALSE
	}else{
		# safe estimate of 2GB per core
		mc.cores <- max(1,min(as.integer(getFreeMemoryMB()/2048)-1,maxCores))
	}

	DataObject <- NULL
	ResultObject <- NULL
	path_to_xml <- NULL
	spia.data.dir <- NULL
	filename_prefix <- NULL
	result_output <- NULL
	uip <- NULL

	pwys <- NULL
	input <- NULL

	DE.method <- NULL
	DE.alpha <- NULL

	roc.partial.se.lowerbound <- 0
	roc.partial.sp.lowerbound <- 0

	do.methods <- list(
						run_SPIA = FALSE,
						run_RONTOTOOLS_PE = FALSE,
						run_RONTOTOOLS_PDIS = FALSE,
						run_PADOG = FALSE,
						run_GSA = FALSE,
						run_SAFE = FALSE,
						run_GSEA = FALSE
						)




	########################################
	############## User Input ##############
	########################################
	if(batch){
		if(is.null(batch_DataObject)||is.null(batch_uip)){
			stop("\nERROR: missing arguments.\n")
		}
		uip <- batch_uip
		DataObject <- batch_DataObject
	}else{
		uip <- getUserInput()
		DataObject <- uip$DataObject
	}
	
	result_output <- uip$result_output
	filename_prefix <- uip$filename_prefix
	do.methods$run_SPIA <- uip$run_SPIA
	do.methods$run_RONTOTOOLS_PE <- uip$run_RONTOTOOLS_PE
	do.methods$run_RONTOTOOLS_PDIS <- uip$run_RONTOTOOLS_PDIS
	do.methods$run_PADOG <- uip$run_PADOG
	do.methods$run_GSA <- uip$run_GSA
	do.methods$run_SAFE <- uip$run_SAFE
	do.methods$run_GSEA <- uip$run_GSEA
	DE.method <- uip$DE.method
	DE.alpha <- uip$DE.alpha
	roc.partial.sp.lowerbound <- uip$roc.partial.sp.lowerbound
	roc.partial.se.lowerbound <- uip$roc.partial.se.lowerbound
	path_to_xml <- uip$path_to_xml
	spia.data.dir <- uip$spia.data.dir
	ResultObject$result_output <- result_output
	ResultObject$filename_prefix <- filename_prefix
	ResultObject$methods <- do.methods



	########################################
	########## KEGG and DE genes ###########
	########################################
	cat("\nINFO: preprocessing...\n")

	# Generate KEGG objects from local version of KEGG xmls
	ret <- generateKEGGInfo(do.parallel=do.parallel, mc.cores=mc.cores, path_to_xml=path_to_xml,DataObject=DataObject)
	pwys <- ret$pwys
	kpg <- ret$kpg
	KEGGgraph <- ret$KEGGgraph
	named_FC <- ret$named_FC
	named_PVAL <- ret$named_PVAL

	# Obtain DE genes
	DE <- GenerateDEGenes(method=DE.method, KEGGgraph=KEGGgraph, DataObject=DataObject, named_FC=named_FC, named_PVAL=named_PVAL, alpha=DE.alpha, result_output=result_output, filename_prefix=filename_prefix)

	# Obtain KO target Pathways
	targetPathways <- getTargetPathways(DataObject,pwys)
	cat("\nINFO: target pathways:\n")
	show(targetPathways)



	########################################
	########### Pathway Analysis ###########
	########################################
	
	collected <- .apply_parallel_pathway_analysis(do.methods=do.methods,DE=DE,kpg=kpg,pwys=pwys,DataObject=DataObject,targetPathways=targetPathways,named_FC=named_FC,named_PVAL=named_PVAL,spia.data.dir=spia.data.dir,path_to_xml=path_to_xml)
	if(is.null(collected)){
		cat("\nWARNING: no pathway analysis performed.\n")
		return(NULL)
	}
	
	########################################
	########### Save PA Results ############
	########################################
	ResultObject$DE <- saveDEInfo(DE=DE,file=paste0(result_output,filename_prefix,"DE_genes.csv"))
	ResultObject$PathwayInfo <- savePathwayInfo(targetPathways=targetPathways,pwys=pwys,file=paste0(result_output,filename_prefix,"target_pathways.csv"))
	ResultObject$PAResults <- savePathwayAnalysis(collected=collected,filePath=paste0(result_output,filename_prefix))


	########################################
	################# ROC ##################
	########################################
	ResultObject$ROC <- suppressWarnings(generateROC(collected=collected,result_output=result_output,filename_prefix=filename_prefix,targetPathways=targetPathways,roc.partial.se.lowerbound=roc.partial.se.lowerbound,roc.partial.sp.lowerbound=roc.partial.sp.lowerbound))
	save(ResultObject, file=paste0(result_output,filename_prefix,"SUM.RData"), compress="xz")

	cat("\nINFO: pathwayko completed.\n")
	if(batch){
		return(ResultObject)
	}else{
		return(TRUE)
	}
}



#' @title batch-performing the KO pathway enrichment analysis on multiple datasets
#' @usage pathwayko_batch()
#' @description It conducts the KO pathway enrichment analysis on multiple datasets in a pipeline manner.
#' @details This function scans all directory under the working directory to identify
#'	data object generated by preprocess function and apply pathwayko to all said objects
#' 	using same sets of parameters obtained from the user
#' @return TRUE if all process completed successfully, FALSE otherwise
#' @export

pathwayko_batch <- function(){
	uip <- getUserInputBatch()
	DataFiles <- uip$DataFiles
	uip$DataFiles <- NULL
	for(i in 1:length(DataFiles)){
		cat(paste("\nINFO: batch job (",i,"/",length(DataFiles),")...\n"))
		dataEnv <- loadToEnvironment(DataFiles[[i]])
		if(length(ls(dataEnv))>1){
			stop("\nERROR: loaded more than one R object. Exit.\n")
		}
		DataObject <- get(ls(dataEnv),dataEnv)
		temp_uip <- uip
		result_output <- paste(DataObject$title,"_result_output",sep="")
		if(!dir.exists(result_output)){
			dir.create(result_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
		}
		result_output <- paste(result_output,"/",sep="")
		filename_prefix <- paste0(DataObject$title,"_",DataObject$normMethod,"_",DataObject$KOgeneSymbol,"_",temp_uip$DE.method,"_")
		temp_uip$result_output <- result_output
		temp_uip$filename_prefix <- filename_prefix

		res <- pathwayko(batch=TRUE,batch_DataObject=DataObject,batch_uip=temp_uip)
		if(is.null(res)){
			cat(paste0("\nERROR: during ",DataObject$title,".\n"))
			return(FALSE)
		}
		sink("pathwayko_batch_log.txt",append=TRUE)
		cat("\n==============================\n")

		cat(paste0("\nINFO: ",res$filename_prefix," (",i,",",length(DataFiles),")\n"))
		cat(paste0("\nINFO: ", length(res$DE)," differentially expressed genes.\n"))
		cat("\nINFO: target pathways:\n")
		show(res$PathwayInfo)
		cat("\n==============================\n\n")
		sink()
		cat(paste("INFO: batch job (",i,"/",length(DataFiles),") completed.\n\n\n"))
	}
	return(TRUE)
}



#' @title for pathwayko_demo on 10 KO GEO datasets
#' @usage pathwayko_demo()
#' @description It generates demo results for testing on the 10 KO GEO datasets
#' @return TRUE if all process completed successfully, FALSE otherwise
#' @export

pathwayko_demo <- function(){
	uip <- list(
					DataFiles=list.files(path=system.file("extdata/geodata",package="pathwayko"),pattern="*PREP.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE),
					run_SPIA=TRUE,
					run_RONTOTOOLS_PE=TRUE,
					run_RONTOTOOLS_PDIS=TRUE,
					run_PADOG=TRUE,
					run_GSA=TRUE,
					run_SAFE=TRUE,
					run_GSEA=TRUE,
					DE.method="HES1",
					DE.alpha=NULL,
					roc.partial.sp.lowerbound=90,
					roc.partial.se.lowerbound=90,
					path_to_xml=system.file("extdata/mmuKEGGxml",package="pathwayko"),
					spia.data.dir=system.file("extdata/",package="pathwayko")
				)

	DataFiles <- uip$DataFiles
	uip$DataFiles <- NULL
	for(i in 1:length(DataFiles)){
		cat(paste("\nINFO: batch job (",i,"/",length(DataFiles),")...\n"))
		dataEnv <- loadToEnvironment(DataFiles[[i]])
		if(length(ls(dataEnv))>1){
			stop("\nERROR: loaded more than one R object. Exit.\n")
		}
		DataObject <- get(ls(dataEnv),dataEnv)
		temp_uip <- uip
		result_output <- paste(DataObject$title,"_result_output",sep="")
		if(!dir.exists(result_output)){
			dir.create(result_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
		}
		result_output <- paste(result_output,"/",sep="")
		filename_prefix <- paste0(DataObject$title,"_",DataObject$normMethod,"_",DataObject$KOgeneSymbol,"_",temp_uip$DE.method,"_")
		temp_uip$result_output <- result_output
		temp_uip$filename_prefix <- filename_prefix

		res <- pathwayko(batch=TRUE,batch_DataObject=DataObject,batch_uip=temp_uip)
		if(is.null(res)){
			cat(paste0("\nERROR: during ",DataObject$title,".\n"))
			return(FALSE)
		}
		sink("pathwayko_batch_log.txt",append=TRUE)
		cat("\n==============================\n")

		cat(paste0("\nINFO: ",res$filename_prefix," (",i,",",length(DataFiles),")\n"))
		cat(paste0("\nINFO: ", length(res$DE)," differentially expressed genes.\n"))
		cat("\nINFO: target pathways:\n")
		show(res$PathwayInfo)
		cat("\n==============================\n\n")
		sink()
		cat(paste("INFO: batch job (",i,"/",length(DataFiles),") completed.\n\n\n"))
	}
	return(TRUE)
}
