#' @title rendering the true positive KO KEGG pathways by the up- and down-regulated DE genes
#' @usage pathwayview()
#' @description It renders the true positive KO KEGG pathways by the up- and down-regulated DE genes, 
#'  which contain and are impacted by the KO genes. 
#'	Reference: Luo, W. and Brouwer C., Pathview: an R/Bioconductor package for pathway-based data integration 
#' 	and visualization. Bioinformatics, 2013, 29(14): 1830-1831, doi: 10.1093/bioinformatics/btt285
#' @return TRUE on success, FALSE otherwise
#' @export

pathwayview <- function(){
	# pathview needs its internal database to function hence has to be attached and not just namespace loaded
	if(!suppressPackageStartupMessages(require(pathview))){
		cat("\nERROR: package 'pathview' is not available, install it and try again.\n")
		return(FALSE)
	}


	cat("\nINFO: summary files in current working directory:\n")
	rdatalist <- list.files(path=getwd(),pattern="*SUM.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter stats data index: ")
	dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
	if(length(ls(dataEnv))>1){
		cat("\nERROR: loaded more than one R object. Exit.\n")
		return(FALSE)
	}
	sumObject <- get(ls(dataEnv),dataEnv)
	dataTables <- sumObject$PAResults[grep("SPIA",names(sumObject$PAResults))]

	cat("\nINFO: preprocessed files in current directory:\n")
	rdatalist <- list.files(path=".",pattern="*PREP.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter preprocessed data index: ")
	dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
	if(length(ls(dataEnv))>1){
		cat("\nERROR: loaded more than one R object. Exit.\n")
		return(FALSE)
	}
	DataObject <- get(ls(dataEnv),dataEnv)

	cat("\nINFO: directories in current directory:\n")
	rdatalist <- list.dirs(path=".",full.names=TRUE,recursive=FALSE)
	rdatalist <- append(x=rdatalist, values=grep("KEGG",list.dirs(path=system.file("extdata",package="pathwayko"),full.names=TRUE,recursive=TRUE),value=TRUE))
	show(rdatalist)
	input <- readline(prompt = "\nEnter KEGG data index: ")
	kegg.dir <- rdatalist[as.numeric(input)]

	result_output <- paste0(DataObject$title,"_pathwayview_output")
	if(!dir.exists(result_output)){
		dir.create(result_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}

	cat("\nINFO: processing...\n")
	for(i in 1:length(dataTables)){
		dataTable <- dataTables[[i]]
		geneID <- dataTable[dataTable$Response==TRUE,"KEGGLINK"]
		geneID <- gsub("^.*show_pathway\\?","",geneID)
		geneID <- strsplit(geneID,"\\+")
		method <- names(dataTables[i])
		res <- lapply(geneID,function(X,a,b){
				tmpTable <- a$exprTable[X[2:length(X)],]
				gene.data <- tmpTable[,1]
				names(gene.data) <- rownames(tmpTable)
				pathway.id = X[1]
				species = substr(X[1],1,3)
				out.suffix <- paste0(DataObject$title,".",DataObject$KOgeneSymbol,".",b)
				pv.out <- tryCatch(
							{
								suppressMessages(
									pathview::pathview(
										gene.data = gene.data, 
										pathway.id = pathway.id, 
										species = species, 
										kegg.dir=kegg.dir, 
										out.suffix = out.suffix, 
										kegg.native = TRUE, 
										na.col = "gray", 
										new.signature=FALSE)
									)
							},
								error = function(cond) {
						            message(paste0("\nINFO: encountered error when processing: ", pathway.id))
						            message("\nINFO: Here's the original error message:") # find an exceptional case & handle it ?????
						            message(cond)
						            message("\n")
						            brokenFilePath <- list.files(path=".",pattern="*.png",full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
						            brokenFilePath <- grep(pathway.id,brokenFilePath,value=TRUE)
						            file.remove(brokenFilePath)
						            return(NULL)
					        	}
					        )
		},a=DataObject,b=method)

		rdatalist <- list.files(path=".",pattern="*.png",full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
		pngs <- unlist(lapply(rdatalist,function(X){
						tryCatch(image_read(X),
							error=function(cond){
								return(NULL)
							})
						}))
		for(i in 1:length(pngs)){
			image_write(pngs[[i]], format="pdf", path=gsub(".png$",".pdf",rdatalist[[i]]))
		}
		if(!dir.exists(paste0(result_output,"/pngs"))) {
			dir.create(paste0(result_output,"/pngs"), showWarnings = TRUE, recursive = FALSE, mode = "0755")
		}
		file.copy(rdatalist,paste0(result_output,"/pngs"))
		file.remove(rdatalist)

		rdatalist <- list.files(path=".",pattern="*.pdf",full.names=TRUE,recursive=FALSE,include.dirs=FALSE)
		if(!dir.exists(paste0(result_output,"/pdfs"))) {
			dir.create(paste0(result_output,"/pdfs"), showWarnings = TRUE, recursive = FALSE, mode = "0755")
		}
		file.copy(rdatalist,paste0(result_output,"/pdfs"))
		file.remove(rdatalist)
	}

	cat("\nINFO: all process completed.\n")
}
