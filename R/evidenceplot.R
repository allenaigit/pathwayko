#' @title generating two-dimensional evidence plot
#' @usage evidenceplot()
#' @description It generates two-dimensional evidence plot for methods based on or derived from SPIA 
#' @details Specifically, this function looks for two columns in raw pathway analysis results:
#'		"pPERT" and "pNDE", to plot the evidence plot. 
#' 		Reference: Tarca A. L. et al. (2009) A novel signaling pathway impact analysis. Bioinformatics, 25, 75–82.)
#' @return NULL
#' @export

evidenceplot <- function(){

	cat("\nINFO: summary files in current working directory:\n")
	rdatalist <- list.files(path=getwd(),pattern="*SUM.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter summary data index: ")
	dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
	if(length(ls(dataEnv))>1){
		stop("\nERROR: loaded more than one R object. Exit.\n")
	}
	DataObject <- get(ls(dataEnv),dataEnv)

	output.dir <- gsub("result","EvidencePlot",DataObject$result_output)
	if(!dir.exists(output.dir)){
		dir.create(output.dir, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}

	cat("\nINFO: choose method:\n")
	rdatalist <- names(DataObject$PAResults)[grep("SPIA",names(DataObject$PAResults))]
	show(rdatalist)
	input <- readline(prompt = "\nEnter method index: ")
	method <- rdatalist[[as.numeric(input)]]
	data <- DataObject$PAResults[[method]]

	data[is.na(data[,"pPERT"]),"pPERT"] <- 1.0
	data[is.na(data[,"pNDE"]),"pNDE"] <- 1.0

	pdf(paste0(output.dir, method, ".evidence_plot.pdf"))
	plotP <- plotP(data,threshold=0.05)
	points(I(-log(pPERT))~I(-log(pNDE)),data=data[data$ID=="05210",],col="green",pch=19,cex=1.5)
	plotP + geom_point(alpha=0.01)
	dev.off()

	cat("\nINFO: completed.\n")
}