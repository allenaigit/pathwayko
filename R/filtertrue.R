#' @title extracting the list of true positive KO KEGG pathways
#' @usage filtertrue()
#' @description It extracts the list of true positive KO KEGG pathways that contain and are impacted by the KO gene from the results of each method
#' @details This function scans user provided directory for any .csv files under the current directory
#'		recursively and reads in all .csv files not tagged with control keywords like "AUC" or "pathway".
#'		Hence this function should only be used on directories with no .csv files from other sources.
#'		After data input, this function filters out from all PA results entries of true positive KEGG pathways
#' 		in that analysis and saves each of filtered pathway analysis result in a separate .csv file.
#' @return NULL
#' @export

filtertrue <- function(){
	cat("\nINFO: directories in current directory:\n")
	rdatalist <- list.dirs(path=".",full.names=TRUE,recursive=FALSE)
	rdatalist <- grep("_result_output",rdatalist, value=TRUE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter data index: ")
	data.dir <- rdatalist[as.numeric(input)]

	cat(paste0("\nINFO: csv files processed in ",data.dir," :\n"))
	rdatalist <- list.files(path=data.dir,pattern="*.csv",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	rdatalist <- rdatalist[!grepl("AUC",rdatalist)]
	rdatalist <- rdatalist[!grepl("genes",rdatalist)]
	rdatalist <- rdatalist[!grepl("stats",rdatalist)]
	rdatalist <- rdatalist[!grepl("pathway",rdatalist)]
	rdatalist <- rdatalist[!grepl("PostProc",rdatalist)]
	show(rdatalist)
	cat("\n")


	lapply(rdatalist,function(X,t){
			fileName <- unlist(strsplit(X,"/"))
			fileName <- gsub(".csv$",".PostProc.TRUE.csv",fileName[[length(fileName)]])
			cat(paste0("INFO: generating ", fileName,"...\n"))
			dataTable <- tryCatch(
					{
						read.csv(X,row.names=1)
					},
					error = function(cond) {
						tmp <- read.csv(X)
						rownames(tmp) <- 1:nrow(tmp)
						tmp$X <- NULL
			            return(tmp)
		        	}
				)
			responseColName <- colnames(dataTable)[grep("response",colnames(dataTable),ignore.case=TRUE)]
			dataTable <- dataTable[dataTable[,responseColName]==TRUE,]
			rankColName <- colnames(dataTable)[grep("rank",colnames(dataTable),ignore.case=TRUE)]
			if(nrow(dataTable)>0){
				dataTable$RankInverse <- 1/as.numeric(dataTable[,rankColName])
			}
			write.csv(dataTable,file=paste0(t,"/",fileName))
		},t=data.dir)

	cat("\nINFO: completed.\n")
}
