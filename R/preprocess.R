#' @title preprocessing GSE data from CEL files and series matrix file
#' @usage preprocess()
#' @description This function reads in relevant CEL files (from *_RAW.tar)
#' 	file and series matrix file (*_series_matrix.txt.gz) typically obtained
#'	from NCBI database for a given GSE. The preprocessing of these data
#'	relies on the user to interactively provide necessary information
#' 	as such information is difficult to obtain prior to execution and hence 
#'	to automate. Reference: Carvalho B.S. and Irizarry R.A. (2010) A framework for 
#'  oligonucleotide microarray preprocessing. Bioinformatics, 16, 2363–2367; 
#'  Ritchie M.E., et al. (2015) limma powers differential expression analyses for 
#'  RNA-sequencing and microarray studies. Nucleic Acids Res, 43, e47.
#' @details To begin, the user is expected to have obtained '*_RAW.tar' and
#' 	'*_series_matrix.txt.gz' files of a given GSE. The user will first be
#'	asked to choose from files under the working directory with matching 
#'	suffix as input files. Then the user will be asked to provide relevant
#'	information such as keywords for case/control and KO gene names to 
#'	build experiment design for packages like 'oligo' and 'limma'. The
#'	final resault will be saved as an compressed R object to be read in
#'	and used by other parts of the package.
#'
#' @return TRUE when successful, FALSE otherwise
#' @export

preprocess <- function(){
	cat("\nINFO: loading libraries...\n")
	# Variables
	preprocess_output <- NULL
	input <- NULL
	dbPkg <- NULL
	pdPkg <- NULL
	KOgeneSymbol <- NULL
	KOgeneEntrez <- NULL
	GSEName <- NULL
	manualAnnotation <- FALSE
	rawData <- NULL


	########################################
	###########  Series Matrix  ############
	########################################
	cat("\nINFO: series matrix files in current directory:\n")
	matrixPath <- list.files(path=".",pattern=".*series_matrix.*",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(matrixPath)
	input <- readline(prompt = "\nEnter data index: ")
	matrixPath <- matrixPath[as.numeric(input)]

	cat("\nINFO: processing series matrix files...\n")

	# load namespace of GEOquery
	if(!suppressPackageStartupMessages(requireNamespace("GEOquery",quietly=TRUE))){
		cat("\nERROR: package 'GEOquery' is not available, install it and try again.\n")
		return(FALSE)
	}

	series_matrix <- suppressMessages(GEOquery::getGEO(filename=matrixPath,getGPL=FALSE))
	platform <- annotation(series_matrix)
	tmp <- pData(series_matrix)
	organism <- as.character(tmp[,grep("organism",colnames(tmp))][[1]])
	pDat <- as.data.frame(tmp[,c(grep("title",colnames(tmp)),grep("organism",colnames(tmp)),grep("platform",colnames(tmp)),grep("source",colnames(tmp)))],tmp[,grep("geo_accession",colnames(tmp))])
	colnames(pDat) <- c("title","organism","platform","source")
	cat("\nINFO: series matrix processed.\n")



	########################################
	######## Set Annotation Package  #######
	########################################
	if(platform=="GPL1261"){
		pdPkg <- "pd.mouse4302.mm.entrezg"
		dbPkg <- "mouse4302mmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL6246"){
		pdPkg <- "pd.mogene10st.mm.entrezg"
		dbPkg <- "mogene10stmmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL11180"){
		pdPkg <- "pd.htmg430pm.mm.entrezg"
		dbPkg <- "htmg430pmmmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL23038"){
		pdPkg <- "pd.clariomsmouse.mm.entrezg"
		dbPkg <- "clariomsmousemmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL81"){
		pdPkg <- "pd.mgu74av2.mm.entrezg"
		dbPkg <- "mgu74av2mmentrezg.db"
		manualAnnotation <- TRUE
	}else if(
				platform=="GPL16570"||
				platform=="GPL17791"||
				platform=="GPL23092"||
				platform=="GPL20710"
			){
		pdPkg <- "pd.mogene20st.mm.entrezg"
		dbPkg <- "mogene20stmmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL17400"){
		pdPkg <- "pd.mogene21st.mm.entrezg"
		dbPkg <- "mogene21stmmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL11533"){
		pdPkg <- "pd.mogene11st.mm.entrezg"
		dbPkg <- "mogene11stmmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL8321"){
		pdPkg <- "pd.mouse430a2.mm.entrezg"
		dbPkg <- "mouse430a2mmentrezg.db"
		manualAnnotation <- TRUE
	}else if(platform=="GPL339"){
		pdPkg <- "pd.moe430a.mm.entrezg"
		dbPkg <- "moe430ammentrezg.db"
		manualAnnotation <- TRUE
	}else if(
				platform=="GPL23535"||
				platform=="GPL21341"||
				platform=="GPL20775"
			){
		pdPkg <- "pd.mta10.mm.entrezg"
		dbPkg <- "mta10mmentrezg.db"
		manualAnnotation <- TRUE
	}else{
		cat(paste0("\n\n\n\n\nWARNING: unknown platform ",platform,".\n\n\n\n\n"))
	}
	if(manualAnnotation){
		if(!(suppressMessages(require(dbPkg,character.only=TRUE))&&suppressMessages(require(pdPkg,character.only=TRUE)))){
			cat(paste0("\nERROR: database ", dbPkg," and/or ", pdPkg, "cannot be loaded. Check and try again.\n"))
			return(FALSE)
		}
	}



	########################################
	############# CEL process  #############
	########################################
	cat("\nINFO: CEL archives in current directory:\n")
	celPath <- list.files(path=".",pattern=".*RAW.tar$",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(celPath)
	input <- readline(prompt = "\nEnter data index: ")
	celPath <- celPath[as.numeric(input)]

	input = readline(prompt = "\nEnter an user defined identifier: (e.g. GSE22873_MKO): ")
	GSEName <- input

	cat("\nINFO: processing CEL files...\n")
	# Make dir for results
	preprocess_output <- paste(GSEName,"_preprocess_output",sep="")
	if(!dir.exists(preprocess_output)){
		dir.create(preprocess_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}
	preprocess_output <- paste(preprocess_output,"/",sep="")

	CEL_list <- untar(celPath, list=TRUE)
	CEL_list <- CEL_list[grepl(".CEL.gz",CEL_list,ignore.case=TRUE)]
	# Limit CEL_list to GSM entries in pData from series matrix
	GSMName <- gsub(".CEL.gz$","",CEL_list)
	GSMName <- gsub(".cel.gz$","",GSMName)
	GSMName <- gsub("_.*$","",GSMName)
	CEL_list <- CEL_list[GSMName%in%rownames(pDat)]
	untar(celPath,files=CEL_list)

	if(manualAnnotation){
		# Read in CEL with specified pdInfo
		rawData <- oligo::read.celfiles(filenames=CEL_list,pkgname=pdPkg,verbose=FALSE)
	}else{
		rawData <- oligo::read.celfiles(filenames=CEL_list,verbose=FALSE)
	}

	file.remove(CEL_list)
	sampleNames(rawData) <- gsub(".CEL.gz","",sampleNames(rawData),ignore.case=TRUE)
	sampleNames(rawData) <- gsub("_.*$","",sampleNames(rawData),ignore.case=TRUE)
	# trim pData to GSM entries in rawData
	pDat <- pDat[rownames(pDat)%in%sampleNames(rawData),]
	write.csv(pDat, file=paste(preprocess_output,GSEName,"_pData.csv",sep=""))

	cat("\nINFO: all CEL files processed.\n")



	########################################
	#######  Case/Control & KO Gene ########
	########################################
	cat("\nINFO: processing experiment designs...\n")
	updatePData <- TRUE
	while(updatePData){
		pDat <- read.csv(file=paste(preprocess_output,GSEName,"_pData.csv",sep=""),row.names=1)
		cat("\nINFO: Reference pData:\n")
		show(pDat)
		input = readline(prompt = "\nUse this pData? (Y/N): ")
		if(grepl("Y",input,ignore.case=TRUE)){
			updatePData <- FALSE
		}else{
			cat(paste("\nINFO: please update pData file in ",preprocess_output,GSEName,"_pData.csv.\n",sep=""))
			input = readline(prompt = "\nEnter anything to refresh pData: ")
		}
	}

	case <- NULL
	control <- NULL
	while(TRUE){
		show(pDat)
		input = readline(prompt = "\nEnter case keyword for pattern matching: (e.g. \"MKO\"): ")
		CaseKey <- input
		input = readline(prompt = "\nEnter control keyword for pattern matching: (e.g. \"WT\"): ")
		ControlKey <- input
		input = readline(prompt = "\nEnter which column to match: title(1), source(2): ")
		if(grepl("1",input)){
			case <- rownames(pDat)[grepl(CaseKey,pDat$title,ignore.case=TRUE)]
			control <- rownames(pDat)[grepl(ControlKey,pDat$title,ignore.case=TRUE)]
		}else{
			case <- rownames(pDat)[grepl(CaseKey,pDat$source,ignore.case=TRUE)]
			control <- rownames(pDat)[grepl(ControlKey,pDat$source,ignore.case=TRUE)]
		}
		cat("\nINFO: selected case samples:\n")
		show(case)
		cat("\nINFO: selected control samples:\n")
		show(control)
		if(any(duplicated(c(case,control)))){
			cat("\nERROR: overlapping samples detected across two groups. Try again\n")
			next
		}else{	
			input = readline(prompt = "\nContinue? (Y/N): ")
			if(grepl("Y",input,ignore.case=TRUE)){
				break
			}else{
				next
			}
		}
	}

	KOgeneSymbol <- NULL
	KOgeneEntrez <- vector()
	while(TRUE){
		input = readline(prompt = "\nEnter KO gene symbol:\n(This is case sensitive and underscore delimited,\n e.g. Myd88_Ager): ")
		KOgeneSymbol <- input
		tempSymbol <- unlist(strsplit(KOgeneSymbol,"_"))
		mm <- get("org.Mm.eg.db")
		queryRes <- tryCatch(
						{
							suppressMessages(
								select(mm,keys=tempSymbol,columns=c("SYMBOL","ENTREZID","GENENAME"),keytype="SYMBOL")
							)
						},
							error = function(cond) {
					            message(paste0("\nWARNING: No gene matched with keyword ", tempSymbol,". Try again.\n"))
					            return(NULL)
				        	}
				        )

		if(is.null(queryRes)){
			next
		}
		cat("\nINFO: Matched genes:\n")
		show(queryRes)
		input = readline(prompt = "\nContinue? (Y/N): ")
		if(grepl("Y",input,ignore.case=TRUE)){
			for(i in 1:length(tempSymbol)){
				KOgeneEntrez[i] <- queryRes[[i,2]]
			}
			cat("\nINFO: matched EntrezID(s):\n")
			show(KOgeneEntrez)
			break
		}
	}
	cat("\nINFO: experiment design processed.\n")
		


	########################################
	########### Quality Control ############
	########################################
	# cat("\nINFO: Generating quality control plots...\n")

	# Removed since Oligo's own PLM method is full of bugs and poorly documented

	# cat("\nINFO: QC plots generated.\n")



	########################################
	##### Correction & Normalization #######
	########################################
	# RMA, quantile, median-polish
	cat("\nINFO: Performing RMA...\n\n")
	rmaData <- oligo::rma(rawData,background=TRUE,normalize=TRUE,subset=NULL)
	cat("\nINFO: RMA done...\n")



	########################################
	######### DataObject Generation ########
	########################################
	makeDataObject <- function(GSEName=NULL,normalized_eset=NULL,pDat=NULL,KOgeneSymbol=NULL,
								KOgeneEntrez=NULL,case=NULL,control=NULL,dbPkg=NULL,
								organism=NULL,normMethod=NULL)
	{
		if(is.null(GSEName)||is.null(normalized_eset)||is.null(pDat)||is.null(KOgeneSymbol)||
			is.null(KOgeneEntrez)||is.null(case)||is.null(control)||is.null(dbPkg)||
			is.null(organism)||is.null(normMethod))
		{
			stop("\nERROR: missing argument in makeDataObject(). Stoping.\n")
		}
		########################################
		################ Limma #################
		########################################
		# Trim Samples from eset
		eset <- normalized_eset[,  c(case,control)]
		# Trim pData
		local_pDat <- pDat[sampleNames(eset),]
		pData(eset) <- local_pDat
		# Clean unmappable probesets
		eset <- eset[!is.na(getEG(rownames(eset),dbPkg)),]
		# Convert from probeID to gene EntrezID
		rownames(eset) <- getEG(rownames(eset),dbPkg)


		# Generate design
		f1 <- factor(rep("case",nrow(local_pDat)), levels=c("case","control"))
		eset$description <- f1
		design <- model.matrix(~ description + 0, eset)
		colnames(design) <- levels(f1)
		design[row.names(design)%in%case, "case"] <- 1
		design[row.names(design)%in%case, "control"] <- 0
		design[row.names(design)%in%control, "control"] <- 1
		design[row.names(design)%in%control, "case"] <- 0
		# Correctly set description in eset
		for(i in 1:length(sampleNames(eset))){
			if(sampleNames(eset)[i]%in%control){
				eset$description[i] <- "control"
			}else{
				eset$description[i] <- "case"
			}
		}
		contrast.matrix <- makeContrasts(case-control, levels = design)
		# Linear fit, contrast fit, Bayes test
		fit <- lmFit(eset,design)
		fit1 <- contrasts.fit(fit,contrast.matrix)
		fit2 <- eBayes(fit1)
		tT <- topTable(fit2, adjust.method = "fdr", sort.by = "B", number=nrow(fit2))


		########################################
		######## Assemble DataObject ###########
		########################################
		DataObject <- NULL
		tT$ID <- rownames(tT)
		mm <- get(dbPkg)
		tT$symbols <- suppressMessages(select(mm,keys=tT$ID,columns="SYMBOL",keytype="ENTREZID")[,"SYMBOL"])



		DataObject$eset <- eset
		DataObject$pData <- pData(DataObject$eset)
		DataObject$exprTable <- tT
		DataObject$KOgeneSymbol <- KOgeneSymbol # gene symbol
		DataObject$KOgeneEntrez <- KOgeneEntrez # gene ID
		DataObject$normMethod <- normMethod
		if(any(grepl("Mus musculus",organism,ignore.case=TRUE))){
			DataObject$Organism <- "mmu"
		}else if(any(grepl("Homo Sapiens",organism,ignore.case=TRUE))){
			DataObject$Organism <- "hsa"
		}else{
			DataObject$Organism <- "Unknown"
		}
		DataObject$db <- dbPkg
		DataObject$case <- case
		DataObject$control <- control
		DataObject$title<-GSEName
		return(DataObject)
	}



	cat("\nINFO: Saving data object...\n")
	rmaData <- makeDataObject(GSEName=GSEName,normalized_eset=rmaData,pDat=pDat,
								KOgeneSymbol=KOgeneSymbol,KOgeneEntrez=KOgeneEntrez,
								case=case,control=control,dbPkg=dbPkg,
								organism=organism,normMethod="RMA")
	fileToSave <- paste(preprocess_output,rmaData$title,"_",rmaData$KOgeneSymbol,"_RMA_PREP.RData",sep="")
	save(rmaData, file=fileToSave, compress="xz")
	cat(paste0("\nINFO: Data object saved to ",fileToSave,".\n"))

	cat("\n\nINFO: preprocess successfully completed. Exiting...\n\n")

	return(TRUE)
}
