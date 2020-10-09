#' @title conducting ROC-curve based statistics analysis
#' @usage roctest()
#' @description It performs roctest on two ROC curves
#' @details Reference: Robin X. et al. (2011) pROC: an open-source package for R and S+ to analyze and 
#'  compare ROC curves. BMC Bioinformatics, 12, 77.
#' @return NULL
#' @export

roctest <- function(){

	cat("\nINFO: summary files in current working directory:\n")
	rdatalist <- list.files(path=getwd(),pattern="*SUM.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter summary data index: ")
	dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
	if(length(ls(dataEnv))>1){
		stop("\nERROR: loaded more than one R object. Exit.\n")
	}
	DataObject <- get(ls(dataEnv),dataEnv)
	output.dir <- gsub("result","roctest",DataObject$result_output)
	if(!dir.exists(output.dir)){
		dir.create(output.dir, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}


	cat("\nINFO: methods run:\n")
	methods <- names(DataObject$PAResults)
	show(methods)
	input <- readline(prompt = "\nEnter index of method one: ")
	indexOne <- as.numeric(input)
	cat(paste0("\nINFO: first method selected: ",methods[[indexOne]],"\n"))
	input <- readline(prompt = "\nEnter index of method two: ")
	indexTwo <- as.numeric(input)
	cat(paste0("\nINFO: second method selected: ",methods[[indexTwo]],"\n"))


	data1 <- DataObject$PAResults[[methods[[indexOne]]]]
	data2 <- DataObject$PAResults[[methods[[indexTwo]]]]
	legendOne <- paste0(gsub("_.*$","",DataObject$filename_prefix),"_",methods[[indexOne]])
	legendTwo <- paste0(gsub("_.*$","",DataObject$filename_prefix),"_",methods[[indexTwo]])
	legendText <- c(legendOne,legendTwo)
	boot.n <- 10000
	colors <- c("#1c61b6","#008600")

	cat(paste0("\nINFO: processing...\n"))



	############# 1-Test on the ROC curves themselves ###############
	pdf(paste0(output.dir,"ROC_test.venkatraman.ROCcurves.pdf"), height=5,width=5)

	rocobj1 <- plot.roc(data1$Response, data1$Predictor, 
		main="Statistical hypothesis testing on ROC curves",
		direction=">", quiet=TRUE, percent=TRUE,col=colors[[1]], grid=TRUE)

	rocobj2 <- plot.roc(data2$Response, data2$Predictor, direction=">", 
		add=TRUE, quiet=TRUE, percent=TRUE,col=colors[[2]], grid=FALSE)

	# Test on the ROC curves themselves
	testobj1 <- roc.test(rocobj1, rocobj2, method="venkatraman", boot.n=boot.n)
	text(50, 50, labels = paste("      ROC curves\n     venkatraman\n p-value =", format.pval(testobj1$p.value)), adj=c(0, .5))
	legend("bottomright", legend=legendText, col=colors, lwd=2)
	dev.off()



	############# 2-Test on the whole AUC with Bootstrap method ###############
	pdf(paste0(output.dir,"ROC_test.bootstrap.AUC.pdf"), height=5,width=5)

	rocobj1 <- plot.roc(data1$Response, data1$Predictor,
		main="Statistical hypothesis testing on AUC",  
		direction=">", quiet=TRUE, percent=TRUE, col=colors[[1]], grid=TRUE)

	rocobj2 <- plot.roc(data2$Response, data2$Predictor, direction=">", 
		add=TRUE, quiet=TRUE, percent=TRUE, col=colors[[2]], grid=FALSE)

	# Test on the whole AUC with Bootstrap method
	testobj2 <- roc.test(rocobj1, rocobj2, method="bootstrap", boot.n=boot.n)
	text(50, 50, labels = paste("      AUC\n     bootstrap\n p-value =", format.pval(testobj2$p.value)), adj=c(0, .5))
	legend("bottomright", legend=legendText, col=colors, lwd=2)
	dev.off()



	############# 3-Test on the pAUC_SP with Bootstrap method ###############
	pdf(paste0(output.dir,"ROC_test.bootstrap.pAUC_SP.pdf"), height=5,width=5)

	rocobj1 <- plot.roc(data1$Response, data1$Predictor, 
		main="Statistical hypothesis testing on pAUC_SP", 
		percent=TRUE, col=colors[[1]], quiet=TRUE,
		partial.auc=c(100,90),partial.auc.correct=TRUE,partial.auc.focus="sp",
		auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=FALSE)

	rocobj2 <- plot.roc(data2$Response, data2$Predictor, 
		add=TRUE, percent=TRUE, col=colors[[2]], quiet=TRUE,
		partial.auc=c(100,90),partial.auc.correct=TRUE,partial.auc.focus="sp",
		auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE, print.auc=FALSE)

	# Test on the partial AUC_SP with Bootstrap method
	testobj3 <- roc.test(rocobj1, rocobj2, method="bootstrap", boot.n=boot.n, partial.auc=c(1, 0.9), partial.auc.focus="sp")
	text(50, 50, labels = paste("      pAUC_SP\n     bootstrap\n p-value =", format.pval(testobj3$p.value)), adj=c(0, .5))
	legend("bottomright", legend=legendText, col=colors, lwd=2)
	dev.off()



	############# 4-Test on the pAUC_SE with Bootstrap method ###############
	pdf(paste0(output.dir,"ROC_test.bootstrap.pAUC_SE.pdf"), height=5,width=5)

	rocobj1 <- plot.roc(data1$Response, data1$Predictor, 
		main="Statistical hypothesis testing on pAUC_SE", quiet=TRUE,
		percent=TRUE, col=colors[[1]], partial.auc=c(100,90), partial.auc.correct=TRUE, partial.auc.focus="sens",
		auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=FALSE)


	rocobj2 <- plot.roc(data2$Response, data2$Predictor, 
		add=TRUE, percent=TRUE, col=colors[[2]], quiet=TRUE,
		partial.auc=c(100,90),partial.auc.correct=TRUE,partial.auc.focus="sens",
		auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE, print.auc=FALSE)

	# Test on the partial AUC_SE with Bootstrap method
	testobj4 <- roc.test(rocobj1, rocobj2, method="bootstrap", boot.n=boot.n, partial.auc=c(1, 0.9), partial.auc.focus="se")
	text(50, 50, labels = paste("      pAUC_SE\n     bootstrap\n p-value =", format.pval(testobj4$p.value)), adj=c(0, .5))
	legend("bottomright", legend=legendText, col=colors, lwd=2)
	dev.off()

	cat("\nINFO: all process completed.\n")
}