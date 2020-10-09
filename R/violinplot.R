#' @title generating violinplots based on the combined PA results
#' @usage violinplot()
#' @description It displays 12 key metrics including the Youden's best p-value threshold, specificity
#'  (i.e., TNR), sensitivity (i.e., TPR), FDR, FPR, FNR, accuracy, precision, recall, AUC, pAUC_SP and 
#'  pAUC_SE across methods under study when benchmarked on a set of KO GEO datasets.
#' @return NULL
#' @export

violinplot <- function(){
	# Preprocessed .RData
	cat("\nINFO: summary files in current working directory:\n")
	rdatalist <- list.files(path=getwd(),pattern="*STATS.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter stats data index: ")
	dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
	if(length(ls(dataEnv))>1){
		stop("\nERROR: loaded more than one R object. Exit.\n")
	}
	DataObject <- get(ls(dataEnv),dataEnv)

	# Output into ./ViolinPlots/
	result_output <- "./ViolinPlots"
	if(!dir.exists(result_output)){
    	dir.create(result_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}
	result_output <- paste(result_output,"/",sep="")


	# Input data from DataObject
	aucTable <- DataObject$AUC
	colnames(aucTable) <- gsub("_AUC$","",colnames(aucTable))

	spTable <- DataObject$pauc_sp_cor
	colnames(spTable) <- gsub("_pAUC.*$","",colnames(spTable))

	seTable <- DataObject$pauc_se_cor
	colnames(seTable) <- gsub("_pAUC.*$","",colnames(seTable))

	statsTable <- DataObject$parameters
	colnames(statsTable) <- gsub("_RES","",colnames(statsTable))

	####################

	thresholdTable <- statsTable[,grep("threshold",colnames(statsTable))]
	colnames(thresholdTable) <- gsub(".threshold$","",colnames(aucTable))

	specificityTable <- statsTable[,grep("specificity",colnames(statsTable))]
	colnames(specificityTable) <- gsub(".specificity$","",colnames(specificityTable))

	sensitivityTable <- statsTable[,grep("sensitivity",colnames(statsTable))]
	colnames(sensitivityTable) <- gsub(".sensitivity$","",colnames(sensitivityTable))

	####################

	accuracyTable <- statsTable[,grep("accuracy",colnames(statsTable))]
	colnames(accuracyTable) <- gsub(".accuracy$","",colnames(accuracyTable))

	precisionTable <- statsTable[,grep("precision",colnames(statsTable))]
	colnames(precisionTable) <- gsub(".precision$","",colnames(precisionTable))

	recallTable <- statsTable[,grep("recall",colnames(statsTable))]
	colnames(recallTable) <- gsub(".recall$","",colnames(recallTable))

	####################

	fdrTable <- statsTable[,grep("fdr",colnames(statsTable))]
	colnames(fdrTable) <- gsub(".fdr$","",colnames(fdrTable))

	fprTable <- statsTable[,grep("fpr",colnames(statsTable))]
	colnames(fprTable) <- gsub(".fpr$","",colnames(fprTable))

	fnrTable <- statsTable[,grep("fnr",colnames(statsTable))]
	colnames(fnrTable) <- gsub(".fnr$","",colnames(fnrTable))

	#####################

	# Plotting
	plotAUCTable <- data.frame(stringsAsFactors=TRUE)
	plotAUCTable <- data.frame()
	for(i in 1:length(colnames(aucTable))){
		tmp <- data.frame(method=colnames(aucTable)[[i]],AUC=aucTable[[colnames(aucTable)[[i]]]],stringsAsFactors=TRUE)
		plotAUCTable <- rbind(plotAUCTable,tmp)
	}

	plotSPTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(spTable))){
		tmp <- data.frame(method=colnames(spTable)[[i]],pAUC_SP=spTable[[colnames(spTable)[[i]]]],stringsAsFactors=TRUE)
		plotSPTable <- rbind(plotSPTable,tmp)
	}

	plotSETable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(seTable))){
		tmp <- data.frame(method=colnames(seTable)[[i]],pAUC_SE=seTable[[colnames(seTable)[[i]]]],stringsAsFactors=TRUE)
		plotSETable <- rbind(plotSETable,tmp)
	}

	####################
	plotThresholdTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(thresholdTable))){
		tmp <- data.frame(method=colnames(thresholdTable)[[i]],Threshold=thresholdTable[[colnames(thresholdTable)[[i]]]],stringsAsFactors=TRUE)
		plotThresholdTable <- rbind(plotThresholdTable,tmp)
	}

	plotSpecificityTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(specificityTable))){
		tmp <- data.frame(method=colnames(specificityTable)[[i]],Specificity=specificityTable[[colnames(specificityTable)[[i]]]],stringsAsFactors=TRUE)
		plotSpecificityTable <- rbind(plotSpecificityTable,tmp)
	}

	plotSensitivityTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(sensitivityTable))){
		tmp <- data.frame(method=colnames(sensitivityTable)[[i]],Sensitivity=sensitivityTable[[colnames(sensitivityTable)[[i]]]],stringsAsFactors=TRUE)
		plotSensitivityTable <- rbind(plotSensitivityTable,tmp)
	}

	####################
	plotAccuracyTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(accuracyTable))){
		tmp <- data.frame(method=colnames(accuracyTable)[[i]],Accuracy=accuracyTable[[colnames(accuracyTable)[[i]]]],stringsAsFactors=TRUE)
		plotAccuracyTable <- rbind(plotAccuracyTable,tmp)
	}

	plotPrecisionTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(precisionTable))){
		tmp <- data.frame(method=colnames(precisionTable)[[i]],Precision=precisionTable[[colnames(precisionTable)[[i]]]],stringsAsFactors=TRUE)
		plotPrecisionTable <- rbind(plotPrecisionTable,tmp)
	}

	plotRecallTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(recallTable))){
		tmp <- data.frame(method=colnames(recallTable)[[i]],Recall=recallTable[[colnames(recallTable)[[i]]]],stringsAsFactors=TRUE)
		plotRecallTable <- rbind(plotRecallTable,tmp)
	}

	####################
	plotFdrTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(fdrTable))){
		tmp <- data.frame(method=colnames(fdrTable)[[i]],FDR=fdrTable[[colnames(fdrTable)[[i]]]],stringsAsFactors=TRUE)
		plotFdrTable <- rbind(plotFdrTable,tmp)
	}

	plotFprTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(fprTable))){
		tmp <- data.frame(method=colnames(fprTable)[[i]],FPR=fprTable[[colnames(fprTable)[[i]]]],stringsAsFactors=TRUE)
		plotFprTable <- rbind(plotFprTable,tmp)
	}


	plotFnrTable <- data.frame(stringsAsFactors=TRUE)
	for(i in 1:length(colnames(fnrTable))){
		tmp <- data.frame(method=colnames(fnrTable)[[i]],FNR=fnrTable[[colnames(fnrTable)[[i]]]],stringsAsFactors=TRUE)
		plotFnrTable <- rbind(plotFnrTable,tmp)
	}

	####################

	Mytheme<- theme_bw() + 
			theme(axis.text.y=element_text(size=6,color="black"),axis.title.y=element_text(face="bold",size=8,color="black")) + 
			theme(axis.text.x=element_text(angle=30,hjust=1,size=6,color="black"),axis.title.x=element_blank())

	p1_AUC <- ggplot(plotAUCTable, aes(x=method,y=AUC)) + 
			geom_violin() + 
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p2_pAUC_SP <- ggplot(plotSPTable, aes(x=method,y=pAUC_SP)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p3_pAUC_SE <- ggplot(plotSETable, aes(x=method,y=pAUC_SE)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme
	####################

	p4_Thre <- ggplot(plotThresholdTable, aes(x=method,y=Threshold)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p5_Spe <- ggplot(plotSpecificityTable, aes(x=method,y=Specificity)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p6_Sen <- ggplot(plotSensitivityTable, aes(x=method,y=Sensitivity)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme
	####################

	p7_Acc <- ggplot(plotAccuracyTable, aes(x=method,y=Accuracy)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p8_Pre <- ggplot(plotPrecisionTable, aes(x=method,y=Precision)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p9_Rec <- ggplot(plotRecallTable, aes(x=method,y=Recall)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme
	####################

	p10_FDR <- ggplot(plotFdrTable, aes(x=method,y=FDR)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p11_FPR <- ggplot(plotFprTable, aes(x=method,y=FPR)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	p12_FNR <- ggplot(plotFnrTable, aes(x=method,y=FNR)) + 
			geom_violin() +  
			geom_boxplot(width=.1,fill="black",outlier.color=NA) + 
			stat_summary(fun=median,geom="point",fill="white",shape=21,size=1.5) +
			Mytheme

	####################

	pdf(paste0(result_output,"plot.AUC.pdf"),width=length(colnames(aucTable))*1.3/2.54)
	plot(p1_AUC)
	dev.off()

	pdf(paste0(result_output,"plot.pAUC_SP.pdf"),width=length(colnames(spTable))*1.3/2.54)
	plot(p2_pAUC_SP)
	dev.off()

	pdf(paste0(result_output,"plot.pAUC_SE.pdf"),width=length(colnames(seTable))*1.3/2.54)
	plot(p3_pAUC_SE)
	dev.off()

	pdf(paste0(result_output,"plot.Threshold.pdf"),width=length(colnames(thresholdTable))*1.3/2.54)
	plot(p4_Thre)
	dev.off()

	pdf(paste0(result_output,"plot.Specificity.pdf"),width=length(colnames(specificityTable))*1.3/2.54)
	plot(p5_Spe)
	dev.off()

	pdf(paste0(result_output,"plot.Sensitivity.pdf"),width=length(colnames(sensitivityTable))*1.3/2.54)
	plot(p6_Sen)
	dev.off()

	pdf(paste0(result_output,"plot.Accuracy.pdf"),width=length(colnames(accuracyTable))*1.3/2.54)
	plot(p7_Acc)
	dev.off()

	pdf(paste0(result_output,"plot.Precision.pdf"),width=length(colnames(precisionTable))*1.3/2.54)
	plot(p8_Pre)
	dev.off()

	pdf(paste0(result_output,"plot.Recall.pdf"),width=length(colnames(recallTable))*1.3/2.54)
	plot(p9_Rec)
	dev.off()

	pdf(paste0(result_output,"plot.FDR.pdf"),width=length(colnames(fdrTable))*1.3/2.54)
	plot(p10_FDR)
	dev.off()

	pdf(paste0(result_output,"plot.FPR.pdf"),width=length(colnames(fprTable))*1.3/2.54)
	plot(p11_FPR)
	dev.off()

	pdf(paste0(result_output,"plot.FNR.pdf"),width=length(colnames(fnrTable))*1.3/2.54)
	plot(p12_FNR)
	dev.off()

	####################
	pdf(paste0(result_output,"plot.AUC.comb.pdf"), height=8/2.54,width=5/2.54)
	grid.arrange(p1_AUC,p2_pAUC_SP,p3_pAUC_SE,nrow=3,ncol=1)
	dev.off()

	pdf(paste0(result_output,"plot.STATS1.comb.pdf"), height=8/2.54,width=5/2.54)
	grid.arrange(p4_Thre,p5_Spe,p6_Sen,nrow=3,ncol=1)
	dev.off()

	pdf(paste0(result_output,"plot.STATS2.comb.pdf"), height=8/2.54,width=5/2.54)
	grid.arrange(p7_Acc,p8_Pre,p9_Rec,nrow=3,ncol=1)
	dev.off()

	pdf(paste0(result_output,"plot.STATS3.comb.pdf"), height=8/2.54,width=5/2.54)
	grid.arrange(p10_FDR,p11_FPR,p12_FNR,nrow=3,ncol=1)
	dev.off()

	####################
	pdf(paste0(result_output,"plot.STATS5.comb.pdf"), height=27/2.54,width=15/2.54)
	grid.arrange(p4_Thre,p10_FDR,
		         p5_Spe,p11_FPR,
		         p6_Sen,p12_FNR,
	    		 p7_Acc,p1_AUC,
	    		 p8_Pre,p2_pAUC_SP,
	    		 p9_Rec,p3_pAUC_SE,
	    		 nrow=6,ncol=2)
	dev.off()
}
