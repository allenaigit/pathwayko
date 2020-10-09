#' @title performing pairwise comparisons on the PA results
#' @usage wilcoxtest()
#' @description It conducts the "paired", pairwise comparisons between each pair of two methods under study 
#'  when benchmarked on the same set of KO datasets. pairwise.wilcox.test is from the R `stats` package, 
#'  with corrections (BH FDR adjustment of p-value) for multiple testing under the two-sided mode.
#' @return NULL
#' @export

wilcoxtest <- function(){
	cat("\nINFO: summary files in current working directory:\n")
	rdatalist <- list.files(path=getwd(),pattern="*STATS.RData",full.names=TRUE,recursive=TRUE,include.dirs=FALSE)
	show(rdatalist)
	input <- readline(prompt = "\nEnter stats data index: ")
	dataEnv <- loadToEnvironment(rdatalist[as.numeric(input)])
	if(length(ls(dataEnv))>1){
		stop("\nERROR: loaded more than one R object. Exit.\n")
	}
	statsObject <- get(ls(dataEnv),dataEnv)


# reading data tables for individual test 
	AUC_table <- statsObject$AUC
	pAUC_SP_ori_table <- statsObject$pauc_sp_ori
	pAUC_SE_ori_table <- statsObject$pauc_se_ori
	pAUC_SP_cor_table <- statsObject$pauc_sp_cor
	pAUC_SE_cor_table <- statsObject$pauc_se_cor

	data1 <- rearrangeTable(AUC_table,"AUC")
	data2 <- rearrangeTable(pAUC_SP_ori_table,"pAUC_SP_ori")
	data3 <- rearrangeTable(pAUC_SE_ori_table,"pAUC_SE_ori")
	data4 <- rearrangeTable(pAUC_SP_cor_table,"pAUC_SP_cor")
	data5 <- rearrangeTable(pAUC_SE_cor_table,"pAUC_SE_cor")

# creating a directory for Wilcoxtest and sink a file

	result_output <- "./Wilcoxtests"
	if(!dir.exists(result_output)){
    	dir.create(result_output, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}
	result_output <- paste(result_output,"/",sep="")
	file <- paste0(result_output,"wilcox.out.txt")

# running each section of individual test - Wilcoxon, Shapiro and t 

	cat("\nINFO: starting...\n")

	sink(file=file,append=FALSE)	# clear for new file wilcox.out.txt
	cat("\n====== INFO_1: Paired-pairwise Wilcoxon Testing for AUC...======\n\n")
	pairwise.wilcox.AUC <- pairwise.wilcox.test(data1$AUC,data1$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.wilcox.test for AUC
		)
	show(pairwise.wilcox.AUC)
	sink()


	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_2: Paired-pairwise Wilcoxon Testing for pAUC_SP_ori...======\n\n")
	pairwise.wilcox.pAUC_SP_ori <- pairwise.wilcox.test(data2$pAUC_SP_ori,data2$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.wilcox.test for pAUC_SP_ori
		)
	show(pairwise.wilcox.pAUC_SP_ori)
	sink()


	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_3: Paired-pairwise Wilcoxon Testing for pAUC_SE_ori...======\n\n")
#	attach(data3)
	pairwise.wilcox.pAUC_SE_ori <- pairwise.wilcox.test(data3$pAUC_SE_ori,data3$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.wilcox.test for pAUC_SE_ori
		)
	show(pairwise.wilcox.pAUC_SE_ori)
	sink()


	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_4: Paired-pairwise Wilcoxon Testing for pAUC_SP_cor...======\n\n")
	pairwise.wilcox.pAUC_SP_cor <- pairwise.wilcox.test(data4$pAUC_SP_cor,data4$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.wilcox.test for pAUC_SP_cor
		)
	show(pairwise.wilcox.pAUC_SP_cor)
	sink()


	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_5: Paired-pairwise Wilcoxon Testing for pAUC_SE_cor...======\n\n")
#	attach(data5)
	pairwise.wilcox.pAUC_SE_cor <- pairwise.wilcox.test(data5$pAUC_SE_cor,data5$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.wilcox.test for pAUC_SE_cor
		)
	show(pairwise.wilcox.pAUC_SE_cor)
	sink()


# Optional t test if Shapiro normality test approved 

	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt	
	cat("\n====== INFO_6: Paired-pairwise t Testing for AUC...======\n\n")
	shapiro_AUC <- tapply(data1$AUC,data1$METHOD,shapiro.test) # Shapiro normality test for AUC
	show(shapiro_AUC)
	# Conclusion: At the significance level (p-value > 0.1), the normality is not true. 

	pairwise.t.AUC <- pairwise.t.test(data1$AUC,data1$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.t.test for AUC
		)
	show(pairwise.t.AUC)
	sink()


	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_7: Paired-pairwise t Testing for pAUC_SP_ori...======\n\n")
	shapiro_pAUC_SP_ori <- tapply(data2$pAUC_SP_ori,data2$METHOD,shapiro.test) # Shapiro normality test for pAUC_SP
	show(shapiro_pAUC_SP_ori)
	# Conclusion: At the significance level (p-value > 0.1), the normality is not true.

	pairwise.t.pAUC_SP_ori <- pairwise.t.test(data2$pAUC_SP_ori,data2$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.t.test for pAUC_SP_ori
		)
	show(pairwise.t.pAUC_SP_ori)
	sink()


	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_8: Paired-pairwise t Testing for pAUC_SE_ori...======\n\n")
#	attach(data3)
	shapiro_pAUC_SE_ori <- tapply(data3$pAUC_SE_ori,data3$METHOD,shapiro.test) # Shapiro normality test for pAUC_SE
	show(shapiro_pAUC_SE_ori)
	# Conclusion: At the significance level (p-value > 0.1), the normality is not true.

	pairwise.t.pAUC_SE_ori <- pairwise.t.test(data3$pAUC_SE_ori,data3$METHOD,
		p.adjust.method = "BH",
		pool.sd = FALSE, 
		paired = TRUE,
		alternative = ("two.sided") # pairwise.t.test for pAUC_SE_ori
		)
	show(pairwise.t.pAUC_SE_ori)
	sink()



	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_9: Paired-pairwise t Testing for pAUC_SP_cor...======\n\n")
	shapiro_pAUC_SP_cor <- tapply(data4$pAUC_SP_cor,data4$METHOD,shapiro.test) # Shapiro normality test for pAUC_SP
	show(shapiro_pAUC_SP_cor)
	# Conclusion: At the significance level (p-value > 0.1), the normality is not true.

	tryCatch(
		{
			pairwise.t.pAUC_SP_cor <- pairwise.t.test(data4$pAUC_SP_cor,data4$METHOD,
				p.adjust.method = "BH",
				pool.sd = FALSE, 
				paired = TRUE,
				alternative = ("two.sided") # pairwise.t.test for pAUC_SP_cor
				)
			show(pairwise.t.pAUC_SP_cor)
		},
			error = function(cond) {
	            cat("Error when processing: pairwise.t.pAUC_SP_cor")
	            cat("Here's the original error message:") # find an exceptional case & handle it ?????
	            print(cond)
	            cat("\n")
	            return(NULL)
        	}
        )	
	sink()



	sink(file=file,append=TRUE) # continue on the file wilcox.out.txt
	cat("\n\n====== INFO_10: Paired-pairwise t Testing for pAUC_SE_cor...======\n\n")
#	attach(data5)
	shapiro_pAUC_SE_cor <- tapply(data5$pAUC_SE_cor,data5$METHOD,shapiro.test) # Shapiro normality test for pAUC_SE
	show(shapiro_pAUC_SE_cor)
	# Conclusion: At the significance level (p-value > 0.1), the normality is not true.


	tryCatch(
		{
			pairwise.t.pAUC_SE_cor <- pairwise.t.test(data5$pAUC_SE_cor,data5$METHOD,
				p.adjust.method = "BH",
				pool.sd = FALSE, 
				paired = TRUE,
				alternative = ("two.sided") # pairwise.t.test for pAUC_SE_cor
				)
			show(pairwise.t.pAUC_SE_cor)
		},
			error = function(cond) {
	            cat("Error when processing: pairwise.t.pAUC_SE_cor")
	            cat("Here's the original error message:") # find an exceptional case & handle it ?????
	            print(cond)
	            cat("\n")
        	}
        )
	sink()

	cat("\nINFO: completed.\n")
}



rearrangeTable <- function(table,header){
	colnames(table) <- gsub(paste0("_",header,"$"),"",colnames(table))
	X <- vector()
	METHOD <- vector()
	NAMES <- vector()
	for(i in 1:ncol(table)){
		t_x <- table[,i]
		t_names <- rownames(table)
		t_method <- rep(x=colnames(table)[i],times=nrow(table))
		X <- c(X,t_x)
		METHOD <- c(METHOD,t_method)
		NAMES <- c(NAMES,t_names)
	}
	data <- data.frame(cbind(NAMES,X,METHOD))
	data$X <- as.numeric(data$X)
	colnames(data) <- c("NAMES",header,"METHOD")
	return(data)
}