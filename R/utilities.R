# https://www.r-bloggers.com/safe-loading-of-rdata-files-2/
loadToEnvironment <- function(RData=NULL, env = new.env()){
	if(is.null(RData)){
		cat("Warning: RData cannot be null in loadToEnvironment.")
    return(NULL)
	}
  	load(RData, env)
  	return(env) 
}



# parse multiple KEGG xml within a directory
parseKEGGxml <- function(path_to_xml = NULL, mc.cores=1){
	if(is.null(path_to_xml)){
		return(NULL)
  	}

  	pwy.files <- list.files(path_to_xml,pattern="*.xml",full.names=TRUE)
  	pwys <- NULL

  	pwys <- mclapply(pwy.files, KEGGgraph::parseKGML, mc.cores = mc.cores, mc.preschedule = TRUE, mc.silent = TRUE)
  	
  	# lapply scrambles rownames, has to be reset to correct names
  	names(pwys) <- pwy.files  	
	names(pwys) <- gsub(path_to_xml,"",names(pwys))
	names(pwys) <- gsub(".xml$","",names(pwys))
	names(pwys) <- gsub("^/","",names(pwys))
	return(pwys)
}



# Convert parsed KEGG xmls into graph objects
keggPathwayGraphs <- function(pwys = NULL, mc.cores=1)
{
  	if(is.null(pwys)){
      return(NULL)
	}
  	

    l <- lapply(pwys, function(path) {
        l <- sapply(edges(path), getType)
        if(length(l) == 0)
            return(0)
        t <- table(l)
        return(t)
    })
  
    allRelTypes <- unique(unlist(lapply(l, names)))
  
    counts <- do.call(rbind,lapply(l, function(x) as.vector(x[allRelTypes])))
    colnames(counts) <- allRelTypes
  
    counts[is.na(counts)] <- 0
  
    accIndex <- rowSums(counts[,c("GErel","PCrel","PPrel")], na.rm=T) / rowSums(counts, na.rm=T) >= 0.9
    accIndex[is.na(accIndex)] <- FALSE
    pwys <- pwys[accIndex]
  
    names(pwys) <- sapply(pwys, getName)
  
    pathwayGraphs <- NULL

	pathwayGraphs <- mclapply(pwys, function(g) 
	{
        g <- KEGGpathway2Graph(g, expandGenes=TRUE)

        if (class(g) != "graphNEL") {
            return(NULL)
        }
    
        kg <- new("graphNEL", nodes(g), edges(g), edgemode = "directed")
    
        if (length(getKEGGedgeData(g)) == 0)
        {
              return(NULL)
        }
        edgeDataDefaults(kg, "subtype") <- NA
    
    
        cf <- cbind(
            do.call(rbind, strsplit(names(getKEGGedgeData(g)), '~')),
            sapply(getKEGGedgeData(g), function(e) 
            paste(lapply(getSubtype(e), getName), collapse=","))
        );
    
        if (nrow(cf) < 2) {
            ucf <- cf
        } else {
            ucf <- cf[unique(rownames(cf)),]
            ucf[,3] <- tapply(cf[,3], rownames(cf), function(ll) return(paste( unique(ll), collapse = ',')))[rownames(ucf)]
        }
    
    
        relGeneTable <- data.frame(ucf, stringsAsFactors = FALSE)
        names(relGeneTable) <- c("from","to","subtype")
    
        edgeData(kg, relGeneTable$from, relGeneTable$to, "subtype") <- relGeneTable$subtype
    
        return(kg)
    },
    mc.cores = mc.cores,
    mc.preschedule = TRUE,
    mc.silent = TRUE)

	pathwayGraphs <- pathwayGraphs[!sapply(pathwayGraphs, is.null)]
	return(pathwayGraphs)
}



# https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
getFreeMemoryMB <- function() {
  osName <- Sys.info()[["sysname"]]
  if (osName == "Windows") {
    x <- system2("wmic", args =  "OS get FreePhysicalMemory /Value", stdout = TRUE)
    x <- x[grepl("FreePhysicalMemory", x)]
    x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
    x <- gsub("\r", "", x, fixed = TRUE)
    return(as.integer(x)/1024)
  } else if (osName == 'Linux') {
    x <- system2('free', args='-k', stdout=TRUE)
    x <- strsplit(x[2], " +")[[1]][4]
    return(as.integer(x)/1024)
  } else {
    stop("Only supported on Windows and Linux")
  }
}