# title Internal - Inverse hyperbolic sine (arcsinh) transformation
#
# @description Computes the inverse hyperbolic sine (arcsinh) of a real number: log(sqrt((x/coeff)^2+x/coeff)), where 'coeff' is a rescaling coefficient (cofactor). The default value of 'coeff' is set to 5.
#
# @details This transformation is an alternative to the log transformation for cytometry data processing.
#
# @param x a numeric vector of values to transform
# @param coeff a numeric value indicating the rescaling coefficient (cofactor)
#
# @return a numeric vector of transformed values
arcsinh <- function(x,coeff=5){
    res <- log((x/coeff)+sqrt((x/coeff)^2+1))
    return(res)
}

# title Internal - hyperbolic sine (sinh) transformation
#
# @description Computes the hyperbolic sine (sinh) of a real number: (exp(x)-exp(-x))/2*coeff, where 'coeff' is a rescaling coefficient (cofactor). The default value of 'coeff' is set to 5.
#
# @details This transformation is used to convert values from an inverse hyperbolic sine (arcsinh) transformation
#
# @param x a numeric vector of values to transform
# @param coeff a numeric value indicating the rescaling coefficient (cofactor)
#
# @return a numeric vector of transformed values
sinh_coeff <- function(x,coeff=5){
    res <- (exp(x)-exp(-x))/2*coeff
    return(res)
}


# title Internal - Creation of a FlowFrame object
#
# @description Creates a FlowFrame object based on a numeric matrix of FCS profiles. 
# 
# @details This function is used internally when writing FCS profiles contained in a FCS object into a FCS file.
#
# @param intensities a numeric matrix corresponding to the FCS profile intensities.
# @param markers a character vector corresponding to marker names.
#
# @return None
createFlowFrame <- function(intensities,markers) {
    
    colnames(intensities) <- markers
    p                     <- c()    
    description           <- list() 
    
    description[["$DATATYPE"]] <- "F"
    
    for (i in 1:ncol(intensities)) {
        name  <- markers[i]
        min   <- min(intensities[,i])
        max   <- max(intensities[,i])
        range <- max-min+1
        
        l           <- matrix(c(name,name,range,min,max),nrow=1)
        colnames(l) <- c("name","desc","range","minRange","maxRange")
        rownames(l) <- paste0("$P",i) 
        p           <- rbind(p,l)
        
        description[[paste("$P",i,"N",sep="")]] <- name;
        description[[paste("$P",i,"S",sep="")]] <- name;
        description[[paste("$P",i,"R",sep="")]] <- toString(range);
        description[[paste("$P",i,"B",sep="")]] <- "32";
        description[[paste("$P",i,"E",sep="")]] <- "0,0";
    }
    
    dataframe <- as(data.frame(p), "AnnotatedDataFrame")
    flowframe <- flowCore::flowFrame(intensities, dataframe, description=description)
    
    return(flowframe)
}


# title Internal - Extract a dictionary from a FCS file
#
# @description This function is used internally to extract the correspondence between the original marker names (first column) and the true marker names (second column)
#
# @param fcs.file a character indicating the location of the fcs file containing the correspondences
#
# @return a two-column data.frame providing the correspondence between the original marker names (first column) and the true marker names (second column)
extract.dictionary <- function(fcs.file){
    flowframe       <- flowCore::read.FCS(fcs.file[1],trans = FALSE)
    dictionary      <- flowframe@parameters@data[, c(1, 2)]
    dictionary[, 1] <- make.names(dictionary[, 1])
    return(dictionary)
}


# title Internal - Renaming cell markers
#
# @description This function is used internally to rename the cell markers based on a dictionary.
#
# @details dictionary is a data.frame used to rename the marker names. The first column must correspond to the original marker names, the second column must correspond to the new marker names. 
#
# @param header a character vector containing the original maker names
# @param dictionary a character vector containing a correspondence between the original and the new marker names
#
# @return a character vector containing the renamed marker names
rename.markers <- function(header,dictionary){
    header         <- make.names(header)
    dictionary[,1] <- as.vector(dictionary[,1])
    dictionary[,2] <- as.vector(dictionary[,2])
    if(length(unique(dictionary[,1]))!=length(dictionary[,1])){
        stop("Duplicate in dictionary 'original marker names'")
    }
    occurences <- table(dictionary[,2])
    occurences <- occurences[occurences>1]
    redondant.names <- names(occurences)
    for (i in 1:nrow(dictionary)) {
        if (any(dictionary[i, 2] %in% redondant.names)) {
            temp             <- gsub("X\\.(.*)\\.Di(_clust$|$)","(\\1)Di\\2",dictionary[i,1])
            dictionary[i,2] <- paste0(dictionary[i,2],"-",temp)
        }
		if(!is.na(dictionary[i,2])){
			header[which(header == dictionary[i,1])[1]] <- dictionary[i,2]
		}
    }
    return(header)
}


# title Internal - Removing of cell markers to exclude from a matrix
#
# @description This function is used internally to remove one or several cell markers from a numeric matrix.
#
# @param data a numeric matrix
# @param exclude a character vector containing the cell markers to be excluded
#
# @return a numeric matrix without the cell markers to exclude
exclude.markers <- function(data,exclude){
    exclude.flags <- exclude %in% colnames(data)
    if(any(!(exclude.flags))){
        warning(paste0("Unknown marker to exclude: ",paste(exclude[!exclude.flags],collapse=", ")))
    }
    data    <- data[,!(colnames(data) %in% exclude)]
    return(data)
}


merge_header <- function(FCS){
	bb              <- FCS@intensities[,grep("BB_",FCS@markers)]
	colnames(bb)    <- FCS@markers[grep("BB_",FCS@markers)]
	res <- c()
	M   <- c()
	for(i in seq(1,ncol(bb),by=2)){
		tmp <- apply(cbind(bb[,i],bb[,i+1]),1,mean)
		res <- cbind(res,tmp)
		M   <- c(M,colnames(bb)[i])
	}
	colnames(res) <- gsub(".FCS1","",M)
	evo             <- FCS@intensities[,grep("BB_",FCS@markers,invert=TRUE)]
	colnames(evo)   <- FCS@markers[grep("BB_",FCS@markers,invert=TRUE)]
	colnames(evo)   <- gsub(".FCS1","",colnames(evo))
	colnames(evo)   <- gsub(".FCS2","",colnames(evo))
	FCS  <- cbind(res,evo)
	FCS <- as.FCS(FCS)
	return(FCS)
}


