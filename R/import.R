#' @title Importation of cytometry profiles from one or several FCS files
#' 
#' @description Imports one or several cytometry profiles from a FCS file or from a set of files into a FCS object.
#'
#' @details If a set of files is specified, then the files are merged in the import procedure.
#'
#' Several transformations can be applied on  marker expressions via the 'trans' parameter:
#' 
#' The transformation functions can be parametrized using the named list 'trans.para'. The scale (cofactor) of the arcsinh transformation function can be parametrized using the 'arcsinh.scale' value. The shift of the log transformation function can be parametrized using the 'log.shift' value and the base of the log transformation function can be parametrized using the 'log.base' value. The 'log.shift' allows the value "auto" which automatically identify the log shift avoiding to apply log transformations on negative values.
#' 
#' @importFrom flowCore exprs read.FCS
#'
#' @param path a character vector indicating the location to a FCS file or to a set of FCS files
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#' @param trans a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation
#' @param trans.para a named list containing parameters for the transformation. Please refer to the details section for more details
#' @param trans.exclude a character vector containing the marker names for which no transformation must be applied on (including the rescaling transformation)
#'
#' @return a S4 object of class FCS
#'
#' @export
import.FCS <- function(path,
                        exclude           = NULL,
                        trans             = "arcsinh",
						trans.para	      = switch(trans,
												   "arcsinh" = list(arcsinh.scale=5),
												   "log"     = list(log.shift="auto",log.base=10),
												   "none"    = NULL),
                        trans.exclude     = "cluster"){
    
    if(length(path)==0){
        stop("Error in import.FCS: path has a length of 0")
    }else if(length(path)==1 && is.na(file.info(path)$isdir)){
        stop(paste0("Error in import.FCS: ",path," do not exist"))
    }else if(length(path)>1 && !any(file.exists(path))){
        stop(paste0("Error in import.FCS: ",file[!any(file.exists(path))]," do not exist"))
    }
    
    if(length(path)==1 && file.info(path)$isdir==TRUE){
        path <- list.files(path,full.names = TRUE,pattern=".+\\.fcs$")
    }
    
    markers    <- c()
    FCS        <- list()
    dictionary <- extract.dictionary(path)
	    
	for(i in 1:length(path)){
        message(paste0("Importing ",path[i]))
        suppressWarnings(data <- flowCore::exprs(flowCore::read.FCS(path[i],trans = FALSE)))
        
		colnames(data)        <- rename.markers(colnames(data),dictionary)
		
		if(!is.null(exclude)){
            data <- exclude.markers(data,exclude)
        }
        markers <- colnames(data)[!(colnames(data) %in% trans.exclude)]
        if(trans=="arcsinh"){
            for(marker in markers){
                data[,marker] <- arcsinh(data[,marker],trans.para$arcsinh.scale)
            }    
        }else if(!is.element(trans,c("log","none"))){
            stop("Unknown transformation to apply")
        }
        
        markers  <- colnames(data)
        
        if(nrow(data)>1)
            profiles <- as.character(1:nrow(data))
        else    
            profiles <- as.character(NA)
        
        rownames(data) <- profiles

        FCS[[i]]     <- as.FCS(data)
        
    }
    FCS  <- do.call("c",FCS)
    
    if(trans=="log"){
        to.transform <- match(setdiff(markers,trans.exclude),markers)
        if(trans.para$log.shift=="auto"){
            min                  <- min(FCS@intensities[,to.transform])
            trans.para$log.shift <- ifelse(min>0,0,abs(floor(min-1)))
        }
        FCS@intensities[,to.transform] <- log(FCS@intensities[,to.transform] + trans.para$log.shift,trans.para$log.base)
    }
	
	FCS@name          <- paste(gsub(".fcs","",basename(path)),sep="",collapse=";")
    FCS@trans         <- trans
    FCS@trans.para    <- trans.para
    FCS@trans.exclude <- trans.exclude
    return(FCS)
}


#' @title Importation of cytometry profiles from a SPADE analysis
#' 
#' @description Imports one cytometry profiles from a SPADE analysis into a FCS object.
#'
#' @details If a set of files is specified, then the files are merged in the import procedure.
#' 
#' @importFrom flowCore exprs read.FCS
#'
#' @param path a character vector indicating the location to a FCS file or to a set of FCS files
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#'
#' @return a S4 object of class FCS
#'
#' @export
import.SPADE <- function(path, exclude = c("density","cluster")){

	fcs.files <- list.files(path,full.names = TRUE,pattern=".+\\.cluster\\.fcs$")
	res <- import.FCS(fcs.files[1])
	
	data <- data.frame(res@intensities)
	data <- sinh_coeff(data,5)
	
	colnames(data)[ncol(data)] <- "cluster"
	
	data <- stats::aggregate(.~cluster, data=data, mean)
	
	data <- data[,c(2:ncol(data),1)]
	
	data <- arcsinh(data,5)
	
	res@intensities <- as.matrix(data)
	res@cell.nb     <- nrow(res@intensities) 
	
	res             <- res[,grep(paste0(exclude,collapse="|"),res@markers,invert=TRUE)]
	
	return(res)

}
