#' @title FCS class definition
#'
#' @description FCS are S4 objects containing marker cell expressions obtained from cytometry profiles.
#'
#' @details This object mainly stores for each cell the intensities of all cell markers.
#'
#' The slot 'trans.para' is a named list contains different parameters depending of the transformation applied on the marker expression intensities. The scale (cofactor) of the arcsinh transformation function is parametrized using the 'arcsinh.scale' value. The shift of the log transformation function is parametrized using the 'log.shift' value and the base of the log transformation function is parametrized using the 'log.base' value. If no transformation function have been applied, the 'trans.para' slot is set to NULL.
#'
#' @slot name a character indicating the internal name of the FCS object
#' @slot profiles a character vector containing the names of the FCS profiles 
#' @slot cell.nb an integer value indicating the number of FCS profiles
#' @slot markers a character vector containing the marker names
#' @slot markers.nb an integer value indicating the number of markers
#' @slot intensities a numeric matrix containing the intensities of each marker for each FCS profile
#' @slot trans a character specifying the name of a transformation function applied on the marker expression intensities. Possible values are "arcsinh" for arc sin hyperbolic transformation, "log" for logarithmic transformation, or "none" for no transformation
#' @slot trans.para a named list containing parameters of the transformation. Please refer to the details section for more details
#' @slot trans.exclude a character vector containing the marker names for which no transformation has been applied on
#'
#' @import methods
#'
#' @name FCS-class
#' @rdname FCS-class
#' @exportClass FCS
FCS <- setClass("FCS",
    slots=c(name          = "character",
		profiles          = "character",
		cell.nb           = "integer",
		markers           = "character",
		markers.nb        = "integer",
		intensities       = "matrix",
		trans             = "character",
		trans.para        = "ANY",
		trans.exclude     = "ANY"),
    validity=function(object){
        if(length(object@profiles)!=length(unique(object@profiles)))
            stop("Error in profiles slot: profile names are not unique")
        if(length(object@markers)!=length(unique(object@markers)))
            stop("Error in markers slot: marker names are not unique")
        if(object@cell.nb!=length(object@profiles))
            stop("Error in cell.nb slot: cell.nb do not correspond to the number of profile names")
        if(object@markers.nb!=length(object@markers))
            stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
        if(!is.element(object@trans,c("arcsinh","log","none")))
            stop("Error in trans slot: trans do not contain allowed value (allowed value are \"arcsinh\", \"log\", \"none\")")
        if(!is.null(object@trans.para) && class(object@trans.para)[1]!="list")
            stop("Error in trans.para slot: trans.para must be a of type list or NULL")
        if(!is.null(object@trans.exclude) && class(object@trans.exclude)[1]!="character")
            stop("Error in trans.exclude slot: trans.exclude must be a of type character or NULL")
        return(TRUE)
    }
)
setMethod("initialize",c("FCS"),
    function(.Object,
            name              = "",
            profiles          = "",
            cell.nb           = 0,
            markers           = "",
            markers.nb        = 0,
            intensities       = as.matrix(0)){
        if(cell.nb==0)
            stop("Error can not create a FCS object with no profile")
        .Object@name              = name     
        .Object@profiles          = profiles
        .Object@cell.nb           = cell.nb
        .Object@markers           = markers
        .Object@markers.nb        = markers.nb
        .Object@trans             = "none"
		.Object@trans.para        = NULL
		.Object@trans.exclude     = NULL
		.Object@intensities       = intensities
        methods::validObject(.Object)
        return(.Object)
    }
)



#' @title Extraction of subsets of data from FCS objects
#'
#' @description Extracts a subset of data from a FCS object.
#'
#' @details For The parameter i represents a vector of cells to extract and the parameter j represents a vector of markers to extract.
#'
#' @param x a FCS object
#' @param i a numeric, logical or character vector
#' @param j a numeric, logical or character vector
#'
#' @return a S4 object of class FCS
#'
#' @name extract
#' @rdname extract-methods
NULL

#' @rdname extract-methods
#' @export
setMethod("[",c("FCS","ANY","ANY"),
    function(x,i,j){

		if(!missing(i) && length(i)>x@cell.nb)
            stop("Too many cluster profiles to extract")
        if(!missing(j) && length(j)>length(x@markers))
            stop("Too many markers to extract")
        if(!missing(i) && !(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(!missing(j) && !(is.logical(j) || is.numeric(j) || is.character(j)))
            stop(paste0("Wrong types in j: ",typeof(j)))
                
        if(missing(i)){
            k <- 1:x@cell.nb
        }else{
            if(is.character(i)){
                k <- which(x@profiles %in% i)
            }else{
                k <- i
            }
        }
        
        if(missing(j)){
            l <- 1:length(x@markers)
        }else{
            if(is.character(j)){
                l <- unlist(sapply(paste0("^",j,"$"),grep,x=x@markers))
            }else{
                l <- j
            }    
        }
        
        FCS <- FCS(name        = x@name,
            profiles           = x@profiles[k],
            cell.nb            = length(x@profiles[k]),
            markers            = x@markers[l],
            markers.nb         = length(x@markers[l]),
            intensities        = x@intensities[k,l,drop=FALSE])
            
        return(FCS)
    }
)


# Combination of FCS objects
#
# @description Combines two or several FCS objects.
#
# @param x a first FCS object
# @param ... further FCS objects to be combined
#
# @return a S4 object of class FCS
#
# @name c
# @rdname c-methods
NULL

# @rdname c-methods
# @export
setMethod("c",c("FCS"),
    function(x,...){
    
        other.FCS  <- list(x,...)
        name        <- c()
        markers     <- x@markers
        cell.nb <- 0
        
        i <- 1
        for(FCS in other.FCS){
            if(class(FCS) != "FCS")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(FCS)))
            name        <- c(name,FCS@name)
            markers     <- union(markers,FCS@markers)
            cell.nb <- cell.nb+FCS@cell.nb
            i <- i+1
        }
        name        <- paste0(name,collapse=";")
        profiles    <- as.character(1:cell.nb)
        
        intensities <- matrix(NA,ncol=length(markers),nrow=cell.nb,dimnames=list(profiles,markers))
        
        i <- 1
        for(FCS in other.FCS){
            nb                                   <- FCS@cell.nb
            intensities[i:(i+nb-1),FCS@markers] <- FCS@intensities
            i <- i+nb
        }
        dimnames(intensities) <- NULL
        
        FCS <- FCS(name = name,
            profiles      = profiles,
            cell.nb   = length(profiles),
            markers       = markers,
            markers.nb    = length(markers),
            intensities   = intensities)
        return(FCS)
    }
)


# @title Coercion to a FCS object
#
# @description Coerces a numeric matrix into a FCS object.
#
# This function transforms a numeric matrix into one FCS object.
#
# @details The matrix must have its column names corresponding to the cell markers.
#
# @param object a numeric matrix
# @param name a character specifying the internal name of the FCS object to create
#
# @return a S4 object of class FCS
#
# @name as.FCS
# @rdname as.FCS-methods
#
# @export
setGeneric("as.FCS", function(object,name="FCS") { standardGeneric("as.FCS") })

# @rdname as.FCS-methods
# @export
setMethod("as.FCS",c("matrix"),
    function(object,name){  
        data           <- object
        dimnames(data) <- NULL
        FCS <- FCS(name = name,
            profiles      = as.character(1:nrow(object)),
            cell.nb       = nrow(object),
            markers       = colnames(object),
            markers.nb    = ncol(object),
            intensities   = data)
        return(FCS)
    }
)

