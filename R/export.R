#' @title Exportation of FCS objects
#'
#' @description Exports a FCS object into a FCS file.
#'
#' @param object a FCS object
#' @param filename a character indicating the location of the FCS file to save
#' @param transform a logical indicating if a sinh transformation must be performed when exporting the FCS
#' @param ... aadditional parameters
#'
#' @return None
#' 
#' @importFrom flowCore write.FCS
#' @import flowUtils
#'
#' @name export
#' @rdname export-methods
setGeneric("export",function(object,filename,...) { standardGeneric("export") })

#' @rdname export-methods
#' @export
setMethod("export",c("FCS"),
    function(object,filename="FCS.fcs",transform=TRUE){
		if(transform==TRUE){
			data                  <- sinh_coeff(object@intensities)
		}
		flowFrame                 <- createFlowFrame(data,object@markers)
        suppressWarnings(filename <- flowCore::write.FCS(flowFrame,filename))
    }
)
