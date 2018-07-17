#' @title Textual preview for FCS objects
#'
#' @description Prints a preview for a FCS object.
#'
#' @param x a FCS object
#' 
#' @return None
#' 
#' @name print
#' @rdname print-methods
NULL

#' @rdname print-methods
#' @export
setMethod("print","FCS",
    function(x){
        cat("Object class: FCS\n")
        cat(paste0("Object name: ",x@name,"\n"))
        cat(paste0("Number of cells: ",format(x@cell.nb,big.mark=","),"\n"))
        cat(paste0("Number of markers: ",length(x@markers),"\n"))
        cat("Markers: ")
        cat(x@markers,sep="; ")
        cat("\n")
    }
)


#' @title Textual preview for FCS objects
#'
#' @description Shows a preview for a FCS object.
#'
#' @param object a FCS object
#' 
#' @return None
#' 
#' @name show
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show","FCS", function(object) print(object))

