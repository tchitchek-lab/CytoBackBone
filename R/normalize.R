#' title Normalization of CELL objects
#'
#' @description Performs a normalization between two of more CELL objects.
#'
#' @param list a list of CELL objects to normalize
#' @param BBmarkers a character vector indicating the marker names on which the normalization will be applied
#'
#' @return a list of normalized CELL objects
#'
#' @export
normalizeQuantile <- function(list,BBmarkers){
	
	i   <- 1
	for(cell in list){
		if(class(cell) != "FCS")
			stop(paste0("Element at position ",i," is not of type FCS"))
		i   <- i+1
	}
	
	min <- Inf
	max <- -Inf
	for(cell in list){
		min <- min(min,cell@cell.nb)
		max <- max(max,cell@cell.nb)
	}
	
	i <- 1
	for(cell in list){
		if(max - cell@cell.nb==0){
			list[[i]] <- cell
		}else{
			tofillNA <- matrix(rep(NA, (max - cell@cell.nb)*length(cell@markers)),ncol=length(cell@markers),nrow=max - cell@cell.nb)
			colnames(tofillNA) <- cell@markers
			tofillNA  <- as.FCS(tofillNA)
			list[[i]] <- c(cell,tofillNA)
		}
		
		i <- i+1
	}
	
	
	for(markers in BBmarkers){
		cat(paste0("processing marker: ",markers,"\n"))
		input <- c()
		for(cell in list){
			input <- cbind(input,cell@intensities[,cell@markers==markers])
		}
		res <- preprocessCore::normalize.quantiles(as.matrix(input))
		i <- 1
		for(cell in list){
			idx = list[[i]]@markers==markers
			list[[i]]@intensities[,idx] <- res[,i]
			i <- i+1
		}
	}
	
	i <- 1
	for(cell in list){
		data <- as.matrix(cell@intensities)
		data <- stats::na.omit(data)
		colnames(data) <- cell@markers
		list[[i]]      <- as.FCS(data)
		i = i+1
	}

	invisible(list)

}

