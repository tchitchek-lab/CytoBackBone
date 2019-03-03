#' @title Merging two FCS objects
#' 
#' @description Perfomrs a merging of two FCS objects.
#' CytoBackBone uses a k Nearest Neighbor distances algorithm, based on k-d tree space partitioning, to identify cells that are acceptable and non-ambiguous nearest neighbors.
#' 
#' To be acceptable neighbors, the phenotypic distance between two cells must be lower than a specific threshold, defined by the user, which corresponds to a Euclidian distance computed as the square root of the sum of the squared expression differences for each pair of common markers.
#' 
#' To be defined as non-ambiguous nearest neighbors, two cells from two different profiles must be reciprocally the closest neighbors. First, the algorithm identifies the closest cells in cytometric profile #2 for each cell from cytometric profile #1. Then, the algorithm identifies the closest cells in cytometric profile #1 for each cell from cytometric profile #2. Merging is possible if two cells from different cytometric profiles are identified as mutual non-ambiguous nearest neighbors.
#'
#' @details More precisely, the algorithm works as follows: (i) cells of each cytometric profile with no acceptable neighbors are first excluded from the two input profiles; (ii) all acceptable and non-ambiguous nearest neighbor cells remaining in the two profiles are merged into a new profile and are discarded from the two input profiles; (iii) the previous step of this procedure is repeated until no more acceptable and non-ambiguous neighbor cells can be found; and (iv) finally, all excluded and remaining cells are stored in supplementary cytometric profiles for information purposes. These successive iterations ensure that the algorithm finds a new set of acceptable and non-ambiguous neighbor cells at each step. The number of iteration performed by the algorithm varies depending on the dataset used, but usually ranges from 2 to 30 interactions.
#'
#' A distance threshold defining acceptable nearest neighbors needs to be specified to avoid merging data between two cells that have very different levels for the backbone markers: the lower the threshold, the more accurate the merging. However, the threshold must not be too stringent as, by definition, cells cannot have exactly the same phenotype. On the contrary, merging can result in a significant loss of cells. Thus, the threshold must be carefully set to allow CytoBackBone to merge data from cells with similar levels of the backbone markers within a reasonable error range. The Euclidean distance between cells increases with the number of backbone markers. Therefore, the threshold value must be adjusted according to the length of this backbone.
#'
#' Cells discarded (having no acceptable and not-ambiguous neighbours) during the merging procedure (i.e., cells specific to the two cytometric profiles to merge) can be extracted using the `leftout` argument. 
#'
#' @param FCS1 a first FCS object to merge
#' @param FCS2 a second FCS object to merge
#' @param BBmarkers a character indicating the markers to use as backbone markers
#' @param th a numerical value indicating the threshold distance to use to consider two cells as acceptable neighbors
#' @param normalize a logical indicating if the expressions of backbone markers must be quantile-normalized
#' @param leftout a logical indicating if left-out cells (cells not merged during the merging procedure) must also be extracted in the output results. If set to TRUE then the output will be a named list containing the merged cytometric profile (named "merged"), the cells specific to the first profile (named "specific.FCS1"), the cells specific to the second profile (named "specific.FCS2"), and the cells specific to both profiles (named "specific").
#'
#' @return a S4 object of class FCS
#'
#' @export
merge <- function(FCS1, FCS2, BBmarkers=NULL, th = length(BBmarkers)*0.20, normalize=TRUE, leftout=FALSE) {

	name1 = FCS1@name
	name2 = FCS2@name
	
	FCS1 <- FCS1[,sort(FCS1@markers)]
	FCS2 <- FCS2[,sort(FCS2@markers)]
	
	cat("=== CytoBackBone ===\n")

	if(normalize==TRUE){
		cat(paste0("Normalizing profiles\n"))
		normalizedProfiles <- normalizeQuantile(list(FCS1,FCS2),BBmarkers)
		FCS1               <- normalizedProfiles[[1]]
		FCS2               <- normalizedProfiles[[2]]
	}
	
	FCS1@markers[FCS1@markers %in% BBmarkers]        <- paste0("BB_",FCS1@markers[FCS1@markers %in% BBmarkers])
	FCS2@markers[FCS2@markers %in% BBmarkers]        <- paste0("BB_",FCS2@markers[FCS2@markers %in% BBmarkers])
	FCS1markersnotbb                                 <- FCS1@markers[!(FCS1@markers %in% paste0("BB_",BBmarkers))]
	FCS2markersnotbb                                 <- FCS2@markers[!(FCS2@markers %in% paste0("BB_",BBmarkers))]
	
	if(length(intersect(FCS1markersnotbb,FCS2markersnotbb))>0){
		message("----------ERROR----------")
		message("backbone markers are not consistent between the two cytometric profiles")
		message("backbone markers: ",appendLF=FALSE)
		message(paste0(BBmarkers,collapse=","))
		message("specific markers for FCS1: ",appendLF=FALSE)
		message(paste0(FCS1markersnotbb,collapse=","))
		message("specific markers for FCS2: ",appendLF=FALSE)
		message(paste0(FCS2markersnotbb,collapse=","))
		message("intersection: ",appendLF=FALSE)
		intersect = intersect(FCS1markersnotbb,FCS2markersnotbb)
		message(paste0(intersect,collapse=","))
		message("----------ERROR----------")
		stop()
	}
		
	FCS1markersnotbb                                 <- FCS1markersnotbb[FCS1markersnotbb %in% FCS2markersnotbb]
	FCS2markersnotbb                                 <- FCS2markersnotbb[FCS2markersnotbb %in% FCS1markersnotbb]
	FCS1@markers[FCS1@markers %in% FCS1markersnotbb] <- paste0(FCS1@markers[FCS1@markers %in% FCS1markersnotbb],".1")
	FCS2@markers[FCS2@markers %in% FCS2markersnotbb] <- paste0(FCS2@markers[FCS2@markers %in% FCS2markersnotbb],".2")
		
	FCS1                        <- FCS1[, FCS1@markers[grep("cluster",FCS1@markers,invert=TRUE)]]
	FCS2                        <- FCS2[, FCS2@markers[grep("cluster",FCS2@markers,invert=TRUE)]]
	cat(paste0("profile 1 contains ",format(FCS1@cell.nb,big.mark=",")," cells\n"))
	cat(paste0("profile 2 contains ",format(FCS2@cell.nb,big.mark=",")," cells\n"))
	cat("===\n")
	FCS1_BB                     <- FCS1[,FCS1@markers[grep("BB_",FCS1@markers)]]@intensities
	colnames(FCS1_BB)           <- FCS1@markers[grep("BB_",FCS1@markers)]
	FCS2_BB                     <- FCS2[,FCS2@markers[grep("BB_",FCS2@markers)]]@intensities
	colnames(FCS2_BB)           <- FCS2@markers[grep("BB_",FCS2@markers)]
	dist                        <- FNN::knnx.dist(FCS2_BB,FCS1_BB,k=1,algorithm="kd_tree")
	FCS1_excluded               <- NULL
	if(sum(dist>th)!=0){
		FCS1_excluded           <- FCS1[dist>th]
	}
	FCS1                        <- FCS1[dist<th]
	dist                        <- FNN::knnx.dist(FCS1_BB,FCS2_BB,k=1,algorithm="kd_tree")
	FCS2_excluded               <- NULL
	if(sum(dist>th)!=0){
		FCS2_excluded           <- FCS2[dist>th]
	}
	FCS2                        <- FCS2[dist<th]
	cat(paste0("profile 1 has ",format(FCS1@cell.nb,big.mark=",")," cells can be potentialy matched\n"))
	cat(paste0("profile 2 has ",format(FCS2@cell.nb,big.mark=",")," cells can be potentialy matched\n"))
	cat("===\n")
	
	max <- min(FCS1@cell.nb,FCS2@cell.nb)
	cat(paste0("maximum of number of cells that can be matched by CytoBackBone = ",format(max,big.mark=",")),"\n")
	cat("===\n")
		
	FCS1_BB                       <- FCS1[,FCS1@markers[grep("BB_",FCS1@markers)]]@intensities
	colnames(FCS1_BB)             <- FCS1@markers[grep("BB_",FCS1@markers)]
	FCS2_BB                       <- FCS2[,FCS2@markers[grep("BB_",FCS2@markers)]]@intensities
	colnames(FCS2_BB)             <- FCS2@markers[grep("BB_",FCS2@markers)]
	knnx                          <- FNN::get.knnx(FCS2_BB,FCS1_BB,k=1,algorithm="kd_tree")
	idx                           <- knnx$nn.index
	dist                          <- knnx$nn.dist
	table_FCS1_to_FCS2            <- cbind(1:nrow(idx),idx)
	table_FCS1_to_FCS2[dist>th,2] <- -1
	colnames(table_FCS1_to_FCS2)  <- c("FCS1","FCS2")
	knnx                          <- FNN::get.knnx(FCS1_BB,FCS2_BB,k=1,algorithm="kd_tree")
	idx                           <- knnx$nn.index
	dist                          <- knnx$nn.dist
	table_FCS2_to_FCS1            <- cbind(1:nrow(idx),idx)
	table_FCS2_to_FCS1[dist>th,2] <- -1
	colnames(table_FCS2_to_FCS1)  <- c("FCS2","FCS1")
	idx                           <- paste0(table_FCS1_to_FCS2[,1],"-",table_FCS1_to_FCS2[,2]) %in% paste0(table_FCS2_to_FCS1[,2],"-",table_FCS2_to_FCS1[,1])
	idx_a                         <- table_FCS1_to_FCS2[idx,"FCS1"]
	idx_b                         <- table_FCS1_to_FCS2[idx,"FCS2"]
	data                          <- cbind(FCS1@intensities[idx_a,],FCS2@intensities[idx_b,])
	colnames(data)                <- c(paste0(FCS1@markers,".FCS1"),paste0(FCS2@markers,".FCS2"))
	cat(paste0("step #",1,": ",format(nrow(data),big.mark=",")," cells matched (",format(max-nrow(data),big.mark=",")," cells remaining)\n"))
	
	if(max-nrow(data)!=0){	
		i <- 2
		while(TRUE){
			
			FCS1                          <- FCS1[-idx_a]
			FCS2                          <- FCS2[-idx_b]
			FCS1_BB                       <- FCS1[,FCS1@markers[grep("BB_",FCS1@markers)]]@intensities
			colnames(FCS1_BB)             <- FCS1@markers[grep("BB_",FCS1@markers)]
			FCS2_BB                       <- FCS2[,FCS2@markers[grep("BB_",FCS2@markers)]]@intensities
			colnames(FCS2_BB)             <- FCS2@markers[grep("BB_",FCS2@markers)]
			knnx                          <- FNN::get.knnx(FCS2_BB,FCS1_BB,k=1,algorithm="kd_tree")
			idx                           <- knnx$nn.index
			dist                          <- knnx$nn.dist
			table_FCS1_to_FCS2            <- cbind(1:nrow(idx),idx)
			table_FCS1_to_FCS2[dist>th,2] <- -1
			colnames(table_FCS1_to_FCS2)  <- c("FCS1","FCS2")
			knnx                          <- FNN::get.knnx(FCS1_BB,FCS2_BB,k=1,algorithm="kd_tree")
			idx                           <- knnx$nn.index
			dist                          <- knnx$nn.dist
			table_FCS2_to_FCS1            <- cbind(1:nrow(idx),idx)
			table_FCS2_to_FCS1[dist>th,2] <- -1
			colnames(table_FCS2_to_FCS1)  <- c("FCS2","FCS1")
			idx                           <- paste0(table_FCS1_to_FCS2[,1],"-",table_FCS1_to_FCS2[,2]) %in% paste0(table_FCS2_to_FCS1[,2],"-",table_FCS2_to_FCS1[,1])
			idx_a                         <- table_FCS1_to_FCS2[idx,"FCS1"]
			idx_b                         <- table_FCS1_to_FCS2[idx,"FCS2"]
			
			if(sum(idx)==0)
				break
			if(sum(idx)==1)
				data_s                    <- cbind(t(FCS1@intensities[idx_a,]),t(FCS2@intensities[idx_b,]))
			else
				data_s                    <- cbind(FCS1@intensities[idx_a,],FCS2@intensities[idx_b,])
			data                          <- rbind(data,data_s)
			
			cat(paste0("step #",i,": ",format(nrow(data),big.mark=",")," cells matched (",format(max-nrow(data),big.mark=",")," cells remaining)","\n"))
			i <- i+1
			
			if((max-nrow(data))==0)
				break
				
		}
	}
	
	if(!is.null(FCS1_excluded)){
		FCS1_excluded          <- c(FCS1_excluded,FCS1)
		FCS1_excluded@markers  <- gsub("BB_","",FCS1_excluded@markers)
		FCS1_excluded@name     <- paste0("cells specific to ",name1)
	}
	
	if(!is.null(FCS2_excluded)){
		FCS2_excluded          <- c(FCS2_excluded,FCS2)
		FCS2_excluded@markers  <- gsub("BB_","",FCS2_excluded@markers)
		FCS2_excluded@name     <- paste0("cells specific to ",name2)
	}
	
	excluded = NULL
	if(!is.null(FCS1_excluded) && !is.null(FCS2_excluded)){
		excluded               <- c(FCS1_excluded,FCS2_excluded)
		excluded@markers       <- gsub("BB_","",excluded@markers)
		excluded@name          <- paste0("cells specific to ",name1," or ",name2)
	}
	
	FCS                    <- as.FCS(data)
	FCS                    <- FCS[,sort(FCS@markers)]
	FCS                    <- merge_header(FCS)
	FCS@markers            <- gsub("BB_","",FCS@markers)
	FCS@name               <- paste0(name1," + ",name2)
	
	cat("====================\n")
	
	if(leftout==FALSE){
		return(FCS)
	}else{
		return(list(merged=FCS,specific.FCS1=FCS1_excluded,specific.FCS2=FCS2_excluded,specific=excluded))
	}
	
}
