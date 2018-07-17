#' @title Plot for FCS objects
#'
#' @description Makes visual representations for FCS objects.
#'
#' @param object1 a FCS to plot
#' @param object2 a FCS to plot
#' @param ... other parameters
#'
#' @return returns a ggplot object
#' 
#' @import ggplot2
#' 
#' @name plot
#' @rdname plot-methods
#' @export 
setGeneric("plot", function(object1,object2,...){ standardGeneric("plot") })

#' @rdname plot-methods
#' @export
setMethod("plot",c("FCS","missing"),
    function(object1,object2){
        				
		intensities <- apply(object1@intensities,2,mean)
		markers     <- object1@markers
			
		frame <- data.frame(profile          = 1,
							markers          = markers,
							intensities      = intensities, 
							stringsAsFactors = FALSE)
		
		frame$markers <- factor(frame$markers,levels = frame$markers)                
		
		ymin <- min(frame$intensities)
		ymax <- max(frame$intensities)
		
		plots <- ggplot2::ggplot(frame,ggplot2::aes_string(x="markers",y="intensities",group="profile",ymin="intensities",ymax="intensities")) +
			ggplot2::ggtitle(paste0(object1@name))+
			ggplot2::geom_point(stat="identity",colour="chartreuse3") + 
			ggplot2::geom_line(stat="identity",colour="chartreuse3") +
			ggplot2::xlab("markers") +
			ggplot2::ylab("expression") +
			ggplot2::scale_y_continuous(limits=c(floor(ymin),ceiling(ymax)),breaks=seq(floor(ymin),ceiling(ymax),by=1)) +
			ggplot2::theme(axis.text.x=ggplot2::element_text(angle=290, hjust=0, vjust=1),
				legend.position="bottom",
				panel.background=ggplot2::element_blank(),
				panel.grid.major.y=ggplot2::element_line(colour="gray75"),
				panel.grid.major.x=ggplot2::element_line(colour="gray75"),
				panel.grid.minor.x=ggplot2::element_line(colour="gray90"),
				plot.title=ggplot2::element_text(hjust = 0.5))

		return(plots)

    }
)


#' @rdname plot-methods
#' @export
setMethod("plot",c("FCS","FCS"),
    function(object1,object2){
	
		markers     <- intersect(object1@markers,object2@markers)
		
		frame.profile1 <- data.frame(profiles    = 1,
								 	 markers     = markers,
									 intensities = apply(object1@intensities[,object1@markers %in% markers],2,mean),
									 stringsAsFactors = FALSE)
									 
		frame.profile2 <- data.frame(profiles    = 2,
									 markers     = markers,
									 intensities = apply(object2@intensities[,object2@markers %in% markers],2,mean),
									 stringsAsFactors = FALSE)
									 
		plots <- ggplot2::ggplot() +
			ggplot2::ggtitle(paste0(object1@name," vs. ",object2@name))+
			ggplot2::geom_point(data=frame.profile1,ggplot2::aes_string(x="markers",y="intensities"),colour="deepskyblue3") + 
			ggplot2::geom_point(data=frame.profile2,ggplot2::aes_string(x="markers",y="intensities"),colour="firebrick3") + 
			ggplot2::geom_line(data=frame.profile1,ggplot2::aes_string(x="markers",y="intensities",group="profiles"),colour="deepskyblue3") +
			ggplot2::geom_line(data=frame.profile2,ggplot2::aes_string(x="markers",y="intensities",group="profiles"),colour="firebrick3") +
			ggplot2::xlab("markers") +
			ggplot2::ylab("expression") +
			ggplot2::scale_color_manual(values=c("deepskyblue3","firebrick3")) +
            ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE)) +
				ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1),
				legend.position="bottom",
				panel.background=ggplot2::element_blank(),
				panel.grid.major.y=ggplot2::element_line(colour="gray75"),
				panel.grid.major.x=ggplot2::element_line(colour="gray75"),
				panel.grid.minor.x=ggplot2::element_line(colour="gray90"),
				plot.title=ggplot2::element_text(hjust = 0.5))
			
		return(plots)
    }
)

