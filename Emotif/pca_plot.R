stat_ellipse <- function(mapping = NULL, data = NULL,
                         geom = "path", position = "identity",
                         ...,
                         type = "t",
                         level = 0.95,
                         segments = 51,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatEllipse,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      level = level,
      segments = segments,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatEllipse <- ggproto("StatEllipse", Stat,
  required_aes = c("x", "y"),

  compute_group = function(data, scales, type = "t", level = 0.95,
                           segments = 51, na.rm = FALSE) {
    calculate_ellipse(data = data, vars = c("x", "y"), type = type,
                      level = level, segments = segments)
  }
)

calculate_ellipse <- function(data, vars, type, level, segments){
  dfn <- 2
  dfd <- nrow(data) - 1

  if (!type %in% c("t", "norm", "euclid")) {
    message("Unrecognized ellipse type")
    ellipse <- rbind(as.numeric(c(NA, NA)))
  } else if (dfd < 3) {
    message("Too few points to calculate an ellipse")
    ellipse <- rbind(as.numeric(c(NA, NA)))
  } else {
    if (type == "t") {
      v <- MASS::cov.trob(data[,vars])
    } else if (type == "norm") {
      v <- stats::cov.wt(data[,vars])
    } else if (type == "euclid") {
      v <- stats::cov.wt(data[,vars])
      v$cov <- diag(rep(min(diag(v$cov)), 2))
    }
    shape <- v$cov
    center <- v$center
    chol_decomp <- chol(shape)
    if (type == "euclid") {
      radius <- level/max(chol_decomp)
    } else {
      radius <- sqrt(dfn * stats::qf(level, dfn, dfd))
    }
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    ellipse <- t(center + radius * t(unit.circle %*% chol_decomp))
  }

  ellipse <- as.data.frame(ellipse)
  colnames(ellipse) <- vars
  ellipse
}



rangescale <- function(X) {

	Xmax <- apply(X, 2, max)
	Xscaled = scale(X, scale=Xmax, center=T)

	return(Xscaled)
}
library('scales')
library("ggplot2")

pca_plot = function(input,output){
	df=read.table(input,sep=",",header=T)
	size_df = dim(df)	
		
	# Change the columns	
	data = data.matrix(df[,1:(size_df[2]-1)])	
	data = rangescale(data)	
	pca = prcomp(data)	
		
	sum_pca = summary(pca)	
		
	# Change the label if necessary	
	a=df[,size_df[2]]	
		
		
	pc_data_frame = data.frame(pc1 = pca$x[,1],pc2 = pca$x[,2],pc3 = pca$x[,3],class=a)	
	b=ggplot(pc_data_frame) +	
	  geom_point(aes(x=pc1, y=pc2, color=factor(class)), size=5, shape=20) +	
	  stat_ellipse(aes(x=pc1,y=pc2,fill=factor(class)),	
	               geom="polygon", level=0.95, alpha=0.2)+	
		xlab(paste("pc1 (",percent(sum_pca$importance[2,1]),")",sep = ""))+ ylab(paste("pc2 (",percent(sum_pca$importance[2,2]),")",sep = ""))	
	ggsave(output)	

}
file <- commandArgs(trailingOnly = TRUE);
pca_plot(file[1],file[2]);

