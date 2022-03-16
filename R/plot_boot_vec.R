#' Plot stochastic vectors
#'
#' Plot median values and confidence intervals (shaded region) for stochastic (e.g. bootstrapped) vectors.
#' @param data a numeric matrix or data frame that can be coerced to a numeric matrix where rows are runs and columns represent a  numeric sequence (e.g. age, year). Note that the function would probably work for any numeric matrix, but it's intended for bootstrapped time- or age-series data.
#' @param CIpct percentage to define the range of the confidence intervals to plot
#' @param runsToKeep optional vector of integer values corresponding to the rows of `data` to include in the plot
#' @param ref_x,ref_y x and y vectors for adding a reference line to the plot (e.g. base run values)
#' @param col_lines color for lines depicting central tendency and confidence limits. Expects hex color
#' @param col_shade color for shading confidence bands
#' @param add logical. Should a new plot be drawn or should it be added to an existing plot?
#' @param ... other arguments to pass to `matplot`
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' bootData <- t(replicate(n=100,random_walk(mu=0.5)))
#' plot_boot_vec(bootData)
#' }

plot_boot_vec <- function(data,
                          CIpct = 95,
                          runsToKeep=NULL,
                          ref_x=NULL,
                          ref_y=NULL,
                          col_lines=rgb(0,0,0,1),
                          col_shade=NULL,
                          add=FALSE,
                          ...){
  if(is.null(runsToKeep)){
    runsToKeep <- 1:nrow(data)
  }
  if(is.null(col_shade)){
    col_lines
    col_shade <- color_trans(col_lines)
  }
  CItail <- ((100-CIpct)/2)/100
  #par(mar=c(2,2,2,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1.5,tck=-0.02)
  data2 <- t(data[runsToKeep,])
  y.median <- apply(X=data2,MARGIN=1,FUN=quantile,probs=0.5,na.rm=TRUE)
  y.lo <- apply(X=data2,MARGIN=1,FUN=quantile,probs=CItail,na.rm=TRUE)
  y.up <- apply(X=data2,MARGIN=1,FUN=quantile,probs=1-CItail,na.rm=TRUE)
  y <- data.frame(y.median,y.lo,y.up)
  x <- as.numeric(rownames(y))
  # draw empty plot
  if(!add){
  matplot(x,y,type="n",...)
  grid(lty=1,col="gray80")
  }

  ## draw confidence bands as shaded region
  # Identify groups of continuous years for plotting polygons correctly
  group <- rep(0,nrow(y))
  for(i in 1:nrow(y))
  {
    if(complete.cases(y[i,])) # If test i is true
    {
      if(i==1) # if i=1 then this is group 1
      {
        group[i] <- 1
      }
      if(i>1){ # if i>1
        if(complete.cases(y[i-1,])) # If test i-1 is true
        {
          group[i] <- group[i-1]
        }else{
          group[i] <- max(group)+1
        }
      }
    }
  }
  y$group <- group

  ycc <- y[complete.cases(y),]
  for(group.i in 1:max(ycc$group)){
    ycc.group.i <-  ycc[ycc$group==group.i,]
    polygon(x=c(rownames(ycc.group.i),rev(rownames(ycc.group.i))),y=c(ycc.group.i[,2],rev(ycc.group.i[,3])),col=col_shade,border=col_shade,lwd=2)
  }

  # draw median and confidence band lines
  matpoints(rownames(y),y[,c(1,2,3)],type="l",lwd=c(2,1,1),lty=c(2,1,1),col=col_lines)
  # draw reference lines (e.g. base run)
  if(!is.null(ref_x)&!is.null(ref_y)){
    points(ref_x,ref_y,type="o",lty=1,lwd=3,pch=16,col=col_lines)
  }
  legend("topright",legend=paste("error bands represent ",CIpct,"% CI",sep=""),bty="n",cex=0.75)
}
