#' Stacked barplot developed for plotting landings data by fleet and year
#'
#' @param data data frame object. Intended for rows (x-axis) to be years, columns (groups) to be fleets, and values (y-axis) to be landings, discards, etc.
#' @param units Text for printing the scale of values in ylab (e.g. fish, lb). character vector
#' @param unit_scale Allows you to provide a numeric value to scale the values for plotting
#' @param ylab y-label for plot. Note that if unit_scale!=0 then the value of unit_scale will be appended to the y-axis label.
#' @param plot_prop_text should text be added to the plot indicating what proportion of the total values (sum(data)) are in each group? logical
#' @param prop_digits Number of digits to display if plot_prop_text==TRUE
#' @param xlabBy Specify frequency for plotting values along x-axis (e.g. 1, 5, 10)
#' @param leg_title Title text for labeling legend. character vector
#' @param leg_bty Pass to argument \code{bty} in \code{legend} (see \code{\link[base]{legend}})
#' ... 
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' barplot_stacked(rdat_BlackSeaBass$N.age,units="fish",unit_scale = 1e6,ylab="N",leg_title = "Age")
#' barplot_stacked(rdat_BlackSeaBass$N.age/1e6,units="million fish",ylab="N",leg_title = "Age")
#' barplot_stacked(rdat_BlackSeaBass$N.age/1e6,units="million fish",ylab="N",plot_prop_text=FALSE,leg_bty = "o")
#' }

# Stacked bar plot
barplot_stacked <- function(data,
                            units,
                            unit_scale=1,
                            ylab="y-axis label",
                            plot_prop_text=TRUE,
                            prop_digits=2,
                            xlabBy=5,
                            leg_title="",
                            leg_bty="n",
                            ...){
  par(mar=c(4,4,1,1),mgp=c(2,.5,0),oma=c(0,0,0,0),tck=-0.01,cex.lab=2,cex.axis=1,xpd=FALSE)
  barcol <- rainbow(ncol(data),start=0,end=0.8)
  data[is.na(data)] <- 0
  data <- data/unit_scale
  ylab <- paste0(ylab," (",ifelse(unit_scale==1,"",paste(unit_scale,"")),units,")")
  mids <- barplot(t(as.matrix(data)),col=barcol,ylab=ylab,space=0,axisnames=FALSE,...)
  xaxislab <- local({
    a <- unique(round(as.numeric(rownames(data))/xlabBy)*xlabBy)
    a[a%in%rownames(data)]
  })
  xaxislab.mids <- mids[match(xaxislab,rownames(data))]
  axis(side=1,at=xaxislab.mids,labels=paste(xaxislab))
  abline(h=pretty(par("usr")[3:4]),v = xaxislab.mids,lty=1,col=rgb(0,0,0,.1))
  
  propByGroup <- signif(colSums(data,na.rm=TRUE)/sum(data,na.rm=TRUE),prop_digits)
  
  if(plot_prop_text){
    legend("topright",bty=leg_bty,legend=paste(names(propByGroup)," (",propByGroup,")",sep=""),cex=1,col=barcol,pch=15,pt.cex=1.5,title=leg_title)
    legend("right",bty=leg_bty,legend="proportions of total\n are in parentheses",cex=0.75,col="blue")
  }else{
    legend("topleft",bty=leg_bty,legend=paste(names(propByGroup),sep=""),cex=1,col=barcol,pch=15,pt.cex=1.5,title=leg_title)
  }
}
