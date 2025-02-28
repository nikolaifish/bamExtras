#' Phase plots for bootstrapped data
#'
#' Bivariate scatter plots, specially formatted for making phase plots for stock and fishery status results.
#' @param x x values from bootstrapping. numeric vector
#' @param y y values from bootstrapping. numeric vector
#' @param xref x reference value (e.g. base run value).
#' @param yref y reference value (e.g. base run value).
#' @param col_point color for (x,y) points
#' @param col_quant color for quantile lines
#' @param col_ref color for (xref,yref) points
#' @param col_text color for percentages plotted in each quadrant
#' @param probs quantile probabilities to compute length of error bars.
#' @param text_x_adj text (percentages in each quadrant) adjustment to x
#' @param text_y_adj text (percentages in each quadrant) adjustment to y
#' @param ... other parameters to pass to \code{plot}
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' nsim <- 1000
#' x <- rnorm(n=nsim,mean=1.1,sd=.2)
#' y <- 1*(1/x)+rnorm(n=nsim,mean=0.1,sd=.1)
#' par(mfrow=c(2,2))
#' plot_boot_phase(x,y)
#' plot_boot_phase(x,y,xref=1.2,yref=1.1)
#' plot_boot_phase(x,y,xref=1.2,yref=1.1,xlab="F/Fmsy",ylab="SSB/SSBmsy")
#' plot_boot_phase(x,y,xref=1.2,yref=1.1,xlab="F/Fmsy",ylab="SSB/SSBmsy",col_point=rainbow(nsim,alpha=0.3))


plot_boot_phase <- function(x,
                            y,
                            xref=NULL,
                            yref=NULL,
                            col_point=rgb(0,0,0,0.2),
                            col_quant="green4",
                            col_ref="black",
                            col_text="blue",
                            probs=c(0.05,0.95),
                            text_x_adj=c(0,0,0,0),
                            text_y_adj=c(0,0,0,0),
                            ...){
  quad1=round(length(x[x>1 & y>1])/length(x),3)*100
  quad2=round(length(x[x<1 & y>1])/length(x),3)*100
  quad3=round(length(x[x<1 & y<1])/length(x),3)*100
  quad4=round(length(x[x>1 & y<1])/length(x),3)*100

  plot(x, y,
       xlim=range(c(0,1,x)), ylim=range(c(0,1,y)),
       pch=16,col=col_point,...)
  abline(h=1,lty=2)
  abline(v=1, lty=2)
  xmed <- quantile(x,0.5,na.rm = TRUE)
  ymed <- quantile(y,0.5,na.rm = TRUE)
  points(rep(xmed,2), quantile(y,probs=probs,na.rm = TRUE), type="l", lty=1, col=col_quant, lwd=2)
  points(quantile(x,probs=probs,na.rm = TRUE),rep(ymed,2), type="l", lty=1, col=col_quant, lwd=2)
  if(!is.null(xref)&!is.null(yref))
  points(xref,yref,pch=16,cex=1.5,col=col_ref)

  usr <- par("usr")

  text(x=c(rep((usr[1]+1)/2,2),rep((usr[2]+1)/2,2))+text_x_adj,
       y=rep(c((usr[3]+1)/2,(usr[4]+1)/2),2)+text_y_adj,
       labels=paste(c(quad3,quad2,quad4,quad1),"%",sep=""),col="blue")
}
