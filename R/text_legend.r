#' Wrapper for legend function to quickly label a plot with text
#'
#' @param legend see \code{\link[base]{legend}}
#' @param x see \code{\link[base]{legend}}
#' @param bty see \code{\link[base]{legend}}
#' @param adj see \code{\link[base]{legend}}
#' @param cex see \code{\link[base]{legend}}
#' @param text.col see \code{\link[base]{legend}}
#' @param ... pass additional arguments to see \code{legend}
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' x <- 1000
#' plot(rnorm(x),rnorm(x),col=rainbow(x),pch=16)
#' text_legend("random rainbow dots")
#' plot(rnorm(x),rnorm(x),col=heat.colors(x),pch=16)
#' text_legend("BAM!!",x="center",text.col=1,cex=4)

text_legend <- function(legend,x="topleft",bty="n",adj=c(0.1,0),cex=1.3,text.col=rgb(0,0,0.3,0.4),...){
  legend(x=x,
         legend=legend,
         adj=adj
         ,bty=bty,
         cex=cex,
         text.col=text.col,
         ...
  )
}
