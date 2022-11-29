#' color function like rainbow() but without yellow (or indigo)
#'
#' color function like rainbow() but without yellow (or indigo)
#' @param n number of colors
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' n <- 7
#' myalpha <- 0.25
#' mycol <- ROGBP(n) # Opaque ROGBP colors
#' mycol_tr2 <- color_trans(mycol,alpha=myalpha) # Transparent ROGBP colors modifing mycol, using alpha argument from color_trans
#'
#' par(cex=2,pch=16)
#' plot(1:n,rep(1,n),col=mycol,ylim=c(0,5))
#' points(1:n,rep(3,n),col=mycol_tr2)
#' }
#'

ROGBP <- colorRampPalette(c("red","orange","green","blue","purple"))
