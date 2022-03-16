#' Geometric mean: logarithm form
#'
#' Simple geometric mean function using only is.finite values of log(x).
#' @param x numeric vector
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' n <- 5
#' mycol <- rainbow(n)
#' mycol_tr <- color_trans(mycol,alpha=0.5)
#' mycol2 <- rgb(1,seq(0,1,length=n),0)
#' mycol2_tr <- color_trans(mycol2)
#'
#' par(cex=2,pch=16)
#' plot(1:n,rep(1,n),col=mycol,ylim=c(0,5))
#' points(1:n,rep(2,n),col=mycol_tr)
#' points(1:n,rep(3,n),col=mycol2)
#' points(1:n,rep(4,n),col=mycol2_tr)
#' }

color_trans <- function(hexcolor,alpha=0.25){
  alpha <- max(min(alpha,1),0) # Constrain alpha between 0 and 1
  hex_alpha <- toupper(as.hexmode(round(alpha*255)))
  hexcolor[nchar(hexcolor)==7] <- paste0(hexcolor[nchar(hexcolor)==7],hex_alpha)
  hexcolor[nchar(hexcolor)==9] <- gsub(".{2}$",hex_alpha,hexcolor[nchar(hexcolor)==9])

  return(hexcolor)
}
