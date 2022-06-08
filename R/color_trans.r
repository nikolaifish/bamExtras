#' Make opaque color transparent
#'
#' This function allows you to modify opaque colors you are already using and make transparent versions. It's handy for adding confidence bands around mean values.
#' @param x numeric vector
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' n <- 7
#' myalpha <- 0.25
#' mycol <- rainbow(n) # Opaque rainbow
#' mycol_tr <- rainbow(n,alpha=myalpha) # Transparent rainbow using alpha argument
#' mycol_tr2 <- color_trans(mycol,alpha=myalpha) # Transparent rainbow modifing mycol, using alpha argument from color_trans

#' par(cex=2,pch=16)
#' plot(1:n,rep(1,n),col=mycol,ylim=c(0,5))
#' points(1:n,rep(2,n),col=mycol_tr)
#' points(1:n,rep(3,n),col=mycol_tr2)

#' }

color_trans <- function(hexcolor,alpha=0.25){
  alpha <- max(min(alpha,1),0) # Constrain alpha between 0 and 1
  hex_alpha <- toupper(as.hexmode(round(alpha*255)))
  if(nchar(hex_alpha)==1){hex_alpha <- paste0("0",hex_alpha)}
  hexcolor[nchar(hexcolor)==7] <- paste0(hexcolor[nchar(hexcolor)==7],hex_alpha)
  hexcolor[nchar(hexcolor)==9] <- gsub(".{2}$",hex_alpha,hexcolor[nchar(hexcolor)==9])

  return(hexcolor)
}
