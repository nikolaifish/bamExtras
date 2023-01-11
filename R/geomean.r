#' Geometric mean: logarithm form
#'
#' Simple geometric mean function using only is.finite values of log(x).
#' @param x numeric vector
#' @param na.rm as in base R functions "a logical value indicating whether NA values should be stripped before the computation proceeds."
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' x <- c(1,2,10)
#' geomean(x)
#' geomean2(x)
#' mean(x)
#' }
#'

geomean <- function(x,na.rm=TRUE) {
  if(na.rm){x <- as.numeric(na.omit(x))}
  log.x <- log(x)
  out <- exp(mean(log.x[is.finite(log.x)]))
  return(out)
}
