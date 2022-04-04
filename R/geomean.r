#' Geometric mean: logarithm form
#'
#' Simple geometric mean function using only is.finite values of log(x).
#' @param x numeric vector
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

geomean <- function(x) {
  log.x <- log(x)
  out <- exp(mean(log.x[is.finite(log.x)]))
  return(out)
}
