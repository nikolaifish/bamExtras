#' Geometric mean: product form
#'
#' Simple geometric mean function.This function will fail when x is a large value.
#' The logarithm form (geomean) works somewhat better but can't handle negative values.
#' @param x numeric vector
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' x <- c(1,2,10)
#' geomean(x)
#' geomean2(x)
#' mean(x)

# Note that this will fail when x is a large value. exp(mean(log(x))) works better but can't handle negative values
geomean2 <- function(x,na.rm=TRUE){
  if(na.rm){x <- as.numeric(na.omit(x))}
  prod(x)^(1/length(x))
}
