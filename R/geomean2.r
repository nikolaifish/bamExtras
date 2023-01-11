#' Geometric mean: product form
#'
#' Simple geometric mean function.
#' The logarithm form (geomean) works somewhat better but can only handle non-zero positive values.
#' @param x numeric vector
#' @param na.rm as in base R functions "a logical value indicating whether NA values should be stripped before the computation proceeds."
#' @param all.na.result Result to return if all values are NA. For some reason prod(numeric(0)) = 1, so the function will return 1 even if all x are NA.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' x <- c(1,2,10)
#' geomean(x)
#' geomean2(x)
#' mean(x)

geomean2 <- function(x,na.rm=TRUE,all.na.result=NA){
  if(all(is.na(x))){
    all.na.result
  }else{
  if(na.rm){x <- as.numeric(na.omit(x))}
  prod(x)^(1/length(x))
  }
}
