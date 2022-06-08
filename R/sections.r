#' Identify sections of a vector of consecutive values whose differences does not exceed a buffer
#'
#' Identify sections of a vector that have the same consecutive value
#' @param x numeric vector
#' @param buffer value compared with difference between consecutive values. If the absolute difference exceeds the buffer, the values are considered different and not part of the same section.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' x <- sample(c(1,0),100,replace=TRUE)
#' sections(x=x)
#' }
#'

sections <- function(x,buffer=0) {
  # Identify vector of changes from previous value. First value is always zero (no change)
  changes <- c(0,(abs(diff(x))>buffer)*1)
  # Assign section numbers
  sct <- 1+cumsum(changes)
  return(sct)
}
