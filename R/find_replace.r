#' Find and replace single or multiple values in a vector. Usage is similar to gsub except that this function will not replace parts of character strings
#'
#' @param pattern vector of values to find
#' @param replacement vector of values to replace
#' @param x a character vector where matches are sought
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' find_replace()

find_replace <- function(pattern,replacement,x)  {
  found   <- match(x, pattern)

  ifelse(is.na(found), x, replacement[found])
  }
