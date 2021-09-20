#' Von Bertalanffy (VB) age to length function
#'
#' Converts age to length with a Von Bertalanffy function
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param a ages to convert
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' vb_len()

vb_len <- function(Linf,K,t0,a)  {Linf*(1-exp(-K*(a-t0))) }
