#' Von Bertalanffy (VB) length to age function
#'
#' Converts age to length with a Von Bertalanffy function
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param L length to convert
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' vb_age()

vb_age <- function(Linf,K,t0,L) {(log(1-L/Linf)/-K)-t0}