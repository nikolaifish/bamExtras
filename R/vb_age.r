#' Von Bertalanffy (VB) length to age function
#'
#' Converts age to length with a Von Bertalanffy function
#' @param L lengths to convert to ages
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' vb_age()

vb_age <- function(L,Linf,K,t0) {(log(1-L/Linf)/-K)-t0}
