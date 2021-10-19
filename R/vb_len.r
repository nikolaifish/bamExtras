#' Von Bertalanffy (VB) age to length function
#'
#' Converts age to length with a Von Bertalanffy function
#' @param a ages to convert to lengths
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param aP proportion of age (value between 0 and 1) at which to compute length (e.g. aP=0.5 to compute length at midyear)
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Length at age for Black Sea Bass
#' vb_len(a=0:11,Linf = 502, K = 0.173, t0 = -0.97)
#' }

vb_len <- function(a,Linf,K,t0,aP=0)  {Linf*(1-exp(-K*(a+aP-t0)))}
