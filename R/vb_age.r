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
#' \dontrun{
#' # Use data for Black Sea Bass
#' rdat <- rdat_BlackSeaBass
#' aser <- rdat$a.series
#' age <- aser$age
#' len <- aser$length
#' parest <- rdat$parm.cons[8,]
#'
#' # Convert age to length
#' len_out <- with(parest,{vb_len(a=age,Linf=Linf,K=K,t0=t0)})
#' len_out
#'
#' # Convert length to age
#' age_out <- with(parest,{vb_age(L=len,Linf=Linf,K=K,t0=t0)})
#' age_out
#'
#' # Note that vb_age is not quite the inverse of vb_len unless t0 = 0
#' # Convert length to age
#' with(parest,{vb_age(L=with(parest,{vb_len(a=age,Linf=Linf,K=K,t0=0)}),Linf=Linf,K=K,t0=0)})

vb_age <- function(L,Linf,K,t0) {(log(1-L/Linf)/-K)-t0}
