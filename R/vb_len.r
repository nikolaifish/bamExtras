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
#'}
#'
vb_len <- function(a,Linf,K,t0,aP=0)  {Linf*(1-exp(-K*(a+aP-t0)))}
