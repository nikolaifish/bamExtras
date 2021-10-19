#' Compute age-varying natural mortality (M) based on Charnov et al (2013) and scale based on constant M
#'
#' @param a ages at which to compute M
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param par_wl_a a parameter from weight~length conversion equation weight = par_wl_a*length^par_wl_b
#' @param par_wl_b b parameter from weight~length conversion equation weight = par_wl_a*length^par_wl_b
#' @param aP proportion of age (value between 0 and 1) at which to compute length (e.g. aP=0.5 to compute length at midyear)
#' @param M_constant constant M value used to scale age-varying M. Defaults to NULL. M will be scaled if a numeric value of M_constant is supplied.
#' @param aMin minimum age to include in computation of scaling factor when scaling M
#' @param par_a a parameter for Lorenzen equation par_a*W^par_b. Defaults to value from manuscript par_a=3.69.
#' @param par_b b parameter for Lorenzen equation par_a*W^par_b. Defaults to value from manuscript par_b=-0.305.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @references Lorenzen, K. 1996. The relationship between body weight and natural mortality in juvenile and adult fish: a comparison of natural ecosystems and aquaculture. Journal of Fish Biology 49:627-642.
#' @export
#' @examples
#' \dontrun{
#' # M for Black Sea Bass
#' M_Lorenzen(a=0:11,Linf = 502, K = 0.173, t0 = -0.97, M_constant=0.38, aMin=2, par_wl_a=5.02e-05, par_wl_b=2.77)
#' }
#'
M_Lorenzen <- function(a,
                       Linf,K,t0,par_wl_a=1, par_wl_b=3, aP=0,
                       M_constant=NULL, aMin=0, par_a=3.69, par_b=-0.305)
{
  L <- vb_len(a=a,Linf=Linf,K=K,t0=t0,aP=aP)
  W <- par_wl_a*L^par_wl_b
  M_Lz <- par_a*W^par_b
  if(is.null(M_constant)){
    M.out <- M_Lz   # unscaled Lorenzen estimates
  }else{
    M_Lzsum <- sum(M_Lz[which(a>=aMin)])
    na_Sc <- length(which(a>=aMin))
    M_LzSc <- M_Lz*(M_constant*na_Sc)/M_Lzsum
    M_out <- M_LzSc  # scaled Lorenzen estimates
  }
  names(M_out) <- paste(a)
  return(M_out)
}
