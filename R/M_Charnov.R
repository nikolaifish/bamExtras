#' Compute age-varying natural mortality (M) based on Charnov et al (2013) and scale based on constant M
#'
#' @param a ages at which to compute M
#' @param Linf length infinity
#' @param K growth coefficient
#' @param t0 age (time) at length zero
#' @param aP proportion of age (value between 0 and 1) at which to compute length (e.g. aP=0.5 to compute length at midyear)
#' @param M_constant constant M value used to scale age-varying M. Defaults to NULL. M will be scaled if a numeric value of M_constant is supplied.
#' @param aMin minimum age to include in computation of scaling factor when scaling M
#' @param par_b b parameter for Charnov equation ((L/Linf)^par_b)*K. Defaults to value from manuscript par_b=-1.5.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @references Charnov, E. L., H. Gislason, and J. G. Pope. 2013. Evolutionary assembly rules for fish life histories. Fish and Fisheries 14:213-224.
#' @export
#' @examples
#' \dontrun{
#' # M for Black Sea Bass
#' M_Charnov(a=0:11,Linf = 502, K = 0.173, t0 = -0.97, M_constant=0.38, aMin=2)
#' }

M_Charnov <- function(a,Linf,K,t0,aP=0,M_constant=NULL,aMin=0,par_b=-1.5){
  L <- vb_len(a=a,Linf=Linf,K=K,t0=t0,aP=aP)
  M_Ch <- ((L/Linf)^par_b)*K
  if(is.null(M_constant)){
    M_out <- M_Ch    # unscaled Charnov estimates
  }else{
    M_ChMean <- mean(M_Ch[which(a>=aMin)])
    M_ChSc <- M_Ch*M_constant/M_ChMean
    M_out <- M_ChSc  # scaled Charnov estimates
  }
  names(M_out) <- paste(a)
  return(M_out)
}
