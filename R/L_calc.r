#' L_calc
#'
#' Compute landings in numbers or weight with Baranov equation. Set wgt_L = NULL to compute numbers instead of weight. Used in run_proj.
#' @param F vector of fishing mortality rate, at age
#' @param Z vector of total mortality rate at age including dead discards
#' @param N vector of numbers of fish, at age, in population
#' @param wgt_L vector of weight of fish, at age, in the landings. units of the result will be in units of wgt_L (e.g. metric tons, 1000 lb)
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' spp <- rdat_Tilefish
#' endyr <- paste(spp$parms$endyr)
#'
#' # Compute landings in numbers for the terminal year
#' sum(L_calc(F = spp$F.age[endyr,], Z = spp$Z.age[endyr,], N = spp$N.age[endyr,]))
#' # The value should be very close to this:
#' sum(spp$CLD.est.mats$Ln.total[endyr,])
#' # Compute landings in weight for the terminal year
#' sum(L_calc(F = spp$F.age[endyr,], Z = spp$Z.age[endyr,], N = spp$N.age[endyr,], wgt_L = spp$a.series$gutwgt.wgted.L.klb))
#' # The value should be very close to this:
#' sum(spp$CLD.est.mats$Lw.total[endyr,])
#' }

L_calc <-  function(F, Z, N, wgt_L=NULL){
  if(is.null(wgt_L)){wgt_L <- rep(1,length(N))}
  wgt_L*F*N*(1.0-exp(-Z))/Z
}
