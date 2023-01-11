#' L2F
#'
#' Approximate level of fishing mortality that achieves a given level of landings.
#' This function is modified from code that was part of the BAM projections. The
#' current function is intended for use prior to running run_proj to compute the
#' F that matches specific values of landings. This is particularly useful if
#' landings estimates corresponding to a year or years between the terminal year of the
#' current stock assessment and the start year of new management is available.
#' @param L_tar target landings value to match to the target fishing mortality value, F_tar. Must be in the same units as wgt_L.
#' @param F_est initial guesses for the F that matches L_tar
#' @param sel_L selectivity for landings
#' @param sel_Z selectivity for total mortality
#' @param M natural mortality rate, at age
#' @param N_a vector of numbers of fish, at age, in population
#' @param wgt_L vector of weight of fish, at age, in the landings
#' @param control A list of control parameters passed to \code{\link[stats]{nlminb}}
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' spp <- rdat_Tilefish
#' endyr <- paste(spp$parms$endyr)
#' # Compute landings in weight for the terminal year
#' L_endyr_klb <- L_calc(F_a = spp$F.age[endyr,], Z_a = spp$Z.age[endyr,],
#' N_a = spp$N.age[endyr,], wgt_L = spp$a.series$gutwgt.wgted.L.klb)
#' L2F(L_tar=L_endyr_klb, F_est=sum(spp$F.age[endyr,]),sel_L = spp$sel.age$sel.v.wgted.L,
#' sel_Z = spp$sel.age$sel.v.wgted.tot,M = spp$a.series$M, N_a = spp$N.age[endyr,], wgt_L = spp$a.series$gutwgt.wgted.L.klb)
#' }

L2F <- function(L_tar, F_est, sel_L, sel_Z , M, N_a, wgt_L, control=list(iter.max=500))
{
  # Objective function used in
  F_obj <- function(F_est, sel_L, sel_Z , M, N_a, wgt_L, L_tar){
    F_a=sel_L*F_est[1]
    Z_a=M+sel_Z*F_est[1]

    L_pre <- L_calc(F_a, Z_a, N_a, wgt_L)
    SS <- (L_pre-L_tar)^2
    return(SS)
  }

  F_fit <- nlminb(start=F_est, objective=F_obj, control=control,
                  sel_L=sel_L, sel_Z=sel_Z, M=M, N_a=N_a, wgt_L=wgt_L, L_tar=L_tar)
  F_tar <- F_fit$par[1]

  return(F_tar)
}
