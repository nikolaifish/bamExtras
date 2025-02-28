#' logistic-exponential function
#' logistic-exponential selectivity function. The function has 4 parameters.
#' When fitting in bam, one of the parameters should be fixed.
#' Note that if the minimim value of \code{ages} is equal to the value of \code{joint},
#' then \code{A50} and \code{slope} have no effect on the result.
#' @param ages vector of ages
#' @param A50  age at 50% selection (ascending limb)
#' @param slope rate of increase
#' @param sigma controls rate of descent (descending limb)
#' @param joint age to join curves
#' @keywords bam stock assessment fisheries population dynamics
#' @author Erik Williams, Kyle Shertzer, and Nikolai Klibansky
#' @export
#' @examples
#'
#' rdat <- rdat_GrayTriggerfish
#' fp <- rdat$parm.cons[8,grepl("rHB\\.D$",names(rdat$parm.cons))]
#' logistic_exponential(ages=rdat$a.series$age,
#' A50 = fp$selpar.L50.rHB.D,
#' slope = fp$selpar.slope.rHB.D,
#' sigma = fp$selpar.sigma.rHB.D,
#' joint = fp$selpar.afull.rHB.D)
logistic_exponential <- function(ages, A50, slope, sigma, joint){
  sel <- rep(1,length(ages)) # this is the value at the joint
  sela <- 1/(1+exp(-1*slope*(ages-A50)))
  selb <- exp(-1.*((ages-joint)/sigma)^2)
  sel[which(ages<joint)] <- sela[which(ages<joint)]
  sel[which(ages>joint)] <- selb[which(ages>joint)]
  sel <- sel/max(sel) # I think that this should be unnecessary, since the max(sel) should be 1

  return(sel)
}
