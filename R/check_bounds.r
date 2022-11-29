#' Check to see if any estimated parameters are near bounds
#'
#' Check to see if any estimated parameters are near bounds
#' @param x a BAM rdat object, i.e. x <- dget(rdat)
#' @param p_cutoff set minimum acceptable value for how close to a bound the estimated parameter can be before considering it too close to a bound, as a proportion of the range of bounds supplied to BAM
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' check_bounds(rdat_BlackSeaBass)
#' }
#'

check_bounds <- function(x,p_cutoff=0.01){
  pc <- x$parm.cons # parameter matrix
  pc <- pc[,pc[4,]>0] # only consider parameters that are being estimated (not fixed)
  pc.range <- pc[3,]-pc[2,] # parameter range (width of bounds)
  pc.mindist <- apply(abs(pc[2:3,,drop=FALSE]-pc[c(8,8),,drop=FALSE]),2,min) # distance from estimate to nearest bound
  pc.pdist <- pc.mindist/pc.range # distance to nearest bound as a proportion of bound range
  pc.nearbound <- pc.pdist<=p_cutoff
  return(list("p_cutoff"=p_cutoff,"pdist"=pc.pdist,"nearbound"=pc.nearbound))
}
