#' Change non-standard names of rdat objects to standard naming conventions. Compute standard objects.
#'
#' @param rdat rdat (list) object read in with dget()
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' names(rdat_RedPorgy$a.series)
#' rdat_RedPorgy_std <- standardize_rdat(rdat_RedPorgy)
#' names(rdat_RedPorgy_std$a.series)
#' }

standardize_rdat <- function(rdat,
                           parms.key=c("R0"="BH.R0"),   
                           a.series.key=c("mat.male.endyr"="mat.male",
                                          "mat.fem.endyr"="mat.female",
                                          "prop.female.endyr"="prop.female"
                                          ),
                           t.series.key=c("total.L.wgt.klb"="total.L.klb")
                           ){
  # info
  info$species <- gsub("SA ","",info$species)
  
  # parms
  parms <- rdat$parms
  parms.names <- names(parms)
  names(parms) <- find_replace(names(parms.key),parms.key,parms.names)
  rdat$parms <- parms
  
  # a.series
  a.series <- rdat$a.series
  a.series.names <- names(a.series)
  names(a.series) <- find_replace(names(a.series.key),a.series.key,a.series.names)
  
  # If a.series includes prop.male instead of prop.female, compute prop.female
  if(!"prop.female"%in%names(a.series)&"prop.male"%in%names(a.series)){
    a.series$prop.female <- 1-a.series$prop.male
  }
  rdat$a.series <- a.series
  
  # t.series
  t.series <- rdat$t.series
  t.series.names <- names(t.series)
  names(t.series) <- find_replace(names(t.series.key),t.series.key,t.series.names)
  rdat$t.series <- t.series
  
  return(rdat)
}