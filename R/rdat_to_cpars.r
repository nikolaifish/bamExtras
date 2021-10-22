#' Build DLMtool cpars list from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param nsim Number of simulations in operating model
#' @param nyears Number of years of historical data
#' @param proyears Number of projection years
#' @param Mat_age1_max Limit maximum value of proportion mature of first age class (usually age-1).
#' @param herm Is the species hermaphroditic? If "gonochoristic", use female maturity. If "protogynous", use a function of male and female maturity.
#' @param V_scaleToMax Should vulnerability be rescaled to max vulnerability?
#' @param ret_scaleToMax Should retention be rescaled to max retention?
#' @keywords bam stock assessment fisheries DLMtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Build DLMtool Stock (Stock-class object)
#' Stock_RedPorgy <- rdat_to_Stock(rdat_RedPorgy)
#'
#' # Build DLMtool operating model (OM-class object)
#' OM_RedPorgy <- new("OM", Stock_RedPorgy, Generic_Fleet, Precise_Unbiased, Perfect_Imp)
#' # Copy OM
#' OM_RedPorgy1 <- OM_RedPorgy2 <- OM_RedPorgy
#' # Build cpars list
#' cpars <- rdat_to_cpars(rdat_RedPorgy, nsim=slot(OM_RedPorgy,"nsim"), nyears=slot(OM_RedPorgy,"nyears"), proyears=slot(OM_RedPorgy,"proyears"))
#'
#' slot(OM_RedPorgy2,"cpars") <- cpars
#'
#' ## Run and plot simple management strategy evaluations (MSE)
#' # Without cpars
#' mse_out1 <- runMSE(OM_RedPorgy1)
#' NOAA_plot(mse_out1)
#' # With cpars
#' mse_out2 <- runMSE(OM_RedPorgy2)
#' NOAA_plot(mse_out2)
#' }

rdat_to_cpars <- function(rdat,nsim,nyears=NULL,proyears=50,
                          Mat_age1_max = 0.49, herm = NULL,
                          V_scaleToMax=FALSE,
                          ret_scaleToMax=FALSE
                          ){
  rdat <- standardize_rdat(rdat)
  a.series <- rdat$a.series
  Name <- gsub(" ","",str_to_title(rdat$info$species))
  if(is.null(nyears)){
    nyears <- length(rdat$parms$styr:rdat$parms$endyr)
  }

  # MSEtool expects age-based data to begin with age 0
  if(min(a.series$age)>0){
    warning(paste(Name,": Minimum age > 0. Age-based data extrapolated to age-0"))
    a.series <- data_polate(a.series,xout=0:max(a.series$age))
    a.series <- data_lim(a.series,xlim=c(0,Inf))
    a.series <- data_lim(a.series,xname=c("prop.female","prop.male","mat.female","mat.male"),xlim=c(0,1))
    a.series <- as.data.frame(a.series)
    rownames(a.series) <- a.series$age
  }
  age <- a.series$age
  M <- a.series$M

  # M_ageArray
  M_ageArray_dim <- c(nsim,length(age),nyears+proyears)
  M_ageArray <- array(data = rep(M,each=nsim), dim = M_ageArray_dim,
                   dimnames=list("sim"=seq_len(M_ageArray_dim[1]),
                                 "M"=age,
                                 "year"=seq_len(M_ageArray_dim[3])))

  # Mat_age
  if(is.null(herm)){herm <- bamStockMisc[Name,"herm"]}
  # Compute proportion mature at age
  pmat <- pmatage(a.series=a.series,Mat_age1_max=Mat_age1_max,herm=herm,age=age)$pmat

  Mat_age_dim <- c(nsim,length(age),nyears+proyears)
  Mat_age <- array(data = rep(pmat,each=nsim), dim = Mat_age_dim,
                   dimnames=list("sim"=seq_len(Mat_age_dim[1]),
                                 "mat"=age,
                                 "year"=seq_len(Mat_age_dim[3])))

  # V (vulnerability at age; i.e. total selectivity)
  Vage <- rdat$sel.age$sel.v.wgted.tot
  if(min(as.numeric(names(Vage)))>0){
    warning(paste(Name,": Minimum age in sel.v.wgted.tot > 0. sel.v.wgted.tot extrapolated to age-0"))
  Vage <- pmin(1,pmax(0,Hmisc::approxExtrap(
    x=as.numeric(names(Vage)),y=Vage,xout=age)$y))
  names(Vage) <- paste(age)
  }
  if(V_scaleToMax){
  Vage <- Vage/max(Vage) # Scaled to a maximum of 1
  }
  Vage <- Vage[paste(age)]
  V_dim <- c(nsim,length(age),nyears+proyears)
  V <- array(data = rep(Vage,each=nsim), dim = V_dim,
             dimnames=list("sim"=seq_len(V_dim[1]),
                           "V"=age,
                           "year"=seq_len(V_dim[3])))

  # retA (retention at age; i.e. landings selectivity)
  retage <- rdat$sel.age$sel.v.wgted.L
  if(min(as.numeric(names(retage)))>0){
    warning(paste(Name,": Minimum age in sel.v.wgted.tot > 0. sel.v.wgted.tot extrapolated to age-0"))
    retage <- pmin(1,pmax(0,Hmisc::approxExtrap(
      x=as.numeric(names(retage)),y=retage,xout=age)$y))
    names(retage) <- paste(age)
  }
  if(ret_scaleToMax){
  retage <-  retage/max(retage) # Scaled to a maximum of 1
  }
  retage <- retage[paste(age)]
  ret_dim <- c(nsim,length(age),nyears+proyears)
  retA <- array(data = rep(retage,each=nsim), dim = ret_dim,
             dimnames=list("sim"=seq_len(ret_dim[1]),
                           "ret"=age,
                           "year"=seq_len(ret_dim[3])))

  # Build cpars
  cpars <- list()
  cpars$M_ageArray <- M_ageArray
  cpars$Mat_age <- Mat_age
  cpars$V <- V
  cpars$retA <- retA

  return(cpars)
}
