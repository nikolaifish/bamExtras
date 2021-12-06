#' Build MSEtool Stock object from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Fleet MSEtool Fleet object to start with
#' @param length_sc Scalar (multiplier) to convert length units. MSEtool examples seem to use cm whereas BAM uses mm.
#' @param Eff_sc Scale Effort vector to compute EffLower and EffUpper vectors. Numeric vector of length 2.
#' @param Spat_targ see \code{\link[MSEtool]{Fleet-class}}
#' @param Esd see \code{\link[MSEtool]{Fleet-class}}. The default for most MSEtool built-in Fleet objects is c(0.1,0.4)
#' @param qinc see \code{\link[MSEtool]{Fleet-class}}
#' @param qcv see \code{\link[MSEtool]{Fleet-class}}
#' @param L5_sc Scalar for L5 (see \code{\link[MSEtool]{Fleet-class}}). Numeric vector of length 2
#' @param LFS_sc Scalar for LFS (see \code{\link[MSEtool]{Fleet-class}}). Numeric vector of length 2
#' @param LR5_sc Scalar for LR5 (see \code{\link[MSEtool]{Fleet-class}}). Numeric vector of length 2
#' @param LFR_sc Scalar for LFR (see \code{\link[MSEtool]{Fleet-class}}). Numeric vector of length 2
#' @param Vmaxlen see \code{\link[MSEtool]{Fleet-class}}
#' @param Rmaxlen see \code{\link[MSEtool]{Fleet-class}}
#' @param DR see \code{\link[MSEtool]{Fleet-class}}
#' @param isRel see \code{\link[MSEtool]{Fleet-class}}
#' @keywords bam stock assessment fisheries MSEtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Build MSEtool Stock and Fleet (Stock-class and Fleet-class objects)
#' Stock_RedPorgy <- rdat_to_Stock(rdat_RedPorgy)
#' Fleet_RedPorgy <- rdat_to_Fleet(rdat_RedPorgy)
#'
#' # Build MSEtool operating model (OM-class object)
#' OM_RedPorgy <- new("OM", Stock_RedPorgy, Fleet_RedPorgy, Precise_Unbiased, Perfect_Imp)
#' # Run and plot simple management strategy evaluation (MSE)
#' mse_out <- runMSE(OM_RedPorgy)
#' NOAA_plot(mse_out)
#' }

rdat_to_Fleet <- function(
  rdat,
  Fleet = new('Fleet'),
  length_sc=0.1,
  Eff_sc = 0.01*c(-1,1)+1,
  Spat_targ = c(1,1),
  Esd = c(0,0),
  qinc = c(0,0),
  qcv = c(0,0),
  L5_sc = c(1,1),
  LFS_sc = c(1,1),
  LR5_sc = c(1,1),
  LFR_sc = c(1,1),
  Vmaxlen = c(1,1),
  Rmaxlen = c(1,1),
  DR = NULL,
  isRel = FALSE
)
{

rdat <- standardize_rdat(rdat)

info <- rdat$info
parms <- rdat$parms
parm.cons <- rdat$parm.cons
years <- paste(parms$styr:parms$endyr)
nyears <- length(years)
a.series <- rdat$a.series
t.series <- rdat$t.series[years,]

age <- a.series$age

Name <- gsub(" ","",str_to_title(info$species))

## Effort
# E=C/qB # Effort calculation from Haddon 2011 page 290 Eq. 11.4
# E=C/B  # If we ignore units we can just use this
E <- local({
  C <- t.series$total.L.klb
  B <- t.series$B # Should be in metric tons
  Bklb <- measurements::conv_unit(B,"metric_ton","lbs")/1000
  E <- C/Bklb
  names(E) <- years
  E
})
# E_scaled <- (E/max(E,na.rm=TRUE))
EffLower <- E*Eff_sc[1]
EffUpper <- E*Eff_sc[2]

# Set slot values
slot(Fleet,"Name") <- Name
slot(Fleet,"nyears") <- nyears
slot(Fleet,"CurrentYr") <- parms$endyr

slot(Fleet,"EffYears") <- 1:nyears # Should be a sequence of integers starting with 1, not actual years
slot(Fleet,"EffLower") <- EffLower
slot(Fleet,"EffUpper") <- EffUpper
slot(Fleet,"Esd") <- Esd

slot(Fleet,"qinc") <- qinc
slot(Fleet,"qcv") <- qcv

# Estimate length at 5% and full (100%) vulnerability (total selectivity) and retention (landings selectivity)
vuln_out <- vulnerability(Vdata = rdat$sel.age$sel.v.wgted.tot, retdata = rdat$sel.age$sel.v.wgted.L, Linf = parm.cons$Linf[8], K = parm.cons$K[8], t0 = parm.cons$t0[8], length_sc=length_sc)

slot(Fleet,"L5") <- vuln_out$L5*L5_sc
slot(Fleet,"LFS") <- vuln_out$LFS*LFS_sc

slot(Fleet,"LR5") <- vuln_out$LR5*LR5_sc
slot(Fleet,"LFR") <- vuln_out$LFR*LFR_sc

slot(Fleet,"isRel") <- isRel

slot(Fleet,"Vmaxlen") <- Vmaxlen
slot(Fleet,"Rmaxlen") <- Rmaxlen

if(is.null(DR)){
  DR <- as.numeric(sort(1-quantile(pmin(1,(vuln_out$retdata/vuln_out$Vdata)[paste(round(vuln_out$AFS):max(age))]),probs=c(0.25,0.75))))
}

slot(Fleet,"DR") <- DR

slot(Fleet,"Spat_targ") <- Spat_targ

return(Fleet)
}
