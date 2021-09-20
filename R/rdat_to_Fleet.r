#' Build DLMtool Stock object from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Fleet DLMtool Fleet object to start with
#' @param Eff_sc Scale Effort vector to compute EffLower and EffUpper vectors. Numeric vector of length 2.
#' @param Spat_targ see \code{\link[DLMtool]{Fleet-class}}
#' @param Esd see \code{\link[DLMtool]{Fleet-class}}. The default for most MSEtool built-in Fleet objects is c(0.1,0.4)
#' @param qinc see \code{\link[DLMtool]{Fleet-class}}
#' @param qcv see \code{\link[DLMtool]{Fleet-class}}
#' @param L5_sc Scalar for L5 (see \code{\link[DLMtool]{Fleet-class}}). Numeric vector of length 2
#' @param LFS_sc Scalar for LFS (see \code{\link[DLMtool]{Fleet-class}}). Numeric vector of length 2
#' @param LR5_sc Scalar for LR5 (see \code{\link[DLMtool]{Fleet-class}}). Numeric vector of length 2
#' @param LFR_sc Scalar for LFR (see \code{\link[DLMtool]{Fleet-class}}). Numeric vector of length 2
#' @param Vmaxlen see \code{\link[DLMtool]{Fleet-class}}
#' @param Rmaxlen see \code{\link[DLMtool]{Fleet-class}}
#' @param isRel see \code{\link[DLMtool]{Fleet-class}}
#' @keywords bam stock assessment fisheries DLMtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Build DLMtool Stock and Fleet (Stock-class and Fleet-class objects)
#' Stock_RedPorgy <- rdat_to_Stock(rdat_RedPorgy)
#' Fleet_RedPorgy <- rdat_to_Fleet(rdat_RedPorgy)
#'
#' # Build DLMtool operating model (OM-class object)
#' OM_RedPorgy <- new("OM", Stock_RedPorgy, Fleet_RedPorgy, Precise_Unbiased, Perfect_Imp)
#' # Run and plot simple management strategy evaluation (MSE)
#' mse_out <- runMSE(OM_RedPorgy)
#' NOAA_plot(mse_out)
#' }

rdat_to_Fleet <- function(
  rdat,
  Fleet = new('Fleet'),
  Eff_sc = 0.1*c(-1,1)+1,
  Spat_targ = c(1,1),
  Esd = c(0.1,0.4),
  qinc = c(0,0),
  qcv = c(0,0),
  L5_sc = c(1,1),
  LFS_sc = c(1,1),
  LR5_sc = c(1,1),
  LFR_sc = c(1,1),
  Vmaxlen = c(1,1),
  Rmaxlen = c(1,1),
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

Name <- gsub(" ","",str_to_title(info$species))

## Effort
# E=C/qB # Effort calculation from Haddon 2011 page 290 Eq. 11.4
# E=C/B  # If we ignore units we can just use this
E <- local({
  C <- t.series$total.L.klb
  B <- t.series$B
  E <- C/B
  names(E) <- years
  E
})

# Set slot values
slot(Fleet,"Name") <- Name
slot(Fleet,"nyears") <- nyears
slot(Fleet,"Spat_targ") <- Spat_targ

slot(Fleet,"EffYears") <- 1:nyears # Should be a sequence of integers starting with 1, not actual years
slot(Fleet,"EffLower") <- E*Eff_sc[1]
slot(Fleet,"EffUpper") <- E*Eff_sc[2]
slot(Fleet,"Esd") <- Esd

slot(Fleet,"qinc") <- qinc
slot(Fleet,"qcv") <- qcv

# Estimate length at 5% and full (100%) selection (landings) and retention (discards)
sel_at_len <- local({
  # Selectivity
  seldata <- rdat$sel.age$sel.v.wgted.L
  seldata <- seldata/max(seldata) # Scaled as a proportion of maximum selectivity

  age <- as.numeric(names(seldata))
  selage <- seldata

  age_pr <- seq(min(age),max(age),length=1000)
  sel_pr <- approx(age,selage,xout = age_pr)$y

  agesel_05 <- age_pr[which.min(abs(sel_pr-0.05))] # age at 5% selection
  agesel_100 <- age_pr[which.min(abs(sel_pr-1.00))] # age at 100% selection

  Linf <- parm.cons$Linf[1]
  K <- parm.cons$K[1]
  t0 <- parm.cons$t0[1]
  lenage <- vb_len(Linf=Linf, K=K, t0=t0, a=age)*length_sc
  len_pr <- vb_len(Linf=Linf, K=K, t0=t0, a=age_pr)*length_sc

  lensel_05 <- vb_len(Linf=Linf, K=K, t0=t0, a=agesel_05)*length_sc # Length at 5% selection
  lensel_100 <- vb_len(Linf=Linf, K=K, t0=t0, a=agesel_100)*length_sc # Length at 100% selection

  # Retention
  retdata <- rdat$sel.age$sel.v.wgted.D
  if(is.null(retdata)){
    lenret_05 <- L5
    lenret_100 <- LFS
  }else{
  retdata <- retdata/max(retdata) # Scaled as a proportion of maximum selectivity

  age <- as.numeric(names(retdata))
  retage <- retdata

  age_pr <- seq(min(age),max(age),length=1000)
  ret_pr <- approx(age,retage,xout = age_pr)$y

  ageret_05 <- age_pr[which.min(abs(ret_pr-0.05))] # age at 5% discard selection
  ageret_100 <- age_pr[which.min(abs(ret_pr-1.00))] # age at 100% discard selection

  Linf <- parm.cons$Linf[1]
  K <- parm.cons$K[1]
  t0 <- parm.cons$t0[1]
  lenage <- vb_len(Linf=Linf, K=K, t0=t0, a=age)*length_sc
  len_pr <- vb_len(Linf=Linf, K=K, t0=t0, a=age_pr)*length_sc

  lenret_05 <- vb_len(Linf=Linf, K=K, t0=t0, a=ageret_05)*length_sc # Length at 5% discard selection
  lenret_100 <- vb_len(Linf=Linf, K=K, t0=t0, a=ageret_100)*length_sc # Length at 100% discard selection
  }
  return(list("L5"=lensel_05,"LFS"=lensel_100,"LR5"=lenret_05,"LFR"=lenret_100))
})

slot(Fleet,"L5") <- sel_at_len$L5*L5_sc
slot(Fleet,"LFS") <- sel_at_len$LFS*LFS_sc

slot(Fleet,"LR5") <- sel_at_len$LR5*LR5_sc
slot(Fleet,"LFR") <- sel_at_len$LFR*LFR_sc

slot(Fleet,"isRel") <- isRel

slot(Fleet,"Vmaxlen") <- Vmaxlen
slot(Fleet,"Rmaxlen") <- Rmaxlen

slot(Fleet,"CurrentYr") <- parms$endyr

return(Fleet)
}
