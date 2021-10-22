#' Build DLMtool Stock object from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Stock DLMtool Stock object to start with
#' @param sc Scalar (multiplier) to compute upper and lower bounds of random uniform distribution from mean value
#' @param M_scLim Scalar for M_constant limits. Numeric vector of length 2
#' @param is_M_age_varying logical. Indicate if age varying M should be used. If TRUE, M and M2 will be upper bounds for age-varying M, set as a function of t.series$M from BAM rdat and M_scLim. If FALSE M will be a function of M.constant and M_scLim
#' @param steep_scLim Scalar for steep limits. Numeric vector of length 2
#' @param rec_sigma_scLim Scalar for rec_sigma limits. Numeric vector of length 2
#' @param rec_AC_scLim Scalar for rec_AC (lag-1 recruitment autocorrelation) limits. Numeric vector of length 2
#' @param Linf_scLim Scalar for Linf (Von Bertalanffy growth function) limits. Numeric vector of length 2
#' @param K_scLim Scalar for K (Von Bertalanffy growth function) limits. Numeric vector of length 2
#' @param t0_scLim Scalar for t0 (Von Bertalanffy growth function) limits. Numeric vector of length 2
#' @param len_cv_val_scLim Scalar for len_cv_val limits. Numeric vector of length 2
#' @param L_50_scLim Scalar for L_50 (length at 50 percent maturity) limits. Numeric vector of length 2
#' @param L_50_95_scLim Scalar for L_50_95 (Length increment from 50 percent to 95 percent maturity) limits. Numeric vector of length 2
#' @param D_scLim Scalar for D_95 limits. Numeric vector of length 2
#' @param Msd see \code{\link[DLMtool]{Stock-class}}
#' @param Ksd  see \code{\link[DLMtool]{Stock-class}}
#' @param Linfsd see \code{\link[DLMtool]{Stock-class}}
#' @param length_sc Scalar (multiplier) to convert length units. MSEtool examples seem to use cm whereas BAM uses mm.
#' @param SRrel see \code{\link[DLMtool]{Stock-class}}
#' @param Size_area_1 see \code{\link[DLMtool]{Stock-class}}
#' @param Frac_area_1 see \code{\link[DLMtool]{Stock-class}}
#' @param Prob_staying see \code{\link[DLMtool]{Stock-class}}
#' @param R0 "The magnitude of unfished recruitment. Single value. Positive real number"  (\code{\link[DLMtool]{Stock-class}})
#' @param use_bam_R0 Should the value of R0 from BAM be used? Note that units may vary (e.g. eggs or Age-1)). logical
#' @param AC Default autocorrelation for rec devs (in log space)
#' @param use_bam_AC Should recruitment autocorrelation be computed from BAM rec devs? If rec devs are not available, default AC value is used. logical
#' @param Fdisc  see \code{\link[DLMtool]{Stock-class}}. This function will try to get a range of values from "D.mort" values in the parms object, but if it can't it sets the value to numeric(0).
#' @param Mat_age1_max Limit maximum value of proportion mature of first age class (usually age-1). Models sometimes fail when maturity of first age class is too high (e.g. >0.5)
#' @param herm Is the species hermaphroditic? If "gonochoristic", use female maturity. If "protogynous", use a function of male and female maturity.
#' @param genus_species Genus and species names separated by a space (e.g. "Centropristis striata").
#' @keywords bam stock assessment fisheries DLMtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Build and plot DLMtool Stock (Stock-class object)
#' Stock_RedPorgy <- rdat_to_Stock(rdat_RedPorgy)
#' plot(Stock_RedPorgy)
#'
#' # Build DLMtool operating model (OM-class object)
#' OM_RedPorgy <- new("OM", Stock_RedPorgy, Generic_Fleet, Precise_Unbiased, Perfect_Imp)
#' # Run and plot simple management strategy evaluation (MSE)
#' mse_out <- runMSE(OM_RedPorgy)
#' NOAA_plot(mse_out)
#' }

rdat_to_Stock <- function(
  rdat, Stock = new('Stock'),
  sc = 0,  scLim = sc*c(-1,1)+1,
  M_scLim = 0.001*c(-1,1)+1,
  is_M_age_varying = FALSE,
  steep_scLim = scLim, rec_sigma_scLim = scLim, rec_AC_scLim = scLim,
  Linf_scLim = scLim, K_scLim = scLim, t0_scLim = scLim, len_cv_val_scLim = scLim, L_50_scLim = scLim,
  L50_95_scLim = scLim, D_scLim  = scLim,
  Msd = c(0,0), Ksd = c(0,0), Linfsd = c(0,0), length_sc=0.1,
  Size_area_1 = c(0.5,0.5), Frac_area_1 = c(0.5,0.5), Prob_staying = c(0.5,0.5),
  SRrel = 1, R0 = 1000, use_bam_R0 = TRUE,
  AC = 0.2, use_bam_AC = TRUE,
  Fdisc=NULL,
  Mat_age1_max = 0.49,
  herm = NULL, genus_species = NULL
){

rdat <- standardize_rdat(rdat)

info <- rdat$info
parms <- rdat$parms
parm.cons <- rdat$parm.cons
a.series <- rdat$a.series
t.series <- rdat$t.series

Name <- gsub(" ","",str_to_title(info$species))
years <- paste(parms$styr:parms$endyr)
nyears <- length(years)

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

t.series <- t.series[years,]

Common_Name <- str_replace_all(Name,"(?<=[a-z])(?=[A-Z])"," ")
if(is.null(genus_species)){genus_species <- bamStockMisc[Name,"Species"]}
if(is.null(herm)){herm <- bamStockMisc[Name,"herm"]}

Linf <- parm.cons$Linf[1]*length_sc
K <- parm.cons$K[1]
t0 <- parm.cons$t0[1]

# Compute R0 in units of catch (klb)
Nage_F0 <- expDecay(age=age,Z=a.series$M,N0=1)
bam_R0 <- parms$BH.R0
bam_age_R <- min(rdat$a.series$age)
N0_F0 <- bam_R0*(Nage_F0/Nage_F0[paste(bam_age_R)]) # Numbers at age zero with no fishing


R0 <- ifelse(use_bam_R0, "yes"=N0_F0, "no" =R0)

# Set slot values
slot(Stock,"Name") <- Name
slot(Stock,"Common_Name") <- Common_Name
slot(Stock,"Species") <- genus_species
slot(Stock,"maxage") <- max(a.series$age)
slot(Stock,"R0") <- R0
if(is_M_age_varying){
slot(Stock,"M") <-  a.series$M*M_scLim[1] # lower bound of age-dependent M
slot(Stock,"M2") <- a.series$M*M_scLim[2] # upper bound of age-dependent M
}else{
  slot(Stock,"M") <- parms[["M.constant"]]*M_scLim # age-dependent M
}
slot(Stock,"Msd") <- Msd
slot(Stock,"h") <- local({
  a <- parm.cons$steep[1]*steep_scLim
  pmax(pmin(a,0.99),0.2) # Constrain to be between 0.2 and 0.99
})
slot(Stock,"SRrel") <- SRrel
slot(Stock,"Perr") <- parm.cons$rec_sigma[1]*rec_sigma_scLim

rec_AC <- local({
  logR.dev <- t.series$logR.dev
  if(all(logR.dev==0)){
    out <- AC}else{
      out <- acf(logR.dev,plot=FALSE)$acf[2,,1] # lag-1 autocorrelation
    }
  out
})
slot(Stock,"AC") <- round(rec_AC*rec_AC_scLim,3) # Upper and lower limits

slot(Stock,"a") <- parms$wgt.a/length_sc^parms$wgt.b # Adjust a parameter for length units
slot(Stock,"b") <- parms$wgt.b

slot(Stock,"Linf") <-  Linf*Linf_scLim
slot(Stock,"K") <- K*K_scLim
slot(Stock,"t0") <- t0*t0_scLim
slot(Stock,"LenCV") <- parm.cons$len_cv_val[1]*len_cv_val_scLim
slot(Stock,"Ksd") <- Ksd
slot(Stock,"Linfsd") <- Linfsd

slot(Stock,"Size_area_1")  <-  Size_area_1
slot(Stock,"Frac_area_1")  <-  Frac_area_1
slot(Stock,"Prob_staying") <-  Prob_staying

# Compute proportion mature at age
pmat <- pmatage(a.series=a.series,Mat_age1_max=Mat_age1_max,herm=herm,age=age)$pmat

# Compute maturity-at-length L50 and L50_95
mat_at_len <- local({
  # Predict proportion mature at from linear interpolation
  age_pr <- seq(min(age),max(age),length=1000)
  pmat_pr <- approx(age,pmat,xout = age_pr)$y

  age50 <- age_pr[which.min(abs(pmat_pr-0.50))] # age at 50% maturity
  age95 <- age_pr[which.min(abs(pmat_pr-0.95))] # age at 95% maturity

  len50 <- vb_len(Linf=Linf, K=K, t0=t0, a=age50) # length at 50% maturity
  len95 <- vb_len(Linf=Linf, K=K, t0=t0, a=age95) # length at 95% maturity
  return(list("L50"=len50,"L50_95"=len95-len50))
})
L50 <- mat_at_len$L50
A50 <- vb_age(L=L50,Linf=Linf,K=K,t0=t0)

slot(Stock,"L50") <- mat_at_len$L50*L_50_scLim
slot(Stock,"L50_95") <- mat_at_len$L50_95*L50_95_scLim

# Compute current level of depletion
SSBend <- parms$SSBmsy*parms$SSBend.SSBmsy
D <- SSBend/parms$SSB0

slot(Stock,"D") <- D*D_scLim
slot(Stock,"Fdisc") <- local({
  if(is.null(Fdisc)){
  a <- grepl("D.mort",names(parms))
  if(any(a)){
  range(unlist(parms[grepl("D.mort",names(parms))]))
  }else{
    numeric(0)
  }
  }else{
    Fdisc
  }
})

slot(Stock,"Source") <- paste0(paste(unlist(info[c("species","title","date")]),collapse = "; "),"; rdat file")

return(Stock)
}
