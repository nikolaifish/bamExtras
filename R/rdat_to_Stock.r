#' Build DLMtool Stock object from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Stock DLMtool Stock object to start with
#' @param sc Scalar (multiplier) to compute upper and lower bounds of random uniform distribution from mean value
#' @param M_constant_sc Scalar for M_constant. Numeric vector of length 2
#' @param steep_sc Scalar for steep. Numeric vector of length 2
#' @param rec_sigma_sc Scalar for rec_sigma. Numeric vector of length 2
#' @param rec_AC_sc Scalar for rec_AC (lag-1 recruitment autocorrelation). Numeric vector of length 2
#' @param Linf_sc Scalar for Linf (Von Bertalanffy growth function). Numeric vector of length 2
#' @param K_sc Scalar for K (Von Bertalanffy growth function). Numeric vector of length 2
#' @param t0_sc Scalar for t0 (Von Bertalanffy growth function). Numeric vector of length 2
#' @param len_cv_val_sc Scalar for len_cv_val. Numeric vector of length 2
#' @param L_50_sc Scalar for L_50 (length at 50 percent maturity). Numeric vector of length 2
#' @param L_50_95_sc Scalar for L_50_95 (Length increment from 50 percent to 95 percent maturity). Numeric vector of length 2
#' @param D_sc Scalar for D_95. Numeric vector of length 2
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
  M_constant_sc = 0.1*c(-1,1)+1, steep_sc = scLim, rec_sigma_sc = scLim, rec_AC_sc = scLim,
  Linf_sc = scLim, K_sc = scLim, t0_sc = scLim, len_cv_val_sc = scLim, L_50_sc = scLim,
  L50_95_sc = scLim, D_sc  = scLim,
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

# MSEtool expects age-based data to begin with age 1
  if(min(rdat$a.series$age)<1){
    warning(paste(Name,": Minimum age <1. Age-based data limited to age >=1"))
    a.series <- a.series[a.series$age%in%1:max(a.series$age),]
    rdat$a.series <- a.series
  }
  age <- rdat$a.series$age

t.series <- t.series[years,]

Common_Name <- str_replace_all(Name,"(?<=[a-z])(?=[A-Z])"," ")
if(is.null(genus_species)){genus_species <- bamStockMisc[Name,"Species"]}
if(is.null(herm)){herm <- bamStockMisc[Name,"herm"]}

Linf <- parm.cons$Linf[1]
K <- parm.cons$K[1]
t0 <- parm.cons$t0[1]

R0 <- ifelse(use_bam_R0, "yes"=parms$BH.R0, "no" =R0)

# Set slot values
slot(Stock,"Name") <- Name
slot(Stock,"Common_Name") <- Common_Name
slot(Stock,"Species") <- genus_species
slot(Stock,"maxage") <- max(a.series$age)
slot(Stock,"R0") <- R0
slot(Stock,"M") <-  a.series$M*M_constant_sc[1] # lower bound of age-dependent M
slot(Stock,"M2") <- a.series$M*M_constant_sc[2] # upper bound of age-dependent M
slot(Stock,"Msd") <- Msd
slot(Stock,"h") <- local({
  a <- parm.cons$steep[1]*steep_sc
  pmax(pmin(a,0.99),0.2) # Constrain to be between 0.2 and 0.99
})
slot(Stock,"SRrel") <- SRrel
slot(Stock,"Perr") <- parm.cons$rec_sigma[1]*rec_sigma_sc

rec_AC <- local({
  logR.dev <- t.series$logR.dev
  if(all(logR.dev==0)){
    out <- AC}else{
      out <- acf(logR.dev,plot=FALSE)$acf[2,,1] # lag-1 autocorrelation
    }
  out
})
slot(Stock,"AC") <- round(rec_AC*rec_AC_sc,3) # Upper and lower limits

slot(Stock,"a") <- parms$wgt.a/length_sc^parms$wgt.b # Adjust a parameter for length units
slot(Stock,"b") <- parms$wgt.b

slot(Stock,"Linf") <-  Linf*Linf_sc*length_sc
slot(Stock,"K") <- K*K_sc
slot(Stock,"t0") <- t0*t0_sc
slot(Stock,"LenCV") <- parm.cons$len_cv_val[1]*len_cv_val_sc
slot(Stock,"Ksd") <- Ksd
slot(Stock,"Linfsd") <- Linfsd

slot(Stock,"Size_area_1")  <-  Size_area_1
slot(Stock,"Frac_area_1")  <-  Frac_area_1
slot(Stock,"Prob_staying") <-  Prob_staying

# Compute proportion mature at age
pmat <- pmatage(rdat=rdat,Mat_age1_max=Mat_age1_max,herm=herm,age=age)$pmat

# Compute maturity-at-length L50 and L50_95
mat_at_len <- local({
  # Predict proportion mature at from linear interpolation
  age_pr <- seq(min(age),max(age),length=1000)
  pmat_pr <- approx(age,pmat,xout = age_pr)$y

  age50 <- age_pr[which.min(abs(pmat_pr-0.50))] # age at 50% maturity
  age95 <- age_pr[which.min(abs(pmat_pr-0.95))] # age at 95% maturity

  len50 <- vb_len(Linf=Linf, K=K, t0=t0, a=age50)*length_sc # length at 50% maturity
  len95 <- vb_len(Linf=Linf, K=K, t0=t0, a=age95)*length_sc # length at 95% maturity
  return(list("L50"=len50,"L50_95"=len95-len50))
})
L50 <- mat_at_len$L50
A50 <- vb_age(Linf,K,t0,L50)

slot(Stock,"L50") <- mat_at_len$L50*L_50_sc
slot(Stock,"L50_95") <- mat_at_len$L50_95*L50_95_sc

# Compute current level of depletion
SSBend <- parms$SSBmsy*parms$SSBend.SSBmsy
D <- SSBend/parms$SSB0

slot(Stock,"D") <- D*D_sc
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
