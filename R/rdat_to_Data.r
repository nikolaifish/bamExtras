#' Build DLMtool Data object from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Data DLMtool Data object to start with
#' @param herm Is the species hermaphroditic? If "gonochoristic", use female maturity. If "protogynous", use a function of male and female maturity.
#' @param nsim see \code{\link[MSEtool]{OM-class}}
#' @param genus_species Genus and species names separated by a space (e.g. "Centropristis striata").
#' @param Region see \code{\link[MSEtool]{Data-class}}
#' @param Fref_name Name of F reference as named in rdat$parms list (e.g. Fmsy, F30)
#' @param Rec see \code{\link[MSEtool]{Data-class}}. Set to "bam_recruits" to set this value to the vector of recruits estimated by bam for years where rec devs were estimated. Set to NULL to leave it empty.
#' @param CAA_abb Abbreviation for identifying the observed age comp matrix for the catch-at-age (CAA) slot. Names of age compositions in the BAM rdat comp.mats list are expected to follow the naming convention "acomp.abb.ob". Examples from SEDAR 53 Red Grouper: "CVT", "HB", "cH". Set to "all" or "none" to use all or none of the age comps, respectively.
#' @param CAL_abb Abbreviation for identifying the observed index for the catch-at-length (CAL) slot. Analogous to CAA_abb. Names of length compositions in the BAM rdat comp.mats list are expected to follow the naming convention "lcomp.abb.ob". Set to "all" or "none" to use all or none of the length comps, respectively.
#' @param Ind_abb Abbreviation for identifying the observed index of abundance for the Ind slot. Names of indices in the BAM rdat t.series matrix are expected to follow the naming convention "U.abb.ob". Examples from SEDAR 53 Red Grouper: "CVT", "HB", "cH". If multiple (valid) abb values are provided, the corresponding indices will be averaged (geomean) and restandardized to a mean of 1. Abbreviations that don't match any indices will be ignored. Set to "all" or "none" to use all or none of the indices, respectively.
#' @param Mat_age1_max Limit maximum value of proportion mature of first age class (usually age-0 or age-1). Models sometimes fail when maturity of first age class is too high (e.g. >0.5)
#' @param length_sc  Scalar (multiplier) to convert length units including wla parameter. MSEtool examples seem to use cm whereas BAM uses mm.
#' @param wla_sc Scalar (multiplier) to convert wla parameter to appropriate weight units. In null, the function will try to figure out if the weight unit of wla was g, kg, or mt based on the range of the exponent. The wla parameter will also be scaled by 1/length_sc^wlb.
#' @param catch_sc  Scalar (multiplier) for catch. BAM catch is usually in thousand pounds (klb). The default 1 maintains that scale.
#' @param CV_vbK see \code{\link[MSEtool]{Data-class}}
#' @param CV_vbLinf see \code{\link[MSEtool]{Data-class}}
#' @param CV_vbt0 see \code{\link[MSEtool]{Data-class}}
#' @param CV_Cat see \code{\link[MSEtool]{Data-class}}
#' @param Units Units of catch. see \code{\link[MSEtool]{Data-class}}
#' @keywords bam stock assessment fisheries DLMtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Write example here
#' }

rdat_to_Data <- function(
  rdat, Data = new('Data'), herm = NULL, nsim=1,
  genus_species=NULL, Region="Southeast US Atlantic",
  Fref_name = "Fmsy",
  Rec="bam_recruits",
  CAA_abb="all",
  CAL_abb="all",
  Ind_abb="all",
  CV_vbK=0.001, CV_vbLinf=0.001, CV_vbt0=0.001,
  CV_Cat=NULL,
  Units="thousand pounds",
  Mat_age1_max = 0.49,
  length_sc=0.1,
  wla_sc=NULL,
  catch_sc=1
){

  rdat <- standardize_rdat(rdat)

  info <- rdat$info
  parms <- rdat$parms
  parm.cons <- rdat$parm.cons
  parm.tvec <- rdat$parm.tvec
  a.series <- rdat$a.series
  t.series <- rdat$t.series
  comp.mats <- rdat$comp.mats
  styr <- parms$styr
  endyr <- parms$endyr

  Name <- gsub(" ","",str_to_title(info$species))

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

  t.series <- t.series[paste(styr:endyr),]

Common_Name <- str_replace_all(Name,"(?<=[a-z])(?=[A-Z])"," ")
if(is.null(genus_species)){genus_species <- bamStockMisc[Name,"Species"]}
if(is.null(herm)){herm <- bamStockMisc[Name,"herm"]}

B.to.klb <- local({
  if(info$units.biomass%in%c("mt","metric tons")){
    out <- 2.204624}else{
      warning("units.biomass not equal to metric tons")
    }
  return(out)
})

catch.raw <- t.series$total.L.klb*catch_sc
cc.yrCat <- complete.cases(t.series$year,catch.raw) # Complete cases for catch data series
year <- t.series$year[cc.yrCat]
nyear <- length(year)
catch <- catch.raw[cc.yrCat]
recruits <- setNames(t.series$recruits,t.series$year)

Linf <- parm.cons$Linf[1]
K <- parm.cons$K[1]
t0 <- parm.cons$t0[1]

LenCV <- parm.cons$len_cv_val[1]

M.constant <- ifelse(!is.null(parms$M.constant),
                     parms$M.constant,
                     tail(a.series$M,1))

if(Fref_name=="Fmsy"){
  Cref <- parms$msy.klb*catch_sc
  Bref <- parms$Bmsy
  Fref <- parms$Fmsy
}
if(Fref_name=="F30"){
  Cref <- parms$L.F30.klb*catch_sc
  Bref <- parms$B.F30
  Fref <- parms$F30
}

Bcurrent <- t.series[paste(endyr),"B"]
B0 <- parms$B0

SSBcurrent <- t.series[paste(endyr),"SSB"]
SSB0 <- parms$SSB0

#An estimate of absolute current vulnerable abundance, converted to klb (needs to be in same units as catch)
Abun <- sum(rdat$B.age[paste(endyr),]*B.to.klb*rdat$sel.age$sel.v.wgted.tot)*catch_sc

FMSY_M <- Fref/M.constant
BMSY_B0 <- Bref/B0
Dep <- SSBcurrent/SSB0
LHYear <- endyr

# Recruitment for years where recruitment deviations were estimated
if(!is.null(Rec)){
  if(Rec=="bam_recruits"){
    Rec <- local({
      year_nodev <- parm.tvec$year[is.na(parm.tvec$log.rec.dev)]
      recruits[paste(year_nodev)] <- NA
      matrix(data=recruits,nrow=nsim,ncol=length(recruits),dimnames=list("sim"=1:nsim,"year"=year))
    })
  }
}else{
  Rec <- Data@Rec
}

# Catch (Cat): Total annual catches (NOTE: DLMtool wants Cat to be a matrix)
Cat <- matrix(data=catch,nrow=nsim,ncol=length(catch),dimnames=list("sim"=1:nsim,"year"=year))

# Catch CV
if(is.null(CV_Cat)){
CV_Cat <- matrix(0.05,nrow=nsim,ncol=nyear)
}

# Abundance Index (Ind): Relative abundance index (NOTE: DLMtool wants Ind to be a matrix)
if(!Ind_abb[1]=="none"){
IndCalc <- local({
  D <- t.series[cc.yrCat,]
  x <- names(D)
  Ind_names_all <- x[grepl(pattern="U.",x=x)&grepl(pattern=".ob",x=x)]
  if(Ind_abb[1]=="all"){
  Ind_names <- Ind_names_all
  }else{
    Ind_names <- paste("U",Ind_abb,"ob",sep=".")
  }
  Ind_names <- Ind_names[Ind_names%in%names(D)] # Identify valid names
  if(length(Ind_names)==0){
    warning(paste("Ind_abb does not match any index names in the rdat t.series. Ind will be the geometric mean of all available indices:",paste(Ind_names_all,collapse=", ")))
    Ind_names <- Ind_names_all
  }
  CV_Ind_names <- paste0("cv.U.",gsub("U.|.ob","",Ind_names))
  CV_Ind_names <- CV_Ind_names[CV_Ind_names%in%names(D)] # Identify valid names

  D[D==-99999] <- NA
  D_Ind <- D[,Ind_names,drop=FALSE]
  D_CV_Ind <- D[,CV_Ind_names,drop=FALSE]

  Ind <- apply(D_Ind,1, geomean) # Calculate geometric mean of all indices
  Ind <- Ind/mean(Ind, na.rm=TRUE) # Restandardize to mean of 1
  Ind <- matrix(data=Ind,nrow=1,ncol=length(Ind))
  Ind[is.nan(Ind)] <- NA

  CV_Ind <- apply(D_CV_Ind,1, geomean) # Calculate geometric mean of all index CVs
  CV_Ind <- matrix(data=CV_Ind,nrow=1,ncol=length(CV_Ind))
  CV_Ind[is.nan(CV_Ind)] <- NA

  return(list("Ind"=Ind,"CV_Ind"=CV_Ind))
})

slot(Data,"Ind") <- IndCalc$Ind
slot(Data,"CV_Ind") <- IndCalc$CV_Ind
}

Dt <- t.series[paste(tail(year,1)),"SSB"]/t.series[paste(year[1]),"SSB"]

# Compute proportion mature at age
pmat <- pmatage(a.series=a.series,Mat_age1_max=Mat_age1_max,herm=herm,age=age)$pmat

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
L95 <- mat_at_len$L50+mat_at_len$L50_95

# Estimate length at 5% and full (100%) vulnerability (total selectivity) and retention (landings selectivity)
vuln_out <- vulnerability(Vdata = rdat$sel.age$sel.v.wgted.tot, retdata = rdat$sel.age$sel.v.wgted.L, Linf = parm.cons$Linf[1], K = parm.cons$K[1], t0 = parm.cons$t0[1], length_sc=length_sc)

L5 <- vuln_out$L5
LFS <- vuln_out$LFS

## Age comps
if(!CAA_abb[1]=="none"){
is_acomp <- grepl(x=names(comp.mats),pattern="acomp.",fixed=TRUE)&grepl(x=names(comp.mats),pattern=".ob",fixed=TRUE)
if(any(is_acomp)){
  acomp_names_all <- names(comp.mats)[is_acomp]
if(CAA_abb[1]=="all"){
  acomp_names <- acomp_names_all
}else{
  acomp_names <- paste("acomp",CAA_abb,"ob",sep=".")
}
acomp_names <- acomp_names[acomp_names%in%names(comp.mats)] # Identify valid names
if(length(acomp_names)==0){
  warning(paste("CAA_abb does not match any names in the rdat comp.mats. CAA will be computed as a cellwise average of all available age compositions:",paste(acomp_names_all,collapse=", ")))
  acomp_names <- acomp_names_all
}
acomp_mats <- comp.mats[acomp_names] # Comps by year

names(acomp_mats) <- gsub(x=names(acomp_mats),pattern=".ob",replacement = "",fixed=TRUE)
acomp_data_n <- t.series[,grepl(x=names(t.series),pattern="acomp.",fixed=TRUE)&grepl(x=names(t.series),pattern=".n",fixed=TRUE)&!grepl(x=names(t.series),pattern=".neff",fixed=TRUE)]  # N by year
acomp_data_n[acomp_data_n==-99999] <- 0
acomp_mats_nfish <- comp_complete(acomp_mats,acomp_data_n,output_type = "nfish", val_rownames=year, val_colnames = age)

acomp_nfish <- comp_combine(acomp_mats_nfish,scale_rows = FALSE)
# Convert CAA to array to appease DLMtool
CAA <- array(data=acomp_nfish, dim=c(nsim,nrow(acomp_nfish),ncol(acomp_nfish)),
             dimnames = list("sim"=nsim,"year"=rownames(acomp_nfish),"age"=colnames(acomp_nfish)))
slot(Data,"CAA") <- CAA
}
}


## Length comps
if(!CAL_abb=="none"){
is_lcomp <- grepl(x=names(comp.mats),pattern="lcomp.",fixed=TRUE)&grepl(x=names(comp.mats),pattern=".ob",fixed=TRUE)
if(any(is_lcomp)){
  lcomp_names_all <- names(comp.mats)[is_lcomp]
  if(CAL_abb=="all"){
    lcomp_names <- lcomp_names_all
  }else{
    lcomp_names <- paste("lcomp",CAL_abb,"ob",sep=".")
  }
  lcomp_names <- lcomp_names[lcomp_names%in%names(comp.mats)] # Identify valid names
  if(length(lcomp_names)==0){
    warning(paste("CAL_abb does not match any names in the rdat comp.mats. CAL will be computed as a cellwise average of all available length compositions:",paste(lcomp_names_all,collapse=", ")))
              lcomp_names <- lcomp_names_all
  }
  lcomp_mats <- comp.mats[lcomp_names] # Comps by year

names(lcomp_mats) <- gsub(x=names(lcomp_mats),pattern=".ob",replacement = "",fixed=TRUE)
lcomp_data_n <- t.series[,grepl(x=names(t.series),pattern="lcomp.",fixed=TRUE)&grepl(x=names(t.series),pattern=".n",fixed=TRUE)]  # N by year
lcomp_data_n[lcomp_data_n==-99999] <- 0
lcomp_mats_nfish <- comp_complete(lcomp_mats,lcomp_data_n,output_type = "nfish", val_rownames=year)

lcomp_nfish <- comp_combine(lcomp_mats_nfish,scale_rows = FALSE)
colnames(lcomp_nfish) <- paste(as.numeric(colnames(lcomp_nfish))*length_sc)

# Mean
ML_LFS_out <- ML_LFS(lcomp_nfish,minL=LFS)
ML <- setNames(ML_LFS_out$mlen,ML_LFS_out$year)
#slot(Data,"ML") <-  matrix(ML,nrow=nsim,ncol=nyear,byrow=TRUE,dimnames=list("sim"=1:nsim,"year"=year))

# Convert CAL to array to appease DLMtool
CAL <- array(data=lcomp_nfish, dim=c(nsim,nrow(lcomp_nfish),ncol(lcomp_nfish)),
             dimnames = list("sim"=nsim,"year"=rownames(lcomp_nfish),"len"=colnames(lcomp_nfish)))

CAL_bins <- as.numeric(dimnames(CAL)$len)

slot(Data,"CAL") <- CAL
slot(Data,"CAL_bins") <- CAL_bins
}
}

slot(Data,"Name") <- Name
slot(Data,"Common_Name") <- Common_Name
slot(Data,"Species") <- genus_species
slot(Data,"Region") <- Region
slot(Data,"Year") <- year
slot(Data,"Cat") <- Cat
slot(Data,"CV_Cat") <- CV_Cat
slot(Data,"Rec") <- Rec
slot(Data,"t") <- nyear
slot(Data,"AvC") <- mean(Cat)
slot(Data,"Dt") <- Dt
slot(Data,"Mort") <- M.constant
slot(Data,"FMSY_M") <- FMSY_M
slot(Data,"BMSY_B0") <- BMSY_B0
slot(Data,"L50") <- L50
slot(Data,"L95") <- L95
slot(Data,"LFC") <- L5
slot(Data,"LFS") <- LFS

# von Bertalanffy K parameter (vbK): growth coefficient
slot(Data,"vbK") <- K
slot(Data,"CV_vbK") <- CV_vbK

# von Bertalanffy Linf parameter (vbLinf): Maximum length
slot(Data,"vbLinf") <- Linf*length_sc
slot(Data,"CV_vbLinf") <- CV_vbLinf

# von Bertalanffy t0 parameter (vbt0): Theoretical age at length zero
slot(Data,"vbt0") <- t0
slot(Data,"CV_vbt0") <- CV_vbt0

# Coefficient of variation of length-at-age (assumed constant for all age classes)
slot(Data,"LenCV") <- LenCV

# Length-weight parameter a (wla)
# Identify weight unit used in weight~length equation and scale parameter to compute weight in kg
wla <- parms$wgt.a
wlb <- parms$wgt.b

if(is.null(wla_sc)){
wla_sc <- local({
  if(wla>=1E-6&wla<=1E-4){
    wla_sc <- 0.001
    message(paste0("For ",Name," weight~length appears to be in grams. Scaling wla parameter by ",wla_sc))
  }
  if(wla>=1E-9&wla<=1E-7){wla_sc <- 1}      # Weight appears to already be in kilograms
  if(wla>=1E-13&wla<=1E-11){
    wla_sc <- 1000
    message(paste0("For ",Name," weight~length appears to be in metric tonnes. Scaling wla parameter by ",wla_sc))
  }
  wla_sc
})
}

# Length-weight parameter b (wlb)
if(wlb<=2|wlb>=4){
  message(paste0("For ",Name," the wlb parameter is outside of the expected range (2-4)"))
}

slot(Data,"wla") <- (wla*wla_sc)/length_sc^wlb # Adjust a parameter for length units
slot(Data,"wlb") <- wlb


slot(Data,"steep") <- parm.cons$steep[1]
slot(Data,"sigmaR") <- parm.cons$rec_sigma[1]

slot(Data,"MaxAge") <- max(age)

slot(Data,"Dep") <- Dep
slot(Data,"Abun") <- Abun
slot(Data,"SpAbun") <- SSBcurrent

slot(Data,"LHYear") <- LHYear

slot(Data,"Cref") <- Cref

slot(Data,"Units") <- Units

return(Data)
}
