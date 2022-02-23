#' Build MSEtool Obs (observation) object from bam rdat object
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Obs MSEtool Obs object to start with
#' @param CAA_abb Abbreviation for identifying the observed age comp matrix for the catch-at-age (CAA) slot. Names of age compositions in the BAM rdat comp.mats list are expected to follow the naming convention "acomp.abb.ob". Examples from SEDAR 53 Red Grouper: "CVT", "HB", "cH". Set to "all" or "none" to use all or none of the age comps, respectively.
#' @param CAL_abb Abbreviation for identifying the observed index for the catch-at-length (CAL) slot. Analogous to CAA_abb. Names of length compositions in the BAM rdat comp.mats list are expected to follow the naming convention "lcomp.abb.ob". Set to "all" or "none" to use all or none of the length comps, respectively.
#' @param Ind_abb Abbreviation for identifying the observed index of abundance for the Ind slot. Names of indices in the BAM rdat t.series matrix are expected to follow the naming convention "U.abb.ob". Examples from SEDAR 53 Red Grouper: "CVT", "HB", "cH". If multiple (valid) abb values are provided, the corresponding indices will be averaged (geomean) and restandardized to a mean of 1. Abbreviations that don't match any indices will be ignored. Set to "all" or "none" to use all or none of the indices, respectively.
#' @param Cobs Range of CVs for landings data. Default set to value used in BAM modelsee \code{\link[MSEtool]{Obs-class}}
#' @param catch_sc  Scalar (multiplier) for catch. BAM catch is usually in thousand pounds (klb). The default 1 maintains that scale.
#' @param CAA_ESS_sc  Scalar multiplied by CAA_nsamp (in numbers of fish) to compute effective sample size of age composition data (CAA_ESS) when empirical estimates of effective sample size are not available. Default of 0.3 is based on MSEtool generic Obs objects.
#' @param CAL_ESS_sc  Scalar multiplied by CAL_nsamp (in numbers of fish) to compute effective sample size of length composition data (CAL_ESS) when empirical estimates of effective sample size are not available. Default of 0.3 is based on MSEtool generic Obs objects.
#' @keywords bam stock assessment fisheries MSEtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Build MSEtool Stock, Fleet, and Obs (Stock-class, Fleet-class, and Obs-class objects)
#' Stock_RedPorgy <- rdat_to_Stock(rdat_RedPorgy)
#' Fleet_RedPorgy <- rdat_to_Fleet(rdat_RedPorgy)
#' Obs_RedPorgy <- rdat_to_Obs(rdat_RedPorgy)
#' 
#' # Build MSEtool operating model (OM-class object)
#' OM_RedPorgy <- new("OM", Stock_RedPorgy, Fleet_RedPorgy, Obs_RedPorgy, Perfect_Imp)
#' # Run and plot simple management strategy evaluation (MSE)
#' mse_out <- runMSE(OM_RedPorgy)
#' NOAA_plot(mse_out)
#' }

rdat_to_Obs <- function(
  rdat,
  Obs = Perfect_Info,
  Ind_abb="all",
  CAA_abb="all",
  CAL_abb="all",
  Cobs = c(0.05, 0.05),
  catch_sc=1,
  CAA_ESS_sc=0.3,
  CAL_ESS_sc=0.3
){
  require(stringr)
  
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

  catch.raw <- t.series$total.L.klb*catch_sc
  cc.yrCat <- complete.cases(t.series$year,catch.raw) # Complete cases for catch data series
  year <- t.series$year[cc.yrCat]
  nyear <- length(year)  
  
  ## Abundance Index (Ind): Relative abundance index (NOTE: DLMtool wants Ind to be a matrix)
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
      
      return(list("Ind"=Ind,"CV_Ind"=CV_Ind,"Ind_names"=Ind_names))
    })
  }
  
  
  ## Age comps
  if(!CAA_abb[1]=="none"){
    is_acomp <- grepl(x=names(comp.mats),pattern="acomp.",fixed=TRUE)&grepl(x=names(comp.mats),pattern=".ob",fixed=TRUE)
    if(any(is_acomp)){
      acomp_names_all <- names(comp.mats)[is_acomp]
      acomp_n_names_all <- names(t.series)[grepl("^acomp.*.n$",names(t.series))]
      acomp_nfish_names_all <- names(t.series)[grepl("^acomp.*.nfish$",names(t.series))]
      acomp_neff_names_all <- names(t.series)[grepl("^acomp.*.neff$",names(t.series))]
      
      if(CAA_abb[1]=="all"){
        acomp_names <- acomp_names_all
        acomp_n_names <- acomp_n_names_all
        acomp_nfish_names <- acomp_nfish_names_all
        acomp_neff_names <- acomp_neff_names_all
      }else{
        acomp_names <- paste("acomp",CAA_abb,"ob",sep=".")
        acomp_n_names <- paste("acomp",CAA_abb,"n",sep=".")
        acomp_nfish_names <- paste("acomp",CAA_abb,"nfish",sep=".")
        acomp_neff_names <- paste("acomp",CAA_abb,"neff",sep=".")
      }
      acomp_names <- acomp_names[acomp_names%in%names(comp.mats)] # Identify valid names
      acomp_n_names <- acomp_n_names[acomp_n_names%in%names(t.series)] # Identify valid names
      acomp_nfish_names <- acomp_nfish_names[acomp_nfish_names%in%names(t.series)] # Identify valid names
      acomp_neff_names <- acomp_neff_names[acomp_neff_names%in%names(t.series)] # Identify valid names
      
      if(length(acomp_names)==0){
        warning(paste("CAA_abb does not match any names in the rdat comp.mats. CAA will be computed as a cellwise average of all available age compositions:",paste(acomp_names_all,collapse=", ")))
        acomp_names <- acomp_names_all
        acomp_n_names <- acomp_n_names_all
        acomp_nfish_names <- acomp_nfish_names_all
        acomp_neff_names <- acomp_neff_names_all
      }
      acomp_mats <- comp.mats[acomp_names] # Comps by year
      
      names(acomp_mats) <- gsub(x=names(acomp_mats),pattern=".ob",replacement = "",fixed=TRUE)
      acomp_data_n <- t.series[,acomp_n_names,drop=FALSE]          # number of trips by year
      acomp_data_n[acomp_data_n<0] <- 0
      
      acomp_data_nfish <- t.series[,acomp_nfish_names,drop=FALSE]  # number of fish by year
      acomp_data_nfish[acomp_data_nfish<0] <- 0
      
      acomp_data_neff <- t.series[,acomp_neff_names,drop=FALSE]    # effective number of trips by year
      acomp_data_neff[acomp_data_neff<0] <- 0
      
      if(!any(dim(acomp_data_neff)==0)){
      # Compute ratio of number of trips to effective number of trips
      acomp_data_n2neff <- local({
        a <- gsub(".n",".neff",names(acomp_data_n))
        ix <-which(a%in%names(acomp_data_neff))
        acomp_data_n[ix]/acomp_data_neff
      })
      
      # Approximate effective sample size for nfish
      acomp_data_nfisheff <- local({
        n_ix <- gsub(".[a-zA-Z]+$","",names(acomp_data_nfish))%in%gsub(".[a-zA-Z]+$","",names(acomp_data_n2neff))
        acomp_data_nfish[n_ix]/acomp_data_n2neff
      })
    }
      # acomp_mats_nfish <- comp_complete(acomp_mats,acomp_data_nfish,n_tag = ".nfish",output_type = "nfish", val_rownames=year, val_colnames = age)
      # 
      # acomp_nfish <- comp_combine(acomp_mats_nfish,scale_rows = FALSE)
      # # Convert CAA to array to appease DLMtool
      # CAA <- array(data=acomp_nfish, dim=c(nsim,nrow(acomp_nfish),ncol(acomp_nfish)),
      #              dimnames = list("sim"=nsim,"year"=rownames(acomp_nfish),"age"=colnames(acomp_nfish)))
    }else{
      warning("No acomp found in comp.mats.")
    }
  }
  
  
  ## Length comps
  if(!CAL_abb[1]=="none"){
    is_lcomp <- grepl(x=names(comp.mats),pattern="lcomp.",fixed=TRUE)&grepl(x=names(comp.mats),pattern=".ob",fixed=TRUE)
    if(any(is_lcomp)){
      lcomp_names_all <- names(comp.mats)[is_lcomp]
      lcomp_n_names_all <- names(t.series)[grepl("^lcomp.*.n$",names(t.series))]
      lcomp_nfish_names_all <- names(t.series)[grepl("^lcomp.*.nfish$",names(t.series))]
      lcomp_neff_names_all <- names(t.series)[grepl("^lcomp.*.neff$",names(t.series))]
      
      if(CAL_abb[1]=="all"){
        lcomp_names <- lcomp_names_all
        lcomp_n_names <- lcomp_n_names_all
        lcomp_nfish_names <- lcomp_nfish_names_all
        lcomp_neff_names <- lcomp_neff_names_all
      }else{
        lcomp_names <- paste("lcomp",CAL_abb,"ob",sep=".")
        lcomp_n_names <- paste("lcomp",CAL_abb,"n",sep=".")
        lcomp_nfish_names <- paste("lcomp",CAL_abb,"nfish",sep=".")
        lcomp_neff_names <- paste("lcomp",CAL_abb,"neff",sep=".")
      }
      
      lcomp_names <- lcomp_names[lcomp_names%in%names(comp.mats)] # Identify valid names
      lcomp_n_names <- lcomp_n_names[lcomp_n_names%in%names(t.series)] # Identify valid names
      lcomp_nfish_names <- lcomp_nfish_names[lcomp_nfish_names%in%names(t.series)] # Identify valid names
      lcomp_neff_names <- lcomp_neff_names[lcomp_neff_names%in%names(t.series)] # Identify valid names
      
      if(length(lcomp_names)==0){
        warning(paste("CAL_abb does not match any names in the rdat comp.mats. CAL will be computed as a cellwise average of all available length compositions:",paste(lcomp_names_all,collapse=", ")))
        lcomp_names <- lcomp_names_all
        lcomp_n_names <- lcomp_n_names_all
        lcomp_nfish_names <- lcomp_nfish_names_all
        lcomp_neff_names <- lcomp_neff_names_all
      }
      lcomp_mats <- comp.mats[lcomp_names] # Comps by year
      
      names(lcomp_mats) <- gsub(x=names(lcomp_mats),pattern=".ob",replacement = "",fixed=TRUE)
      lcomp_data_n <- t.series[,lcomp_n_names,drop=FALSE]          # number of trips by year
      lcomp_data_n[lcomp_data_n<0] <- 0
      
      lcomp_data_nfish <- t.series[,lcomp_nfish_names,drop=FALSE]  # number of fish by year
      lcomp_data_nfish[lcomp_data_nfish<0] <- 0
      
      lcomp_data_neff <- t.series[,lcomp_neff_names,drop=FALSE]    # effective number of trips by year
      lcomp_data_neff[lcomp_data_neff<0] <- 0
      
      if(!any(dim(lcomp_data_neff)==0)){
      # Compute ratio of number of trips to effective number of trips
      lcomp_data_n2neff <- local({
        a <- gsub(".n",".neff",names(lcomp_data_n))
        ix <-which(a%in%names(lcomp_data_neff))
        lcomp_data_n[ix]/lcomp_data_neff
      })
      
      # Approximate effective sample size for nfish
    lcomp_data_nfisheff <- local({
      n_ix <- gsub(".[a-zA-Z]+$","",names(lcomp_data_nfish))%in%gsub(".[a-zA-Z]+$","",names(lcomp_data_n2neff))
      lcomp_data_nfish[n_ix]/lcomp_data_n2neff
      })

      }
      # lcomp_nfish <- comp_combine(lcomp_mats_nfish,scale_rows = FALSE)
      # colnames(lcomp_nfish) <- paste(as.numeric(colnames(lcomp_nfish))*length_sc)
      # 
      # # Convert CAL to array to appease DLMtool
      # CAL <- array(data=lcomp_nfish, dim=c(nsim,nrow(lcomp_nfish),ncol(lcomp_nfish)),
      #              dimnames = list("sim"=nsim,"year"=rownames(lcomp_nfish),"len"=colnames(lcomp_nfish)))

    }else{
      warning("No lcomp found in comp.mats.")
    }
  }
  
  slot(Obs,"Iobs")      <- range(IndCalc$CV_Ind,na.rm=TRUE)
  
  if(exists("acomp_data_nfish")){
  slot(Obs,"CAA_nsamp") <- range(acomp_data_nfish[acomp_data_nfish>0],na.rm=TRUE)
  }
  
  if(exists("acomp_data_neff")){
  if(!any(dim(acomp_data_neff)==0)){
    slot(Obs,"CAA_ESS") <- round(range(acomp_data_nfisheff[acomp_data_nfisheff>0],na.rm=TRUE))
  }}else{
    warning(paste("Cannot estimate CAA_ESS from empirical data. Scaling CAA_nsamp by CAA_ESS_sc to compute CAA_ESS."))
    slot(Obs,"CAA_ESS") <- round(slot(Obs,"CAA_nsamp")*CAA_ESS_sc)
  }

  slot(Obs,"CAL_nsamp") <- range(lcomp_data_nfish[lcomp_data_nfish>0],na.rm=TRUE)  

  if(exists("lcomp_data_neff")){
  if(!any(dim(lcomp_data_neff)==0)){
    slot(Obs,"CAL_ESS") <- round(range(lcomp_data_nfisheff[lcomp_data_nfisheff>0],na.rm=TRUE))
  }}else{
    warning(paste("Cannot estimate CAL_ESS from empirical data. Scaling CAL_nsamp by CAL_ESS_sc to compute CAL_ESS."))
    slot(Obs,"CAL_ESS") <- round(slot(Obs,"CAL_nsamp")*CAL_ESS_sc)
  }
  
return(Obs)
}
