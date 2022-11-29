#' Run MCBE uncertainty analysis
#'
#' This function reruns the base model, from user supplied input files or objects, creating all standard BAM ADMB files
#' most notably the executable (\code{exe}) file used in the MCBE process, stored in \code{dir_bam_base}.
#' Simulated input data sets are generated for specified fixed parameters and data sets. New BAM data input (\code{dat}) files
#' are created and the BAM base executable is rerun with each \code{dat} file, using parallel processing. After each simulation
#' is run, \code{run_BAM} checks if the objective function value (i.e. total likelihood) for a numeric value; if the value is non numeric
#' (e.g. nan), then the run is considered to have failed. The \code{dat} and results (\code{rdat}) files are for successful
#' runs are moved to \code{dir_bam_sim}. If any runs fail (which is not common) a folder is created (named \code{dir_bam_sim}
#' with suffix '_fail') and failed runs are moved there. \code{run_BAM}
#' and results file (\code{rdat}) are copied to \code{dir_bam_sim}.
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @param fileName Name given to BAM files, not including file extensions.
#' @param dir_bam_sim Name of directory to write MCBE files to, relative to the working directory.
#' @param dir_bam_base Name of directory to write bam base model files to, relative to the working directory.
#' @param bam Output of \code{bam2r}.
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param nsim number of simulations to run
#' @param data_sim list for supplying optional data not available in input data or rdat (e.g. cvs for simulating time series of landings, discards, and cpue)
#' @param par_default list for supplying default values for particular parameters to use if values cannot be found by in the usual locations (e.g. if a time series of landings does not have a corresponding time series of cvs, the default cv_L will be used)
#' @param standardize Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM
#' @param sc Scalar (multiplier) to compute upper and lower bounds of random uniform distribution from mean value
#' @param M_scLim Scalar for M (constant M) limits. Numeric vector of length 2. If M_constant is not available in the base model output, M_scLim values will be used to scale M at age.
#' @param steep_scLim Scalar for steep limits. Numeric vector of length 2
#' @param Linf_scLim Scalar for Linf (Von Bertalanffy growth function) limits. Numeric vector of length 2
#' @param K_scLim Scalar for K (Von Bertalanffy growth function) limits. Numeric vector of length 2
#' @param t0_scLim Scalar for t0 (Von Bertalanffy growth function) limits. Numeric vector of length 2
#' @param subset_rdat  Subset objects in sim rdat files. Value should be a list of the object names and either a numeric value indicating how many evenly spaced rows to include in the subset of a matrix, or NULL to set the object to NULL
#' @param Dmort_scLim Scalar for Dmort limits. Numeric vector of length 2
#' @param coresUse number of cores to use for parallel processing
#' @param ndigits number of digits to round simulated values to
#' @param unlink_dir_bam_base Should dir_bam_base be deleted after this function is run?
#' @param run_bam_base If FALSE, the function will look for an executable named fileName.exe in dir_bam_base and use it as the base model.
# If TRUE and overwrite_bam_base=TRUE, the function will call run_bam.
#' @param overwrite_bam_base If FALSE, the files in dir_bam_base will not be overwritten if run_bam_base=TRUE
#' @param admb_switch_base Character string pasted to fileName to build \code{run_command} for the base model when running BAM with \code{shell(run_command)}
#' (i.e. \code{run_command <- paste(fileName, admb_switch)})
#' @param run_sim If FALSE, the simulated data will be generated but won't be used in new BAM runs
#' @param admb_switch_sim ADMB code snippet used in shell script when running bam
#' @param prompt_me Turn on/off prompts that ask for user input before deleting files.
#' @param subset_rdat list of rdat objects to subset to decrease rdat file size
#' @param random_seed random seed value. If NULL, random seed is not set within the function.
#' @return Invisibly returns a data frame, MCBE_out, containing basic results of sim runs, including total likelihood and maximum gradient values. This data frame is also written to a csv file in \code{dir_bam_sim}.
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' Run MCBE, writing files to dir_bam_sim
#' run_MCBE("AtlanticMenhaden", dir_bam_base="AtMe_base", dir_bam_sim="AtMe_sim")
#' run_MCBE("BlackSeaBass", dir_bam_base="BlSB_base", dir_bam_sim="BlSB_sim")
#' run_MCBE("BluelineTilefish", dir_bam_base="BlTi_base", dir_bam_sim="BlTi_sim")
#' run_MCBE("Cobia", dir_bam_base="Cobi_base", dir_bam_sim="Cobi_sim")
#' run_MCBE("GagGrouper", dir_bam_base="GaGr_base", dir_bam_sim="GaGr_sim")
#' run_MCBE("GrayTriggerfish", dir_bam_base="GrTr_base", dir_bam_sim="GrTr_sim")
#' run_MCBE("GreaterAmberjack", dir_bam_base="GrAm_base", dir_bam_sim="GrAm_sim")
#' run_MCBE("RedGrouper", dir_bam_base="ReGr_base", dir_bam_sim="ReGr_sim")
#' run_MCBE("RedPorgy", dir_bam_base="RePo_base", dir_bam_sim="RePo_sim")
#' run_MCBE("RedSnapper", dir_bam_base="ReSn_base", dir_bam_sim="ReSn_sim")
#' run_MCBE("ScampGrouper", dir_bam_base="ScGr_base", dir_bam_sim="ScGr_sim")
#' run_MCBE("SnowyGrouper", dir_bam_base="SnGr_base", dir_bam_sim="SnGr_sim")
#' run_MCBE("Tilefish", dir_bam_base="Tile_base", dir_bam_sim="Tile_sim")
#' run_MCBE("VermilionSnapper", dir_bam_base="VeSn_base", dir_bam_sim="VeSn_sim")
#'
#' MCBE_ReSn <- run_MCBE("RedSnapper",steep_scLim = c(1,1))
#' }


run_MCBE <- function(CommonName = NULL,
                     fileName     = "bam",
                     dir_bam_sim  = "sim",
                     dir_bam_base = "base",
                     bam=NULL,
                     dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                     dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                     data_sim =  list(cv_U=NULL,
                                      cv_L=NULL,
                                      cv_D=NULL),
                     par_default = list(cv_U=0.2,
                                        cv_L=0.2,
                                        cv_D=0.2),
                     standardize=TRUE,
                     nsim=10,
                     sc = 0.1,  scLim = sc*c(-1,1)+1,
                     M_scLim = 0.1*c(-1,1)+1,
                     steep_scLim = scLim,
                     Linf_scLim = scLim, K_scLim = scLim, t0_scLim = scLim,
                     Dmort_scLim = scLim,
                     # parallel=TRUE, # Right now it has to be in parallel
                     coresUse=NULL,
                     ndigits=4, # number of digits to round simulated values to
                     unlink_dir_bam_base=FALSE,
                     run_bam_base=TRUE, # If FALSE, the function will look for an executable named fileName.exe in dir_bam_base and use it as the base model.
                     # If TRUE and overwrite_bam_base=TRUE, the function will call run_bam.
                     overwrite_bam_base=TRUE, # If FALSE, the files in dir_bam_base will not be overwritten if run_bam_base=TRUE
                     admb_switch_base = '-nox',
                     run_sim=TRUE, # If FALSE, the simulated data will be generated but won't be used in new BAM runs
                     admb_switch_sim = '-est -nox -ind', # -ind changes the name of the data input file each bootstrap iteration
                     prompt_me=FALSE, # Turn on/off prompts that ask for user input before deleting files.
                     subset_rdat=list("eq.series"=101,"pr.series"=101),
                     random_seed=12345


                     # admb2r_obj = admb2r.cpp,
                     # cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                     #           "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                     #           "variance","*.dep","*.hes","*.tmp"))
){
#######################
  library(doParallel)
  library(foreach)
  library(msm)

  dir_bam_sim_fail <- paste0(dir_bam_sim,"_fail")

  if(is.numeric(random_seed)){
    set.seed(random_seed)
  }

# parallel setup
if(is.null(coresUse)){
coresAvail <- detectCores()
coresUse <- coresAvail-1
}
coresUse  <- min(c(coresUse,coresAvail))
cl <- makeCluster(coresUse)
registerDoParallel(cl)

nm_sim <- sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim)

## Get base model stuff

# Identify working directory
  wd <- getwd()
  message(paste("working directory:",wd))

# Define bam base model objects
  if(!is.null(CommonName)){
    dat <- get(paste0("dat_",CommonName))
    tpl <- get(paste0("tpl_",CommonName))
    cxx <- get(paste0("cxx_",CommonName))
  }

  if(!is.null(bam)){
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
  }

  if(!is.null(dat_obj)&!is.null(tpl_obj)&!is.null(cxx_obj)){
    dat <- dat_obj
    tpl <- tpl_obj
    cxx <- cxx_obj
  }

  # Read in dat, tpl, and cxx files
  if(!is.null(dat_file)&!is.null(tpl_file)&!is.null(cxx_file)){
    dat <- readLines(con=dat_file)
    tpl <- readLines(con=tpl_file)
    cxx <- readLines(con=cxx_file)
  }

  if(standardize){
    message("Running bamExtras::standardize_bam()")
    bam <- standardize_bam(dat_obj=dat, tpl_obj=tpl,cxx_obj=cxx)
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
  }else{
      message("Running bamExtras::bam2r()")
      bam <- bam2r(dat_obj=dat, tpl_obj=tpl,cxx_obj=cxx)
  }
  init <- bam$init

  # Run bam base
  if(run_bam_base){
    if(dir.exists(dir_bam_base)){
      message(paste("The folder",paste0("'",dir_bam_base,"'"),"already exists."))
      if(overwrite_bam_base){
        message(paste("Since overwrite_bam_base = TRUE, files in ",dir_bam_base,"will be overwritten when run_bam is called"))
        spp <- run_bam(bam=bam, dir_bam = dir_bam_base, unlink_dir_bam=unlink_dir_bam_base,admb_switch=admb_switch_base)
      }else{
        message(paste("Since overwrite_bam_base = FALSE, run_bam will not be called to rerun the base run"))
        spp <- dget(file.path(dir_bam_base,paste0(fileName,".rdat")))
      }
    }else{
      dir.create(dir_bam_base)
      spp <- run_bam(bam=bam, dir_bam = dir_bam_base, unlink_dir_bam=unlink_dir_bam_base,admb_switch=admb_switch_base)

    }
  }else{
    spp <- dget(file.path(dir_bam_base,paste0(fileName,".rdat")))
  }

## Identify objects from bam base output
comp.mats <- spp$comp.mats
t.series <- spp$t.series

##############################
## Conduct Monte Carlo draws and bootstrap data

  # Check for values
    M_constant_is <- "set_M_constant"%in%names(init)
    Dmort_is <- length(grep("^set_Dmort",names(init),value=TRUE))>0


  # get base parameter values
  agebins <- as.numeric(init$agebins)

  Linf <- as.numeric(init$set_Linf[1])
  K <- as.numeric(init$set_K[1])
  t0 <- as.numeric(init$set_t0[1])

  if(M_constant_is){
    M <- as.numeric(init$set_M_constant[1])
  }else{
    M <- NA
    message("M_constant NOT found in names(init)")
  }
  Ma <- as.numeric(init$set_M)
  steep <- as.numeric(init$set_steep[1])
  rec_sigma <- as.numeric(init$set_rec_sigma[1])


  if(Dmort_is){
  Dmort <- local({
    a <- init[grep("^set_Dmort",names(init),value=TRUE)] # list
    b <- lapply(a,as.numeric)
    array(unlist(b),dim=c(1,length(b)),dimnames = list(NULL,names(b)))
    })
  }else{
    message("Discard mortality rates NOT found in the base model (i.e. no names(init) beginning with set_Dmort).")
  }

  # compute limits of parameter values
  Linf_Lim <- Linf*Linf_scLim
  K_Lim    <- K*K_scLim
  t0_Lim   <- t0*t0_scLim
  M_Lim    <- M*M_scLim
  steep_Lim <- steep*steep_scLim

  if(Dmort_is){Dmort_Lim <- Dmort[c(1,1),,drop=FALSE]*Dmort_scLim}


  ## Conduct Monte Carlo draws
  # Vectors (sim)
  sim_Linf <- runif(nsim,min(Linf_Lim),max(Linf_Lim))
  sim_K <- runif(nsim,min(K_Lim),max(K_Lim))
  sim_t0 <- runif(nsim,min(t0_Lim),max(t0_Lim))
  if(M_constant_is){
    sim_M <- runif(nsim,min(M_Lim),max(M_Lim))
    sim_Masc <- sim_M/M # Multiply by M-at-age to rescale appropriately with sim M
  }else{
    sim_M <- rep(NA,nsim)
    sim_Masc <- runif(nsim,M_scLim[1],M_scLim[2])
  }
  sim_steep <- runif(nsim,min(steep_Lim),max(steep_Lim))
  sim_rec_sigma <- round(rtnorm(n=nsim,mean=0.6,sd=0.15,lower=0.3,upper=1.0),ndigits) # values from meta-analysis. HARDCODED VALUES!!!  MAP TO FUNCTION ARGUMENTS !!!

  if(Dmort_is){
    sim_Dmort <- apply(Dmort_Lim,2,function(x){
      runif(nsim,min(x),max(x))
      })
    }

  # Matrices (sim,age)
  sim_Ma <- round(t(matrix(rep(Ma,nsim),ncol=nsim,dimnames=list(age=agebins,sim=nm_sim)))*sim_Masc,ndigits)


  sim_vectors <- local({
    a <- data.frame(Linf=sim_Linf,
                         K=sim_K,
                         t0=sim_t0,
                         M=sim_M,
                         steep=sim_steep,
                         rec_sigma=sim_rec_sigma
                         )
    if(Dmort_is) {b <- cbind(a,sim_Dmort)
    }else{
      b <- a
    }
    round(b,ndigits)
  })

  # reproductive parameters

  # Bootstrap data


  #%%%% bootstrap setup %%%%#
  #%% Indices %%#
  obs_cpue_nm <- names(init)[grepl(pattern="obs_cpue_",names(init))]
  cv_cpue_nm <- names(init)[grepl(pattern="cv_cpue_",names(init))]
  cpue_root <- gsub("obs_cpue_","",obs_cpue_nm)

  # if(!is.null(data_sim$cv_cpue)){
  #   cv_cpue <-
  # }

  sim_cpue <- list()
  for(nm_i in obs_cpue_nm){
    cpue_root_i <- gsub("obs_cpue_","",nm_i)
    cv_cpue_nm_i <- gsub("obs_cpue_","obs_cv_cpue_",nm_i)

    yrs_cpue_nm_i <- paste0("yrs_cpue_",cpue_root_i)
    styr_cpue_nm_i <- gsub("obs","styr",nm_i)
    endyr_cpue_nm_i <- gsub("obs","endyr",nm_i)
    if(yrs_cpue_nm_i%in%names(init)){ # If the specific years of the index are given, use them. Otherwise look for range of years given and use those
      yrs_cpue_i <- init[[yrs_cpue_nm_i]]
    }else{
      styr_cpue_i  <- init[[styr_cpue_nm_i]]
      endyr_cpue_i <- init[[endyr_cpue_nm_i]]
      yrs_cpue_i <- paste(styr_cpue_i:endyr_cpue_i)
    }

    cpue_i <- as.numeric(init[[nm_i]])

    # Set default cvs
    cv_cpue_i <- if(cv_cpue_nm_i%in%names(init)){
      as.numeric(init[[cv_cpue_nm_i]])}else{
        message(paste0(cv_cpue_nm_i," not found in names(init). Using cvs of par_default$cv_U = ",par_default$cv_U, " to simulate ",nm_i))
        rep(par_default$cv_U,length(cpue_i))
      }

    # Override and replace cvs if values are provided in data_sim
    if(any(!is.na(data_sim$cv_U[,cpue_root_i]))){
      message(paste("using values from data_sim for",nm_i))
      ci <- local({
        ai <- data_sim$cv_U[,cpue_root_i,drop=FALSE]
        ai[complete.cases(ai),,drop=FALSE]
      })
      cv_cpue_i <- ci[,cpue_root_i]
      if(any(yrs_cpue_i != rownames(ci))){warning(paste("years of",paste0("data_sim$cv_U$",cpue_root_i),"do not match years of",nm_i,"in the base model"))}
    }



    cv_cpue_i <- as.numeric(init[[cv_cpue_nm_i]])

    # Override cvs if values are provided in data_sim
    if(any(!is.na(data_sim$cv_U[,cpue_root_i]))){
      message(paste("using values from data_sim for",nm_i))
      ci <- local({
        ai <- data_sim$cv_U[,cpue_root_i,drop=FALSE]
        ai[complete.cases(ai),,drop=FALSE]
      })
      cv_cpue_i <- ci[,cpue_root_i]
      yrs_cpue_i <- rownames(ci)
    }

    sim_cpue_i <- lnorm_vector_boot(cpue_i,cv_cpue_i,nsim,
                                 standardize=TRUE,digits=ndigits)
    dimnames(sim_cpue_i) <- list("year"=yrs_cpue_i,"sim"=nm_sim)
    sim_cpue[[cpue_root_i]] <- sim_cpue_i
  }

  #%% Landings %%#
  obs_L_nm <- names(init)[grepl(pattern="obs_L_",names(init))]
  sim_L <- list()
  cv_L_nm <- names(init)[grepl(pattern="cv_L_",names(init))]

  for(nm_i in obs_L_nm){
    L_root_i <- gsub("obs_L_","",nm_i)
    cv_L_nm_i <- gsub("obs_L_","obs_cv_L_",nm_i)

    yrs_L_nm_i <- paste0("yrs_L_",L_root_i)
    styr_L_nm_i <- gsub("obs","styr",nm_i)
    endyr_L_nm_i <- gsub("obs","endyr",nm_i)
    if(yrs_L_nm_i%in%names(init)){ # If the specific years are given, use them. Otherwise look for range of years given and use those.
      yrs_L_i <- init[[yrs_L_nm_i]]
    }else{
      styr_L_i  <- init[[styr_L_nm_i]]
      endyr_L_i <- init[[endyr_L_nm_i]]
      yrs_L_i <- paste(styr_L_i:endyr_L_i)
    }

    L_i <- as.numeric(init[[nm_i]])

    # Set default cvs
    cv_L_i <- if(cv_L_nm_i%in%names(init)){
      as.numeric(init[[cv_L_nm_i]])}else{
        message(paste0(cv_L_nm_i," not found in names(init). Using cvs of par_default$cv_L = ",par_default$cv_L, " to simulate ",nm_i))
        rep(par_default$cv_L,length(L_i))
      }

    # Override and replace cvs if values are provided in data_sim
    if(any(!is.na(data_sim$cv_L[,L_root_i]))){
      message(paste("using values from data_sim for",nm_i))
      ci <- local({
        ai <- data_sim$cv_L[,L_root_i,drop=FALSE]
        ai[complete.cases(ai),,drop=FALSE]
      })
      cv_L_i <- ci[,L_root_i]
      if(any(yrs_L_i != rownames(ci))){warning(paste("years of",paste0("data_sim$cv_L$",L_root_i),"do not match years of",nm_i,"in the base model"))}
    }

    sim_L_i <- lnorm_vector_boot(L_i,cv_L_i,nsim,
                                    standardize=FALSE,digits=ndigits)
    dimnames(sim_L_i) <- list("year"=yrs_L_i,"sim"=nm_sim)
    sim_L[[L_root_i]] <- sim_L_i
  }

  #%% Discards (released) %%#
  obs_released_nm <- names(init)[grepl(pattern="obs_released_",names(init))]
  sim_D <- list()

  if(length(obs_released_nm)>0){
  cv_D_nm <- names(init)[grepl(pattern="obs_cv_D_",names(init))]

  for(nm_i in obs_released_nm){
    D_root_i <- gsub("obs_released_","",nm_i)
    cv_D_nm_i <- gsub("obs_released_","obs_cv_D_",nm_i)

    yrs_D_nm_i <- paste0("yrs_D_",D_root_i)
    styr_D_nm_i <- gsub("obs_released","styr_D",nm_i)
    endyr_D_nm_i <- gsub("obs_released","endyr_D",nm_i)
    if(yrs_D_nm_i%in%names(init)){ # If the specific years are given, use them. Otherwise look for range of years given and use those
      yrs_D_i <- init[[yrs_D_nm_i]]
    }else{
      styr_D_i  <- init[[styr_D_nm_i]]
      endyr_D_i <- init[[endyr_D_nm_i]]
      yrs_D_i <- paste(styr_D_i:endyr_D_i)
    }

    D_i <- as.numeric(init[[nm_i]])

    # Set default cvs
    cv_D_i <- if(cv_D_nm_i%in%names(init)){
      as.numeric(init[[cv_D_nm_i]])}else{
        message(paste0(cv_D_nm_i," not found in names(init). Using cvs of par_default$cv_D = ",par_default$cv_D, " to simulate ",nm_i))
        rep(par_default$cv_D,length(D_i))
      }

    # Override and replace cvs if values are provided in data_sim
    if(any(!is.na(data_sim$cv_D[,D_root_i]))){
      message(paste("using values from data_sim for",nm_i))
      ci <- local({
        ai <- data_sim$cv_D[,D_root_i,drop=FALSE]
        ai[complete.cases(ai),,drop=FALSE]
      })
      cv_D_i <- ci[,D_root_i]
      if(any(yrs_D_i != rownames(ci))){warning(paste("years of",paste0("data_sim$cv_D$",D_root_i),"do not match years of",nm_i,"in the base model"))}
    }

    sim_D_i <- lnorm_vector_boot(D_i,cv_D_i,nsim,
                                 standardize=FALSE,digits=ndigits)
    dimnames(sim_D_i) <- list("year"=yrs_D_i,"sim"=nm_sim)
    sim_D[[D_root_i]] <- sim_D_i
  }
  }else{
  message("No discard time series found in the base model (i.e. no names(init) beginning with obs_released).")
  }

  #%% Age composition %%#
  obs_agec_nm <- names(init)[grepl(pattern="obs_agec",names(init))]
  sim_acomp <- list()

  if(length(obs_agec_nm)>0){
  agec_root <- gsub("^obs_agec_","",obs_agec_nm)

  for(agec_root_i in agec_root){
    acomp_ob_nm_i <- gsub("\\_","\\.",paste0("acomp.",agec_root_i,".ob"))
    x_i <- comp.mats_i <- comp.mats[[acomp_ob_nm_i]]

    agebins_acomp_i <- factor(colnames(x_i),levels=colnames(x_i))

    yrs_i <- as.numeric(rownames(x_i))

    # Get nfish from init instead of spp. If years of comps are excluded from
    # fitting by bam based on minSS then values of -99999 will appear for
    # values of nfish in the rdat. This won't affect results, but the -99999
    # values cause an error when trying to simulate comps.
    nfish_i <- setNames(as.numeric(init[[paste0("nfish_agec_",agec_root_i)]]),paste(yrs_i))
    dimnames(x_i) <- list("year"=yrs_i,"age"=agebins_acomp_i)

    x_i <- as.data.frame.matrix(x_i)
    x_i$nfish <- nfish_i


    sim_acomp_i <- foreach(i=1:nsim) %dopar% {
      out <- apply(x_i,1,function(y){
        if(sum(y[names(y)!="nfish"])>0&y["nfish"]>0){
          ages_y <- sample(x=agebins_acomp_i, size=y["nfish"], replace=TRUE, prob=y[paste(agebins_acomp_i)])
          round(table(factor(ages_y,agebins_acomp_i))/y["nfish"],4)
        }else{
          ages_y <- rep(0,length(agebins_acomp_i))
        }
      })
      t(out)
    }
    sim_acomp_i <- array(unlist(sim_acomp_i),dim=c(dim(comp.mats_i),length(sim_acomp_i)),
                         dimnames=list(year=yrs_i,age=levels(agebins_acomp_i),sim=nm_sim))

    sim_acomp[[agec_root_i]] <- sim_acomp_i
  }
  }else{
    message("No age compositions found in the base model (i.e. no names(init) beginning with obs_agec).")
  }

  #%% Length composition %%#
  obs_lenc_nm <- names(init)[grepl(pattern="obs_lenc",names(init))]
  sim_lcomp <- list()

  if(length(obs_lenc_nm)>0){
    lenc_root <- gsub("^obs_lenc_","",obs_lenc_nm)

    for(lenc_root_i in lenc_root){
      lcomp_ob_nm_i <- gsub("\\_","\\.",paste0("lcomp.",lenc_root_i,".ob"))
      x_i <- comp.mats_i <- comp.mats[[lcomp_ob_nm_i]]

      lenbins_lcomp_i <- factor(colnames(x_i),levels=colnames(x_i))

      yrs_i <- as.numeric(rownames(x_i))

      # Get nfish from init instead of spp. If years of comps are excluded from
      # fitting by bam based on minSS then values of -99999 will appear for
      # values of nfish in the rdat. This won't affect results, but the -99999
      # values cause an error when trying to simulate comps.
      nfish_i <- setNames(as.numeric(init[[paste0("nfish_lenc_",lenc_root_i)]]),paste(yrs_i))

      dimnames(x_i) <- list("year"=yrs_i,"len"=lenbins_lcomp_i)

      x_i <- as.data.frame.matrix(x_i)
      x_i$nfish <- nfish_i

      sim_lcomp_i <- foreach(i=1:nsim) %dopar% {
        out <- apply(x_i,1,function(y){
          if(sum(y[names(y)!="nfish"])>0&y["nfish"]>0){
            lens_y <- sample(x=lenbins_lcomp_i, size=y["nfish"], replace=TRUE, prob=y[paste(lenbins_lcomp_i)])
            round(table(factor(lens_y,lenbins_lcomp_i))/y["nfish"],4)
          }else{
            lens_y <- rep(0,length(lenbins_lcomp_i))
          }
        })
        t(out)
      }
      sim_lcomp_i <- array(unlist(sim_lcomp_i),dim=c(dim(comp.mats_i),length(sim_lcomp_i)),
                           dimnames=list(year=yrs_i,len=levels(lenbins_lcomp_i),sim=nm_sim))

      sim_lcomp[[lenc_root_i]] <- sim_lcomp_i
    }
  }else{
    message("No length compositions found in the base model (i.e. no names(init) beginning with obs_lenc).")
  }

############################
## Run the bam base model to build the executable and establish a reference?
  #%%%%%%%%%%%%%%%%%%%%%
  #### Run MCBE (maybe this should be a separate function?) ####
  #%%%%%%%%%%%%%%%%%%%%%
  if(run_sim){
    if(dir.exists(dir_bam_sim)){
      if(prompt_me){
        delete_files <- readline(prompt=paste0(paste("The folder",paste0("'",dir_bam_sim,"'"),"already exists. "),"Are you sure you want to delete all files in ",dir_bam_sim,"? (Enter TRUE or FALSE)"))
      }else{
        delete_files <- TRUE
      }
      if(delete_files){
        unlink(paste0(dir_bam_sim,"/*"))
        message(paste("The contents of the folder",paste0("'",dir_bam_sim,"'"),"have been deleted."))
        if(dir.exists(dir_bam_sim_fail)){
          unlink(paste0(dir_bam_sim_fail,"/*"))
          message(paste("The contents of the folder",paste0("'",dir_bam_sim_fail,"'"),"have been deleted."))
        }
        if(length(list.files(dir_bam_sim))>0){
          stop(paste("The contents of the folder",paste0("'",dir_bam_sim,"'"),"were NOT actually deleted."))
        }

      }else{
        message(paste("If you do not want to delete the contents of ",dir_bam_sim,", please specify a new value of dir_bam_sim and rerun."))
      }
    }else{
      dir.create(dir_bam_sim)
      message(paste("Created the folder",paste0("'",dir_bam_sim,"'.")))
    }

  message(paste("Running MCBE for", nsim,"sims in parallel on",coresUse,"cores at",Sys.time()))

  MCBE_out <- foreach(i=1:nsim,
          #, .export = ls()[grepl("xboot.obs",ls())] # Pass additional objects to foreach
          .packages=c("bamExtras")
  ) %dopar% {
    nm_sim_i <- nm_sim[i]

    init_i <- init # Copy init from base for each bootstrap run, to initialize

    #%%%%  Modify init_i %%%%#

    #%% Create Monte Carlo/bootstrap data set %%#

    ##### Monte Carlo #####
    # M
    if(M_constant_is){
      init_i$set_M_constant[c(1,5)] <- paste(sim_vectors$M[i])
    }
    init_i$set_M <- paste(sim_Ma[i,])

    # Dmort
    if(Dmort_is){
    for(j in colnames(Dmort)){
      set_Dmort_ij <- sim_vectors[i,j]
      init_i[[j]] <- paste(set_Dmort_ij)
    }
    }

    # steepness
    init_i$set_steep[c(1,5)] <- paste(sim_vectors$steep[i])

    # rec_sigma
    init_i$set_rec_sigma[c(1,5)] <- paste(sim_vectors$rec_sigma[i])

    ##### Bootstrap #####
    spf <- paste0("%01.",ndigits,"f")
    # Add indices (cpue)
    for(nm_j in obs_cpue_nm){
      cpue_root_j <- gsub("obs_cpue_","",nm_j)
      init_i[[nm_j]] <- sprintf(spf,sim_cpue[[cpue_root_j]][,i])
    }

    # Add landings (L)
    for(nm_j in obs_L_nm){
      L_root_j <- gsub("obs_L_","",nm_j)
      init_i[[nm_j]] <- sprintf(spf,sim_L[[L_root_j]][,i])
    }

    # Add discards (D)
    if(length(obs_released_nm)>0){
      for(nm_j in obs_released_nm){
        released_root_j <- gsub("obs_released_","",nm_j)
        init_i[[nm_j]] <- sprintf(spf,sim_D[[released_root_j]][,i])
      }
    }

    # Add age comps
    if(length(obs_agec_nm)>0){
      for(agec_root_j in agec_root){
        # acomp_ob_nm_j <- paste0("acomp.",agec_root_j,".ob") # rdat naming convention
        obs_agec_nm_j <- gsub("\\.","\\_",paste0("obs_agec_",agec_root_j))    # tpl naming convention
        agec_ij <- sim_acomp[[agec_root_j]][,,i,drop=FALSE]
        dim_agec_ij <- dim(agec_ij)
        agec_ij <- matrix(apply(agec_ij,2,function(x){sprintf(spf,x)}),nrow=dim_agec_ij[1])
        dimnames(agec_ij) <- NULL
        init_i[[obs_agec_nm_j]] <- agec_ij
      }
    }

    if(length(obs_lenc_nm)>0){
      for(lenc_root_j in lenc_root){
        # lcomp_ob_nm_j <- paste0("lcomp.",lenc_root_j,".ob") # rdat naming convention
        obs_lenc_nm_j <- gsub("\\.","\\_",paste0("obs_lenc_",lenc_root_j))    # tpl naming convention
        lenc_ij <- sim_lcomp[[lenc_root_j]][,,i,drop=FALSE]
        dim_lenc_ij <- dim(lenc_ij)
        lenc_ij <- matrix(apply(lenc_ij,2,function(x){sprintf(spf,x)}),nrow=dim_lenc_ij[1])
        dimnames(lenc_ij) <- NULL
        init_i[[obs_lenc_nm_j]] <- lenc_ij
      }
    }

    #%% Incorporate changes to init_i back into BAM dat %%#
    bam_i <- bam2r(
      dat_obj = dat,
      tpl_obj = tpl,
      cxx_obj = cxx,
      init = init_i)

    #%%  File management stuff
    sim_dir_i <- paste0("sim_",nm_sim_i)
    dir.create(sim_dir_i)
    setwd(sim_dir_i)

    fileName_exe_base <- paste0(fileName,".exe")

    fileName_dat_i <- paste(nm_sim_i,'-',fileName,'.dat',sep="") # Name of dat file for i
    fileName_rdat_i <- paste(nm_sim_i,'-',fileName,'.rdat',sep="") # Name of rdat file for i
    fileName_exe_i <- paste(nm_sim_i,'-',fileName,'.exe',sep="") # Name of exe file for i
    fileName_par_i <- paste(nm_sim_i,'-',fileName,'.par',sep="") # Name of par file for i

    # Copy executable to directory for sim i
    file.copy(file.path("..",dir_bam_base,fileName_exe_base), fileName_exe_i, overwrite=TRUE)

    # Rewrite dat file incorporating modified data
    writeLines(text=bam_i$dat, con=fileName_dat_i)

    #######Run bam for sim_i
    shell(paste(fileName_exe_i, admb_switch_sim, fileName_dat_i, sep=" "))

    par_i <- readLines(fileName_par_i)
    lk_total_i <- gsub("^.*Objective function value = (.*)  Maximum.*$","\\1",par_i[1])
    grad_max_i <- gsub("^.*component = ","\\1",par_i[1])

    if(!is.na(as.numeric(lk_total_i))){ # If the total likelihood of the model was a numeric result
    # Subset rdat to decrease size on disk
    rdat_i <- dget(file=fileName_rdat_i)
    for(nm_k in names(subset_rdat)){
      k <- rdat_i[[nm_k]]
      if(is.null(subset_rdat[[nm_k]])){
        rdat_i[[nm_k]] <- NULL
      }else{
        rdat_i[[nm_k]] <- k[seq(1,nrow(k),length=subset_rdat[[nm_k]]),,drop=FALSE]
      }
    }
    # Save file over previous
    dput(x=rdat_i,file=fileName_rdat_i)
    file_to <- file.path("..",dir_bam_sim)
    }else{
      file_to <- file.path("..",dir_bam_sim_fail)
      if(!dir.exists(file_to)){
        dir.create(file_to)
      }
    }

    ########## Copy data files to folder
    file.copy(from=c(fileName_dat_i,fileName_rdat_i),
              to=file_to,
              overwrite = TRUE)

    #######Remove individual processing folders
    setwd(wd)
    unlink(sim_dir_i, recursive=T)

  return(setNames(c(lk_total_i,grad_max_i),c("lk_total","grad_max")))
  } #end MCBE foreach loop


  message(paste("Finished running MCBE for", nsim,"sims in parallel on",coresUse,"cores at",Sys.time()))

  MCBE_out <- apply(do.call(rbind,MCBE_out),2,as.numeric)
  MCBE_out <- cbind("sim"=nm_sim,as.data.frame(MCBE_out),
                    "Linf"=sim_Linf,"K"=sim_K, "t0"=sim_t0, "M"=sim_M, "Mage_scale"=sim_Masc,
                    "steep"=sim_steep,"rec_sigma"=sim_rec_sigma)

  write.csv(x=MCBE_out,file=file.path(dir_bam_sim,"MCBE_results.csv"),row.names = FALSE)

### Return stuff
invisible(MCBE_out)

  } # end if(run_sim)

  ## parallel shut down
  stopCluster(cl)
}
