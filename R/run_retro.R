#' Simulate bam models for running retrospective analysis
#'
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @param fileName Name given to BAM files, not including file extensions.
#' @param dir_bam_sim Name of directory to write retrospective analysis files to, relative to the working directory.
#' @param dir_bam_base Name of directory to write bam base model files to, relative to the working directory.
#' @param bam Output of \code{bam2r}.
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param standardize Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM
#' @param nyr_remove number of years to remove in the retrospective analysis
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
#' @param admb2r_obj Character string containing admb2r C++ code, which is written with \code{base::writeLines} to \code{dir_bam}
#' @param cleanup List object written to \code{cleanup.bat} file in \code{dir_bam}.
#'
#' @return Invisibly returns a data frame, sim_out, containing basic results of sim runs, including total likelihood and maximum gradient values. This data frame is also written to a csv file in \code{dir_bam_sim}.
#' @keywords bam retrospective analysis stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' Run retrospective analysis, writing files to dir_bam_sim
#' run_retro("AtlanticMenhaden", dir_bam_base="AtMe_base", dir_bam_sim="AtMe_sim")
#' run_retro("BlackSeaBass", dir_bam_base="BlSB_base", dir_bam_sim="BlSB_sim",nyr_remove=1:4) # removing 5 years results in an error likely due to closed season commercial pot discards running out of data
#' run_retro("BluelineTilefish", dir_bam_base="BlTi_base", dir_bam_sim="BlTi_sim")
#' run_retro("Cobia", dir_bam_base="Cobi_base", dir_bam_sim="Cobi_sim")
#' run_retro("GagGrouper", dir_bam_base="GaGr_base", dir_bam_sim="GaGr_sim")
#' run_retro("GrayTriggerfish", dir_bam_base="GrTr_base", dir_bam_sim="GrTr_sim")
#' run_retro("GreaterAmberjack", dir_bam_base="GrAm_base", dir_bam_sim="GrAm_sim",nyr_remove=1:4) # removing 5 years results in an error because you run out of rGN discard length comps which start in 2013
#' run_retro("RedGrouper", dir_bam_base="ReGr_base", dir_bam_sim="ReGr_sim")
#' run_retro("RedPorgy", dir_bam_base="RePo_base", dir_bam_sim="RePo_sim")
#' run_retro("RedSnapper", dir_bam_base="ReSn_base", dir_bam_sim="ReSn_sim")
#' run_retro("ScampGrouper", dir_bam_base="ScGr_base", dir_bam_sim="ScGr_sim")
#' run_retro("SnowyGrouper", dir_bam_base="SnGr_base", dir_bam_sim="SnGr_sim")
#' run_retro("Tilefish", dir_bam_base="Tile_base", dir_bam_sim="Tile_sim")
#' run_retro("VermilionSnapper", dir_bam_base="VeSn_base", dir_bam_sim="VeSn_sim")
#'
#' }


run_retro <- function(CommonName = NULL,
                     fileName     = "bam",
                     dir_bam_sim  = "sim",
                     dir_bam_base = "base",
                     bam=NULL,
                     dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                     dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                     standardize=FALSE,
                     nyr_remove=1:5,
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
                     random_seed=12345,
                     admb2r_obj = admb2r.cpp,
                     cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                               "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                               "variance","*.dep","*.hes","*.tmp"))
){
######################
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

  endyr <- as.numeric(init$endyr)

  endyr_sim <- local({
    a <- as.numeric(endyr)-nyr_remove
    setNames(a,paste0("endyr",a))
  })

  nm_sim <- names(endyr_sim)
  nsim <- length(nm_sim)

  inits <- list()

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

  ### Identify temporal objects in init
  # styr
  init_styr <- init[grepl("^styr",names(init))]

  # endyr
  init_endyr <- init[grepl("^endyr",names(init))]

  # nyr
  init_nyr <- init[grepl("^nyr",names(init))]

  # yrs
  init_yrs <- init[grepl("^yrs",names(init))]

  # obs
  init_obs <- init[grepl("^obs_(cpue|L|cv|released|lenc|agec|maturity)",names(init))]

  # tv (time varying objects currently in Atlantic Menhaden model)
  init_tv <- init[grepl("_tv$",names(init))]

  # set_log_dev_vals (not Nage)
  init_set_log_dev_vals <- init[grepl("^set_log_dev_vals(?!.*(Nage))",names(init),perl=TRUE)]

  # nsamp
  init_nsamp <- init[grepl("^nsamp",names(init))]

  # nfish
  init_nfish <- init[grepl("^nfish",names(init))]

# For each value in endyr_sim..
  for(i in seq_along(endyr_sim)){
    init_i <- init
    endyr_i <- endyr_sim[i]
    nyr_remove_i <- nyr_remove[i]

    #### Trim temporal objects
    # styr values
    # (few styr values will be affected other than styr_regs which is used in the projection code)
    init_styr_i <- lapply(init_styr,function(x){
      if(as.numeric(x)>endyr_i){
        paste(as.numeric(x)-nyr_remove_i)
      }else{
        paste(x)
      }
    })

    # endyr values
    init_endyr_i <- setNames(lapply(seq_along(init_endyr),function(j){
      x <- init_endyr[j]
      xd <- as.numeric(x)-as.numeric(endyr_i)
      # For endyr, endyr_rec_, or endyr_proj values, always reduce the values by nyr_remove_i
      if(grepl("^endyr$|^endyr_rec_(?!.*(phase1))|endyr_proj",names(init_endyr[j]),perl=TRUE)){
        paste(as.numeric(x)-nyr_remove_i) # Reduce it by nyr_remove_i
        # or if the current _endyr value is greater than the new endyr_i..
        }else if(as.numeric(x)>endyr_i){
          as.character(as.numeric(endyr_i)) # Set it equal to endyr_i
        }else{
        paste(x)
        }
    }),
    names(init_endyr))

      # # Make sure endyr_rec values aren't larger than endyr_rec_dev
      # init_endyr_i[grepl("endyr_rec",names(init_endyr_i))] <-
      #   lapply(init_endyr_i[grepl("endyr_rec",names(init_endyr_i))],function(x){paste(min(init_endyr_i$endyr_rec_dev,as.numeric(x)))})

    # yrs values
    init_yrs_i <- lapply(init_yrs,function(x){
      x[which(x<=endyr_i)]
    })

    # nyr values
    init_nyr_i <- init_nyr # initialize
    for(j in seq_along(init_nyr)){
        nyr_nm_j <- names(init_nyr)[j]
        init_nyr_i[[j]] <- paste(length(init_yrs_i[[gsub("^nyr","yrs",names(init_nyr)[j])]]))
    }

    # obs values
    init_obs_i <- lapply(init_obs,function(x){
      if(is.vector(x)){
      y <- as.numeric(names(x))
      # If the names are years within the model data range..
      if(all(y%in%init$styr:init$endyr)){
        x[which(y<=endyr_i)] # ..truncate the object to the desired set of years
      }else{ # ..otherwise do nothing because the object doesn't have a time dimension
        x
      }
      }else if(is.matrix(x)){
        y <- as.numeric(rownames(x))
        x[paste(y)[which(y<=endyr_i)],,drop=FALSE]

      }else{
        warning("element x in init_obs_i is neither a vector nor a matrix")
      }

    })

    # set_log_dev_vals values
    init_set_log_dev_vals_i <- setNames(lapply(seq_along(init_set_log_dev_vals),function(j){
      name_j <- names(init_set_log_dev_vals)[j]
      x <- init_set_log_dev_vals[[j]]
      y <- as.numeric(names(x))
      if(name_j=="set_log_dev_vals_rec"){
        x[which(y<=init_endyr_i$endyr_rec_dev)]
      }else{
      x[which(y<=endyr_i)]
      }
    }),
    names(init_set_log_dev_vals)
    )

    # nsamp values
    init_nsamp_i <- lapply(init_nsamp,function(x){
      y <- as.numeric(names(x))
      x[which(y<=endyr_i)]
    })

    # nfish values
    init_nfish_i <- lapply(init_nfish,function(x){
      y <- as.numeric(names(x))
      x[which(y<=endyr_i)]
    })

    #### Add trimmed temporal objects back into init_i

    init_i[names(init_styr_i)] <- init_styr_i
    init_i[names(init_endyr_i)] <- init_endyr_i
    init_i[names(init_yrs_i)] <- init_yrs_i
    init_i[names(init_nyr_i)] <- init_nyr_i
    init_i[names(init_obs_i)] <- init_obs_i
    init_i[names(init_set_log_dev_vals_i)] <- init_set_log_dev_vals_i
    init_i[names(init_nsamp_i)] <- init_nsamp_i
    init_i[names(init_nfish_i)] <- init_nfish_i

    inits[[i]] <- init_i
  }
  names(inits) <- names(endyr_sim)


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

    message(paste("Running retrospective analysis for", nsim,"sims in parallel on",coresUse,"cores at",Sys.time()))


  sim_out <- foreach(i=1:nsim,
          #, .export = ls()[grepl("xboot.obs",ls())] # Pass additional objects to foreach
          .packages=c("bamExtras")
  ) %dopar% {
    nm_sim_i <- nm_sim[i]
    init_i <- inits[[nm_sim_i]]

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
    writeLines(text=bam_i$tpl, con=gsub("dat$","tpl",fileName_dat_i))
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
  } #end retrospective analysis foreach loop


  message(paste("Finished running retrospective analysis for", nsim,"sims in parallel on",coresUse,"cores at",Sys.time()))

  sim_out <- apply(do.call(rbind,sim_out),2,as.numeric)
  sim_out <- cbind("sim"=nm_sim,as.data.frame(sim_out)
                   # ,"Linf"=sim_Linf,"K"=sim_K, "t0"=sim_t0, "M"=sim_M, "Mage_scale"=sim_Masc,
                   #  "steep"=sim_steep,"rec_sigma"=sim_rec_sigma
                   )

  write.csv(x=sim_out,file=file.path(dir_bam_sim,"retro_results.csv"),row.names = FALSE)

### Return stuff
invisible(sim_out)

  } # end if(run_sim)

  ## parallel shut down
  stopCluster(cl)
}
