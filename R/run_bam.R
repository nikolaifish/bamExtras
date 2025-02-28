#' Run bam model
#'
#' Creates a folder in a directory \code{dir_bam} to write BAM files including model-specific files
#' (.dat, .tpl, and .cxx) a standard \code{admb2r.cpp} file stored in \code{bamExtras}, and a standard
#' but customizable \code{cleanup.bat} file. Model-specific files are supplied by the user in any one of several different ways
#' (see Arguments). It temporarily changes the working directory to \code{dir_bam}
#' then calls a shell script to run the BAM model, and changes the working directory back to the previous path.
#' The function invisibly returns the bam file objects (indicated by return_obj) in a list.
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @param fileName Name given to BAM files, not including file extensions.
#' @param dir_bam Name of directory to write BAM files to, relative to the working directory.
#' @param bam Output of \code{bam2r}.
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @param exe_file exe file path. If you specify only exe_file and dat_file or dat_obj,
#' the function will copy the executable and will not try to compile bam. It will just
#' run teh executable with the dat file you give it. This option is used with
#' run_MCBE.
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param standardize Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM
#' @param subset_rdat list of rdat objects to subset in different ways. For eq.series and pr.series
#' specify number of evenly spaced values to retain. For t.series specify a series
#' of years to retain, either as a numeric vector, or as a character string that will be
#' evaluated as code internally with \code{eval(parse(text=my_string))}.
#' When subsetting eq.series or pr.series, this option can
#' substantially decrease rdat file size, without affecting precision of reference
#' point calculations. This is particularly helpful for the MCBE runs.
#' @param unlink_dir_bam Should \code{dir_bam} be deleted after this function is run?
#' @param admb_options Character string pasted to fileName to build \code{run_command} when running BAM with \code{shell(run_command)}. See \href{https://www.admb-project.org/docs/refcards/admb-additional-reference-card.pdf}{ADMB reference card} for more options.
#' (i.e. \code{run_command <- paste(fileName, admb_options)})
#' @param admb2r_obj Character string containing admb2r C++ code, which is written with \code{base::writeLines} to \code{dir_bam}
#' @param cleanup List object written to \code{cleanup.bat} file in \code{dir_bam}.
#' @param run_cleanup Should the \code{cleanup.bat} be run?
#' @param return_obj names of objects to return from the function. May include one or more
#' of the following values which refer to ADMB file extensions:
#' \code{dat}, \code{tpl}, \code{cxx}, \code{rdat} or \code{par}. These files will be
#' read with \code{readLines} and passed to the result as character vectors. May
#' also include \code{admb} which returns a list of basic run results including
#' \code{lk_total} (total likelihood), \code{grad_max} (maximum gradient), and
#' \code{file_names} (list of file names in \code{dir_bam} after running
#' \code{\link[base]{shell}}).
#' Character vector. If NULL, all objects are returned.
#' @keywords bam stock assessment fisheries
#' @returns See \code{return_obj}
#' @export
#' @examples
#' \dontrun{
#' # Run a bam model and assign rdat output to object
#' rdat_AtMe <- run_bam("AtlanticMenhaden")$rdat
#' rdat_BlSB <- run_bam("BlackSeaBass")$rdat
#' rdat_BlTi <- run_bam("BluelineTilefish")$rdat
#' rdat_Cobi <- run_bam("Cobia")$rdat
#' rdat_GagG <- run_bam("GagGrouper")$rdat
#' rdat_GrTr <- run_bam("GrayTriggerfish")$rdat
#' rdat_GrAm <- run_bam("GreaterAmberjack")$rdat
#' rdat_ReGr <- run_bam("RedGrouper")$rdat
#' rdat_RePo <- run_bam("RedPorgy")$rdat
#' rdat_ReSn <- run_bam("RedSnapper")$rdat
#' rdat_ScGr <- run_bam("ScampGrouper")$rdat
#' rdat_SnGr <- run_bam("SnowyGrouper")$rdat
#' rdat_Tile <- run_bam("Tilefish")$rdat
#' rdat_VeSn <- run_bam("VermilionSnapper")$rdat
#'
#' ### Change a value in the init object, and run bam
#' bam_RePo <- bam_RePo_raw <- bam2r("RedPorgy")
#' bam_RePo$init$obs_cpue_sCT[1] <- paste(2*as.numeric(bam_RePo$init$obs_cpue_sCT[1]))
#' bam_RePo <- bam2r(init=bam_RePo$init,dat_obj = bam_RePo$dat,tpl_obj = bam_RePo$tpl,cxx_obj = bam_RePo$cxx)
#' bamout_RePo <- run_bam(bam=bam_RePo,return_obj = NULL,
#'                        shell_args_compile = list("intern"=TRUE),
#'                        shell_args_run = list("intern"=TRUE)
#' )
#'
#' ##### Create and run examples where errors occur #####
#' ### silly value in dat file. Data file NOT READ CORRECTLY.
#' # Note that although this run fails, it is not a shell error (i.e. shell error code = 0)
#' bam_RePo_err1 <- bam_RePo_raw
#' bam_RePo_err1$init$endyr <- "pickles"
#' bam_RePo_err1 <- bam2r(init=bam_RePo_err1$init,dat_obj = bam_RePo_err1$dat,tpl_obj = bam_RePo_err1$tpl,cxx_obj = bam_RePo_err1$cxx)
#' bamout_RePo_err1 <- run_bam(bam=bam_RePo_err1,return_obj = NULL,
#'                             shell_args_compile = list("intern"=TRUE),
#'                             shell_args_run = list("intern"=TRUE)
#' )
#' # look at results to see the Data File reading error
#' bamout_RePo_err1$admb$run_out[grepl("Data File",bamout_RePo_err1$admb$run_out)]
#'
#' ### Negative value in cpue. Data file read correctly, but shell(run_command) fails with error code 1
#' bam_RePo_err2 <- bam_RePo_raw
#' bam_RePo_err2$init$obs_cpue_sCT[1] <- paste(-10)
#' bam_RePo_err2 <- bam2r(init=bam_RePo_err2$init,dat_obj = bam_RePo_err2$dat,tpl_obj = bam_RePo_err2$tpl,cxx_obj = bam_RePo_err2$cxx)
#' bamout_RePo_err2 <- run_bam(bam=bam_RePo_err2,return_obj = NULL,
#'                             shell_args_compile = list("intern"=TRUE),
#'                             shell_args_run = list("intern"=TRUE)
#' )
#' # look at results to see that the Data File was read correctly
#' bamout_RePo_err2$admb$run_out[grepl("Data File",bamout_RePo_err2$admb$run_out)]
#' # get run error code
#' bamout_RePo_err2$admb$error_code_run
#' # also note that lk_total is NaN
#' bamout_RePo_err2$admb$lk_total
#'
#' ### Crazy values of M. Data file read correctly, no shell error, but Hessian doesn't converge
#' bam_RePo_err3 <- bam_RePo_raw
#' bam_RePo_err3$init$set_M <- setNames(paste(as.numeric(bam_RePo_err3$init$set_M)*10),names(bam_RePo_err3$init$set_M))
#' bam_RePo_err3 <- bam2r(init=bam_RePo_err3$init,dat_obj = bam_RePo_err3$dat,tpl_obj = bam_RePo_err3$tpl,cxx_obj = bam_RePo_err3$cxx)
#' bamout_RePo_err3 <- run_bam(bam=bam_RePo_err3,return_obj = NULL,
#'                             shell_args_compile = list("intern"=TRUE),
#'                             shell_args_run = list("intern"=TRUE)
#' )
#' # look at results to see the indefinite Hessian error
#' bamout_RePo_err3$admb$run_out[grepl("Hessian",bamout_RePo_err3$admb$run_out)]
#' }

run_bam <- function(CommonName = NULL, fileName = "bam", dir_bam = NULL,
                    bam=NULL,
                    dat_file=NULL,tpl_file=NULL,cxx_file=NULL,exe_file=NULL,
                    dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                    standardize=TRUE,
                    subset_rdat=list("eq.series"=101,"pr.series"=101,"t.series"="styr:endyr"),
                    unlink_dir_bam=TRUE,
                    admb_options = '-nox',
                    admb2r_obj = admb2r.cpp,
                    cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                              "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                              "variance","*.dep","*.hes","*.tmp")),
                    run_cleanup=TRUE,
                    shell_args_compile = list("intern"=FALSE),
                    shell_args_run = list("intern"=FALSE),
                    return_obj = "all"
){
  dat <- tpl <- cxx <- NULL # initialize as NULL
  is_exe <- !is.null(exe_file)
  if(is_exe){ # Read in dat and copy exe to working directory
    if(!is.null(dat_obj)){
      dat <- dat_obj
      message("exe_file and dat_obj specified. exe and dat will be used to run bam")
    }else if(!is.null(dat_file)){
      message("exe_file and dat_file specified. exe and dat will be used to run bam")
      dat <- readLines(con=dat_file)
    }else{
      stop("Can't run bam: exe_file specified but dat not found")
    }
    if(!is.null(tpl_obj)){
      tpl <- tpl_obj
      message("tpl_obj specified. tpl object will be stored in output, but will not be used to run bam")
    }else if(!is.null(tpl_file)){
      tpl <- readLines(con=tpl_file)
      message("tpl_file specified. tpl object will be stored in output, but will not be used to run bam")
    }else{
      tpl <- NULL
    }
    if(!is.null(cxx_obj)){
      cxx <- cxx_obj
      message("cxx_obj specified. cxx object will be stored in output, but will not be used to run bam")
    }else if(!is.null(cxx_file)){
      cxx <- readLines(con=cxx_file)
      message("cxx_file specified. cxx object will be stored in output, but will not be used to run bam")
    }else{
      cxx <- NULL
    }
  }else if(!is.null(CommonName)){ # Use stored files associated with CommonName
    if(is.null(fileName)) {fileName <- CommonName}
    nm_dat <- paste0("dat_",CommonName)
    nm_tpl <- paste0("tpl_",CommonName)
    nm_cxx <- paste0("cxx_",CommonName)
    dat <- get(nm_dat)
    tpl <- get(nm_tpl)
    cxx <- get(nm_cxx)
    message(paste0("CommonName specified. ",nm_dat,", ",nm_tpl,", and ",nm_cxx," will be used to run bam"))
  }else if(!is.null(bam)){ # Use parts from bam object
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
    message(paste0("bam object specified. dat, tpl, and cxx will be used to run bam"))
  }else if(!is.null(dat_obj)&!is.null(tpl_obj)&!is.null(cxx_obj)){ # Use dat, tpl, and cxx objects
    dat <- dat_obj
    tpl <- tpl_obj
    cxx <- cxx_obj
    message(paste0("dat_obj, tpl_obj, and cxx_obj specified. dat, tpl, and cxx will be used to run bam"))
  }else if(!is.null(dat_file)&!is.null(tpl_file)&!is.null(cxx_file)){ # Read in dat, tpl, and cxx files
    dat <- readLines(con=dat_file)
    tpl <- readLines(con=tpl_file)
    cxx <- readLines(con=cxx_file)
    message(paste0("dat_file, tpl_file, and cxx_file specified. dat, tpl, and cxx will be used to run bam"))
  }else{
    stop("Can't run bam: insufficient information supplied")
  }

wdt <- getwd()
message(paste("working directory:",wdt))

  if(standardize&!is_exe){
    message("Running bamExtras::standardize_bam()")
    bam <- standardize_bam(dat_obj=dat, tpl_obj=tpl,cxx_obj=cxx)
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
  }

  if(is.null(dir_bam)){
    dir_bam <- file.path(wdt,fileName)
  }

  if(!dir.exists(dir_bam)) {
    dir.create(dir_bam)
    message(paste("Created directory:",dir_bam))
  }
  if(is_exe){
    # Copy specified executable
    fileName_exe  <- paste0(fileName,".exe")
    fileName_exe_to <- file.path(dir_bam,fileName_exe)
    message(paste("exe_file specified. \nCopying",exe_file,"to:\n",fileName_exe_to))
    file.copy(exe_file, fileName_exe_to, overwrite=TRUE)
  }

  setwd(dir_bam)
  message(paste("Working directory temporarily changed to:\n",dir_bam))

  if(run_cleanup){
    cleanup_bat <- unlist(lapply(1:length(cleanup),function(x){
      xname <- names(cleanup)[x]
      xx <- cleanup[[x]]
      paste(xname,xx)
    }))
    message("cleanup.bat written to current directory")
    writeLines(text=cleanup_bat, con="cleanup.bat")
  }

  fileName_dat  <- paste0(fileName,".dat")
  fileName_tpl  <- paste0(fileName,".tpl")
  fileName_cxx  <- paste0(fileName,".cxx")
  fileName_rdat <- paste0(fileName,".rdat")
  fileName_par  <- paste0(fileName,".par")

  # Change name of cxx file included in tpl to fileName_cxx
  if(!is_exe){
    lineno_tpl_cxx <- which(grepl("^\\s.*(#include).*(.cxx)",tpl))
    tpl[lineno_tpl_cxx] <- gsub("([a-zA-Z_0-9.])*\\.cxx",fileName_cxx, tpl[lineno_tpl_cxx])
  }

  if(!is.null(dat)){writeLines(text=dat, con=fileName_dat)}
  if(!is.null(tpl)){writeLines(text=tpl, con=fileName_tpl)}
  if(!is.null(cxx)){writeLines(text=cxx, con=fileName_cxx)}

     if(!is_exe){
       # Compile bam if an executable wasn't supplied
       writeLines(text=admb2r_obj,con="admb2r.cpp")
       compile_command <- paste("admb", fileName)
       message(paste("Compiling",fileName,"model"))
       compile_out <- do.call(shell,c(list("cmd"=compile_command),shell_args_compile))
     }else{
       compile_out <- NULL
     }

  #######Run admb
  run_command <- paste(fileName, admb_options)
  message(paste("Running",fileName,"model"))
  run_out <- do.call(shell,c(list("cmd"=run_command),shell_args_run))

  # get error_code_run for run command
  error_code_run <- if(!is.null(shell_args_run$intern)&shell_args_run$intern){
    # If shell_args_run$intern==TRUE, then shell() returns the output of cmd..
    attr_ro <- attributes(run_out)
    if(!is.null(attr_ro)){
    # If there is an attribute, it is probably 'status' and it's an error (see base::shell help)
      attr_ro$status
    }else{
      # If there is no 'status' attribute, there was no error.
      0
    }
  }else{
    # ..if not, it returns a straight error code
    run_out
  }

  nm_files <- list.files()
  is_par_file <- fileName_par%in%nm_files
  is_rdat_file <- fileName_rdat%in%nm_files

  if(is_par_file){
  par_obj <- readLines(fileName_par)
  lk_total <- gsub("^.*Objective function value = (.*)  Maximum.*$","\\1",par_obj[1])
  grad_max <- gsub("^.*component = ","\\1",par_obj[1])
  }else{
    par_obj <- lk_total <- grad_max <- NA
  }

  admb_out <- list("lk_total"=as.numeric(lk_total),"grad_max"=as.numeric(grad_max),
                   "file_names"=nm_files,"compile_out"=compile_out, "run_out"=run_out,
                   "error_code_run"=error_code_run)

  if(is_rdat_file&is.finite(as.numeric(lk_total))&is.finite(as.numeric(grad_max))&error_code_run==0){
    rdat <- dget(fileName_rdat)
    styr <- rdat$parms$styr
    endyr <- rdat$parms$endyr

    # Subset rdat to decrease size on disk
    for(nm_k in names(subset_rdat)){
      k_sub_val <- subset_rdat[[nm_k]]
      if(!is.null(k_sub_val)){
        k <- rdat[[nm_k]]
        if(nm_k%in%c("eq.series","pr.series")){
          rdat[[nm_k]] <- k[seq(1,nrow(k),length=k_sub_val),,drop=FALSE]
        }else if(nm_k=="t.series"){
          if(is.character(k_sub_val)){ # if it's a character string, try to evaluate it as code
            ts_yr_sub <- eval(parse(text=k_sub_val))
          }else if(is.numeric(k_sub_val)){ # if it's a numeric vector, try to use it directly to subset the years
            ts_yr_sub <- k_sub_val
          }
          rdat[[nm_k]] <- k[paste(ts_yr_sub),]
        }else{
          message(paste("no subsetting method is currently specified for",nm_k,"so no change was made"))
        }
      }
    }
    # Save file over previous
    dput(x=rdat,file=fileName_rdat)
  }else{
    rdat <- NULL
  }

  if(run_cleanup){
    shell("cleanup.bat")
  }

  bam_obj <- list("dat"=dat,"tpl"=tpl,"cxx"=cxx,"rdat"=rdat,"par"=par_obj,"admb"=admb_out)

  out <- if(!is.null(return_obj)){
    if("all"%in%return_obj){
      bam_obj
    }else{
      bam_obj[return_obj]
    }
  }else{
    bam_obj
  }


  setwd(wdt)
  message(paste("Working directory changed back to:\n",wdt))

  if(unlink_dir_bam){
    # delete a directory -- must add recursive = TRUE
    unlink(dir_bam, recursive = TRUE)
  }

  rm(wdt)

  invisible(out)
}
