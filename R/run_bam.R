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
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param standardize Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM
#' @param subset_rdat list of rdat objects to subset and number of values to retain.
#' This option can substantially decrease rdat file size, without affecting precision of
#' reference point calculations.
#' @param unlink_dir_bam Should \code{dir_bam} be deleted after this function is run?
#' @param admb_options Character string pasted to fileName to build \code{run_command} when running BAM with \code{shell(run_command)}. See \href{https://www.admb-project.org/docs/refcards/admb-additional-reference-card.pdf}{ADMB reference card} for more options.
#' (i.e. \code{run_command <- paste(fileName, admb_options)})
#' @param admb2r_obj Character string containing admb2r C++ code, which is written with \code{base::writeLines} to \code{dir_bam}
#' @param cleanup List object written to \code{cleanup.bat} file in \code{dir_bam}.
#' @param return_obj objects to return from the function. May include one or more of the following "dat", "tpl", "cxx", "rdat". character vector.
#' @keywords bam stock assessment fisheries
#' @returns See \code{return_obj}
#' @export
#' @examples
#' \dontrun{
#' Run a bam model and assign rdat output to object
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
#' }

run_bam <- function(CommonName = NULL, fileName = "bam", dir_bam = NULL,
                    bam=NULL,
                    dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                    dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                    standardize=TRUE,
                    subset_rdat=list("eq.series"=101,"pr.series"=101),
                    unlink_dir_bam=TRUE,
                    admb_options = '-nox',
                    admb2r_obj = admb2r.cpp,
                    cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                              "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                              "variance","*.dep","*.hes","*.tmp")),
                    return_obj = "rdat"
){
  wdt <- getwd()
  message(paste("working directory:",wdt))

  if(!is.null(CommonName)){
    if(is.null(fileName)) {fileName <- CommonName}
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
  }

  if(is.null(dir_bam)){
    dir_bam <- file.path(wdt,fileName)
  }

  if(!dir.exists(dir_bam)) {
    dir.create(dir_bam)
    message(paste("Created directory:",dir_bam))
  }
  setwd(dir_bam)

  message(paste("Working directory temporarily changed to:",dir_bam))

  cleanup_bat <- unlist(lapply(1:length(cleanup),function(x){
    xname <- names(cleanup)[x]
    xx <- cleanup[[x]]
    paste(xname,xx)
  }))

  fileName_dat  <- paste0(fileName,".dat")
  fileName_tpl  <- paste0(fileName,".tpl")
  fileName_cxx  <- paste0(fileName,".cxx")
  fileName_rdat <- paste0(fileName,".rdat")

  # Change name of cxx file included in tpl to fileName_cxx
  lineno_tpl_cxx <- which(grepl("^\\s.*(#include).*(.cxx)",tpl))
  tpl[lineno_tpl_cxx] <- gsub("([a-zA-Z_0-9.])*\\.cxx",fileName_cxx, tpl[lineno_tpl_cxx])

  writeLines(text=dat, con=fileName_dat)
  writeLines(text=tpl, con=fileName_tpl)
  writeLines(text=cxx, con=fileName_cxx)
  writeLines(text=cleanup_bat, con="cleanup.bat")
  writeLines(text=admb2r_obj,con="admb2r.cpp")

  compile_command <- paste("admb", fileName)
  message(paste("Compiling",fileName,"model"))
  shell(compile_command)

  #######Run admb
  run_command <- paste(fileName, admb_options)
  message(paste("Running",fileName,"model"))
  shell(run_command)

  shell("cleanup.bat")

  rdat <- dget(fileName_rdat)

  # Subset rdat to decrease size on disk
  for(nm_k in names(subset_rdat)){
    k <- rdat[[nm_k]]
    if(is.null(subset_rdat[[nm_k]])){
      rdat[[nm_k]] <- NULL
    }else{
      rdat[[nm_k]] <- k[seq(1,nrow(k),length=subset_rdat[[nm_k]]),,drop=FALSE]
    }
  }

  setwd(wdt)
  message(paste("Working directory changed back to:",wdt))

  if(unlink_dir_bam){
    # delete a directory -- must add recursive = TRUE
    unlink(dir_bam, recursive = TRUE)
  }
  rm(wdt)

  bam_obj <- list("dat"=dat,"tpl"=tpl,"cxx"=cxx,"rdat"=rdat)
  out <- bam_obj[return_obj]
#
#   if(length(out)==1){
#     out <- out[[1]]
#   }

  invisible(out)
}
