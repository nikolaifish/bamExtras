#' Build cxx from tpl
#'
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
#' @param unlink_dir_bam Should \code{dir_bam} be deleted after this function is run?
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' Run a bam model and assign rdat output to object
#' rdat_AtMe <- run_bam("AtlanticMenhaden")
#' rdat_BlSB <- run_bam("BlackSeaBass")
#' rdat_BlTi <- run_bam("BluelineTilefish")
#' rdat_Cobi <- run_bam("Cobia")
#' rdat_GagG <- run_bam("GagGrouper")
#' rdat_GrTr <- run_bam("GrayTriggerfish")
#' rdat_GrAm <- run_bam("GreaterAmberjack")
#' rdat_ReGr <- run_bam("RedGrouper")
#' rdat_RePo <- run_bam("RedPorgy")
#' rdat_ReSn <- run_bam("RedSnapper")
#' rdat_SnGr <- run_bam("SnowyGrouper")
#' rdat_Tile <- run_bam("Tilefish")
#' rdat_VeSn <- run_bam("VermilionSnapper")
#'
#' }

run_bam <- function(CommonName = NULL, fileName = "bam", dir_bam = NULL,
                    bam=NULL,
                    dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                    dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                    standardize=TRUE,
                    unlink_dir_bam=TRUE
){
  # wd <- getwd()
  # message(paste("working directory:",wd))

  if(!is.null(CommonName)){
    fileName <- CommonName
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

  # if(is.null(dir_bam)){
  #   dir_bam <- ifelse(fileName=="bam",
  #                     fileName,
  #                     paste0("bam_",fileName)
  #   )
  # }
  #
  # if(!dir.exists(dir_bam)) {
  #   dir.create(dir_bam)
  #   message(paste("Created directory:",dir_bam))
  # }
  # setwd(dir_bam)
  #
  # message(paste("Working directory temporarily changed to:",dir_bam))

  #### A bunch of regex ####
  ## Get various types of object names and also their dimensions
  tpl_parNames <- c()

  ## Build cxx code and return a cxx character vector

  regex <- "(?<=const double )(.*?)(?==)"
  regex_A <- "(?<="
  regex_Z <- " )(.*?)(?==)"

  tpl_object_classes <- c("const double")

  tpl_names <- list()

  for(class_i in tpl_object_classes){
    tpl_names[[class_i]] <- unlist(str_extract_all(tpl,paste0(regex_A,"const double", regex_Z)))
  }




  # setwd(wd)
  # message(paste("Working directory changed back to:",wd))
  #
  # if(unlink_dir_bam){
  #   # delete a directory -- must add recursive = TRUE
  #   unlink(dir_bam, recursive = TRUE)
  # }

  invisible(rdat)
}
