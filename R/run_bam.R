#' Run bam model
#'
#' Calls a shell script to run a bam model
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#'
#' }

run_bam <- function(CommonName = NULL, fileName = "bam", dir_bam = NULL,
                    bam=NULL,
                    dat_file=NULL,tpl_file=NULL,cxx_file=NULL,
                    dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,
                    admb_switch = '-nox',
                    admb2r_obj = admb2r.cpp,
                    cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                              "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                              "variance","*.dep","*.hes","*.tmp"))
){
  wd <- getwd()
  message(paste("working directory:",wd))

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

  if(is.null(dir_bam)){
    dir_bam <- ifelse(fileName=="bam",
                      fileName,
                      paste0("bam_",fileName)
    )
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
  run_command <- paste(fileName, admb_switch)
  message(paste("Running",fileName,"model"))
  shell(run_command)

  shell("cleanup.bat")

  rdat <- dget(fileName_rdat)
  setwd(wd)
  message(paste("Working directory changed back to:",wd))

  invisible(rdat)
}
