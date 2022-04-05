#' Run bam model
#'
#' Calls a shell script to run a bam model
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Read in any of the current BAM models
#' bam_AtMe <- bam2r("AtlanticMenhaden")
#' bam_BlSB <- bam2r("BlackSeaBass")
#' bam_BlTi <- bam2r("BluelineTilefish")
#' bam_Cobi <- bam2r("Cobia")
#' bam_GagG <- bam2r("GagGrouper")
#' bam_GrTr <- bam2r("GrayTriggerfish")
#' bam_GrAm <- bam2r("GreaterAmberjack")
#' bam_ReGr <- bam2r("RedGrouper")
#' bam_RePo <- bam2r("RedPorgy")
#' bam_ReSn <- bam2r("RedSnapper")
#' bam_SnGr <- bam2r("SnowyGrouper")
#' bam_Tile <- bam2r("Tilefish")
#' bam_VeSn <- bam2r("VermilionSnapper")
#'
#' # Run a bam model and assign rdat output to object
#' rdat_RePo <- run_bam(bam=bam_RePo,fileName="RePo")
#'
#' # Modify data input from a BAM model and incorporate it back into the dat file object
#' L_init2 <- bam_RePo$L_init
#' L_init2$set_steep[c(1,5)] <- paste(0.6) # Change steepness
#' L_init2$set_steep[4] <- paste(-abs(as.numeric(L_init2$set_steep[4]))) # Fix steepness
#' bam_RePo2 <- bam2r("RedPorgy",L_init_user=L_init2)
#' rdat_RePo2 <- run_bam(bam=bam_RePo2,fileName="RePo2")
#'
#' # Compare models
#' plot(rdat_RePo$t.series$year, rdat_RePo$t.series$SSB.msst,type="o", ylim=range(c(rdat_RePo$t.series$SSB.msst,rdat_RePo2$t.series$SSB.msst),na.rm=TRUE), xlab="year",ylab="SSB.msst")
#' points(rdat_RePo2$t.series$year,rdat_RePo2$t.series$SSB.msst,type="o",col="blue")
#' legend("topright",legend=c("RePo","RePo2"),pch=1,lty=1,col=c("black","blue"),bty="n")
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
