#' Standardize object names in bam files (dat, tpl, and cxx)
#'
#' This script reads in BAM dat, tpl, and cxx files, and converts them to R objects (character vectors of each line of the code). It then runs a set of gsub() functions to
#' replace object names with preferred naming conventions and returns dat, tpl, and cxx character vectors to use to rewrite the corresponding files with writeLines)
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' standardize_bam()

standardize_bam <- function(dat_file,tpl_file,cxx_file){
  # Read in dat, tpl, and cxx files
  dat <- readLines(con=dat_file)
  dat <- trimws(dat) # Remove leading and/or trailing whitespace from character strings
  tpl <- readLines(con=tpl_file)
  cxx <- readLines(con=cxx_file)

  bamr <- bam_to_r(
    dat_file = dat_file,
    tpl_file = tpl_file,
    cxx_file = cxx_file)

  dat <- bamr$dat
  dat <- trimws(dat) # Remove leading and/or trailing whitespace from character strings
  tpl <- bamr$tpl
  cxx <- bamr$cxx
  L_init <- bamr$L_init

  # Restructure variable names so that fleet name is at the end (this makes it easier to identify similar variables in search)
  initNames <- names(L_init)

  ## Landings
  # styr_L
  nm1_styr_L <- initNames[grepl("^styr_.*_L$",initNames)]
  nm2_styr_L <- paste0("styr_L_",gsub("^styr_|_L$","",nm1_styr_L))

  # endyr_L
  nm1_endyr_L <- initNames[grepl("^endyr_.*_L$",initNames)]
  nm2_endyr_L <- paste0("endyr_L_",gsub("^endyr_|_L$","",nm1_endyr_L))

  # obs_L
  nm1_obs_L <- initNames[grepl("^obs_.*_L$",initNames)]
  nm2_obs_L <- paste0("obs_L_",gsub("^obs_|_L$","",nm1_obs_L))

  # obs_cv_L (landings cvs)
  nm1_obs_cv_L <- initNames[grepl("_L_cv",initNames)]
  nm2_obs_cv_L <- paste0("obs_cv_L_",gsub("_L_cv","",nm1_obs_cv_L))

  ## Discards
  # styr_D
  nm1_styr_D <- initNames[grepl("^styr_.*_D$",initNames)]
  nm2_styr_D <- paste0("styr_D_",gsub("^styr_|_D$","",nm1_styr_D))

  # endyr_D
  nm1_endyr_D <- initNames[grepl("^endyr_.*_D$",initNames)]
  nm2_endyr_D <- paste0("endyr_D_",gsub("^endyr_|_D$","",nm1_endyr_D))

  # obs_D
  nm1_obs_D <- initNames[grepl("^obs_.*_D$",initNames)]
  nm2_obs_D <- paste0("obs_D_",gsub("^obs_|_D$","",nm1_obs_D))

  # obs_cv_D (landings cvs)
  nm1_obs_cv_D <- initNames[grepl("_D_cv",initNames)]
  nm2_obs_cv_D <- paste0("obs_cv_D_",gsub("_D_cv","",nm1_obs_cv_D))

  ## CPUE
  # styr_cpue
  nm1_styr_cpue <- initNames[grepl("^styr_.*_cpue$",initNames)]
  nm2_styr_cpue <- paste0("styr_cpue_",gsub("^styr_|_cpue$","",nm1_styr_cpue))

  # endyr_cpue
  nm1_endyr_cpue <- initNames[grepl("^endyr_.*_cpue$",initNames)]
  nm2_endyr_cpue <- paste0("endyr_cpue_",gsub("^endyr_|_cpue$","",nm1_endyr_cpue))


  # nyr_cpue (not found in most assessments)
  nm1_nyr_cpue <- initNames[grepl("^nyr_.*_cpue$",initNames)]
  nm2_nyr_cpue <- paste0("nyr_cpue_",gsub("^nyr_|_cpue$","",nm1_nyr_cpue))

  # yrs_cpue (not found in most assessments)
  nm1_yrs_cpue <- initNames[grepl("^yrs_.*_cpue$",initNames)]
  nm2_yrs_cpue <- paste0("yrs_cpue_",gsub("^yrs_|_cpue$","",nm1_yrs_cpue))

  # obs_cpue
  nm1_obs_cpue <- initNames[grepl("^obs_.*_cpue$",initNames)]
  nm2_obs_cpue <- paste0("obs_cpue_",gsub("^obs_|_cpue$","",nm1_obs_cpue))

  # obs_cv_cpue (cpue cvs)
  nm1_obs_cv_cpue <- initNames[grepl("_cpue_cv",initNames)]
  nm2_obs_cv_cpue <- paste0("obs_cv_cpue_",gsub("_cpue_cv","",nm1_obs_cv_cpue))


  ## Length comps
  # nyr_lenc
  nm1_nyr_lenc <- initNames[grepl("^nyr_.*_lenc$",initNames)]
  nm2_nyr_lenc <- paste0("nyr_lenc_",gsub("^nyr_|_lenc$","",nm1_nyr_lenc))

  # yrs_lenc
  nm1_yrs_lenc <- initNames[grepl("^yrs_.*_lenc$",initNames)]
  nm2_yrs_lenc <- paste0("yrs_lenc_",gsub("^yrs_|_lenc$","",nm1_yrs_lenc))

  # nsamp_lenc
  nm1_nsamp_lenc <- initNames[grepl("^nsamp_.*_lenc$",initNames)]
  nm2_nsamp_lenc <- paste0("nsamp_lenc_",gsub("^nsamp_|_lenc$","",nm1_nsamp_lenc))

  # nfish_lenc
  nm1_nfish_lenc <- initNames[grepl("^nfish_.*_lenc$",initNames)]
  nm2_nfish_lenc <- paste0("nfish_lenc_",gsub("^nfish_|_lenc$","",nm1_nfish_lenc))

  # obs_lenc
  nm1_obs_lenc <- initNames[grepl("^obs_.*_lenc$",initNames)]
  nm2_obs_lenc <- paste0("obs_lenc_",gsub("^obs_|_lenc$","",nm1_obs_lenc))

  # pred_lenc (not in dat, but in tpl and cxx)
  nm1_pred_lenc <- initNames[grepl("^pred_.*_lenc$",initNames)]
  nm2_pred_lenc <- paste0("pred_lenc_",gsub("^pred_|_lenc$","",nm1_pred_lenc))


  # lenc (rearrange and replace remaining variable names)
  nm1_lenc <- paste(gsub("obs_","",nm1_obs_lenc))
  nm2_lenc <- paste("lenc",gsub("_lenc","",nm1_lenc),sep="_")

  ## Length comps pooled (not in most assessments)
  # nyr_lenc_pool
  nm1_nyr_lenc_pool <- initNames[grepl("^nyr_.*_lenc_pool$",initNames)]
  nm2_nyr_lenc_pool <- paste0("nyr_lenc_pool_",gsub("^nyr_|_lenc_pool$","",nm1_nyr_lenc_pool))

  # yrs_lenc_pool
  nm1_yrs_lenc_pool <- initNames[grepl("^yrs_.*_lenc_pool$",initNames)]
  nm2_yrs_lenc_pool <- paste0("yrs_lenc_pool_",gsub("^yrs_|_lenc_pool$","",nm1_yrs_lenc_pool))

  # nsamp_lenc_pool
  nm1_nsamp_lenc_pool <- initNames[grepl("^nsamp_.*_lenc_pool$",initNames)]
  nm2_nsamp_lenc_pool <- paste0("nsamp_lenc_pool_",gsub("^nsamp_|_lenc_pool$","",nm1_nsamp_lenc_pool))

  # nfish_lenc_pool
  nm1_nfish_lenc_pool <- initNames[grepl("^nfish_.*_lenc_pool$",initNames)]
  nm2_nfish_lenc_pool <- paste0("nfish_lenc_pool_",gsub("^nfish_|_lenc_pool$","",nm1_nfish_lenc_pool))

  # obs_lenc_pool
  nm1_obs_lenc_pool <- initNames[grepl("^obs_.*_lenc_pool$",initNames)]
  nm2_obs_lenc_pool <- paste0("obs_lenc_pool_",gsub("^obs_|_lenc_pool$","",nm1_obs_lenc_pool))

  # pred_lenc_pool (not in dat, but in tpl and cxx)
  nm1_pred_lenc_pool <- initNames[grepl("^pred_.*_lenc_pool$",initNames)]
  nm2_pred_lenc_pool <- paste0("pred_lenc_pool_",gsub("^pred_|_lenc_pool$","",nm1_pred_lenc_pool))

  ## Age comps
  # nyr_agec
  nm1_nyr_agec <- initNames[grepl("^nyr_.*_agec$",initNames)]
  nm2_nyr_agec <- paste0("nyr_agec_",gsub("^nyr_|_agec$","",nm1_nyr_agec))

  # yrs_agec
  nm1_yrs_agec <- initNames[grepl("^yrs_.*_agec$",initNames)]
  nm2_yrs_agec <- paste0("yrs_agec_",gsub("^yrs_|_agec$","",nm1_yrs_agec))

  # nsamp_agec
  nm1_nsamp_agec <- initNames[grepl("^nsamp_.*_agec$",initNames)]
  nm2_nsamp_agec <- paste0("nsamp_agec_",gsub("^nsamp_|_agec$","",nm1_nsamp_agec))

  # nfish_agec
  nm1_nfish_agec <- initNames[grepl("^nfish_.*_agec$",initNames)]
  nm2_nfish_agec <- paste0("nfish_agec_",gsub("^nfish_|_agec$","",nm1_nfish_agec))

  # obs_agec
  nm1_obs_agec <- initNames[grepl("^obs_.*_agec$",initNames)]
  nm2_obs_agec <- paste0("obs_agec_",gsub("^obs_|_agec$","",nm1_obs_agec))

  # pred_agec (not in dat, but in tpl and cxx)
  nm1_pred_agec <- initNames[grepl("^pred_.*_agec$",initNames)]
  nm2_pred_agec <- paste0("pred_agec_",gsub("^pred_|_agec$","",nm1_pred_agec))

  # agec (rearrange and replace remaining variable names)
  nm1_agec <- paste(gsub("obs_","",nm1_obs_agec))
  nm2_agec <- paste("agec",gsub("_agec","",nm1_agec),sep="_")

  ## Age comps pooled (not in most assessments)
  # nyr_agec_pool
  nm1_nyr_agec_pool <- initNames[grepl("^nyr_.*_agec_pool$",initNames)]
  nm2_nyr_agec_pool <- paste0("nyr_agec_pool_",gsub("^nyr_|_agec_pool$","",nm1_nyr_agec_pool))

  # yrs_agec_pool
  nm1_yrs_agec_pool <- initNames[grepl("^yrs_.*_agec_pool$",initNames)]
  nm2_yrs_agec_pool <- paste0("yrs_agec_pool_",gsub("^yrs_|_agec_pool$","",nm1_yrs_agec_pool))

  # nsamp_agec_pool
  nm1_nsamp_agec_pool <- initNames[grepl("^nsamp_.*_agec_pool$",initNames)]
  nm2_nsamp_agec_pool <- paste0("nsamp_agec_pool_",gsub("^nsamp_|_agec_pool$","",nm1_nsamp_agec_pool))

  # nfish_agec_pool
  nm1_nfish_agec_pool <- initNames[grepl("^nfish_.*_agec_pool$",initNames)]
  nm2_nfish_agec_pool <- paste0("nfish_agec_pool_",gsub("^nfish_|_agec_pool$","",nm1_nfish_agec_pool))

  # obs_agec_pool
  nm1_obs_agec_pool <- initNames[grepl("^obs_.*_agec_pool$",initNames)]
  nm2_obs_agec_pool <- paste0("obs_agec_pool_",gsub("^obs_|_agec_pool$","",nm1_obs_agec_pool))

  # pred_agec_pool (not in dat, but in tpl and cxx)
  nm1_pred_agec_pool <- initNames[grepl("^pred_.*_agec_pool$",initNames)]
  nm2_pred_agec_pool <- paste0("pred_agec_pool_",gsub("^pred_|_agec_pool$","",nm1_pred_agec_pool))

  ## set_ parameters
  # set_log_dm_lenc
  nm1_set_log_dm_lenc <- initNames[grepl("^set_log_dm_.*_lc$",initNames)]
  nm2_set_log_dm_lenc <- paste0("set_log_dm_lenc_",gsub("^set_log_dm_|_lc$","",nm1_set_log_dm_lenc))

  # set_log_dm_agec
  nm1_set_log_dm_agec <- initNames[grepl("^set_log_dm_.*_ac$",initNames)]
  nm2_set_log_dm_agec <- paste0("set_log_dm_agec_",gsub("^set_log_dm_|_ac$","",nm1_set_log_dm_agec))

  # set_log_q_cpue
  nm1_set_log_q_cpue <- initNames[grepl("^set_log_q_|set_logq_",initNames)]
  nm2_set_log_q_cpue <- paste0("set_log_q_cpue_",gsub("^set_log_q_|set_logq_","",nm1_set_log_q_cpue))

  # F pars for landings (L)
  # set_log_avg_F_L
  nm1_set_log_avg_F_L <- initNames[grepl("^set_log_avg_F_.*[^D]$",initNames)]
  nm2_set_log_avg_F_L <- paste0("set_log_avg_F_L_",gsub("^set_log_avg_F_","",nm1_set_log_avg_F_L))

  # set_log_F_dev_vals_L
  nm1_set_log_F_dev_vals_L <- initNames[grepl("set_log_F_dev_.*[^_D]_vals$",initNames)]
  nm2_set_log_F_dev_vals_L <- paste0("set_log_F_dev_vals_L_",gsub("^set_log_F_dev_|_vals$","",nm1_set_log_F_dev_vals_L))

  # set_log_F_dev_L
  nm1_set_log_F_dev_L <- initNames[grepl("^set_log_F_dev_.*[^_vals|_D]$",initNames)]
  nm2_set_log_F_dev_L <- paste0("set_log_F_dev_L_",gsub("^set_log_F_dev_","",nm1_set_log_F_dev_L))

  # F pars for discards (D)
  # set_log_avg_F_D
  nm1_set_log_avg_F_D <- initNames[grepl("^set_log_avg_F_.*_D$",initNames)]
  nm2_set_log_avg_F_D <- paste0("set_log_avg_F_D_",gsub("^set_log_avg_F_|_D$","",nm1_set_log_avg_F_D))

  # set_log_F_dev_vals_D
  nm1_set_log_F_dev_vals_D <- initNames[grepl("set_log_F_dev_.*_D_vals$",initNames)]
  nm2_set_log_F_dev_vals_D <- paste0("set_log_F_dev_vals_D_",gsub("^set_log_F_dev_|_D_vals$","",nm1_set_log_F_dev_vals_D))

  # set_log_F_dev_D
  nm1_set_log_F_dev_D <- initNames[grepl("^set_log_F_dev_.*_D$",initNames)]
  nm2_set_log_F_dev_D <- paste0("set_log_F_dev_D_",gsub("^set_log_F_dev_|_D$","",nm1_set_log_F_dev_D))

  ## weight parameters
  # set_w_cpue
  nm1_set_w_cpue <- initNames[grepl("^set_w_I_",initNames)]
  nm2_set_w_cpue <- paste0("set_w_cpue_",gsub("^set_w_I_","",nm1_set_w_cpue))

  # set_w_lenc
  nm1_set_w_lenc <- initNames[grepl("^set_w_lc_",initNames)]
  nm2_set_w_lenc <- paste0("set_w_lenc_",gsub("^set_w_lc_","",nm1_set_w_lenc))

  # set_p_lenc (not in most assessments. See Black Sea Bass)
  nm1_set_p_lenc <- initNames[grepl("^set_p_lenc_",initNames)]
  nm2_set_p_lenc <- paste0("set_p_lenc_",gsub("^set_p_lenc_","",nm1_set_p_lenc))

  # set_w_agec
  nm1_set_w_agec <- initNames[grepl("^set_w_ac_",initNames)]
  nm2_set_w_agec <- paste0("set_w_agec_",gsub("^set_w_ac_","",nm1_set_w_agec))

  # Combine all nm1 and nm2
  nm1 <- list(nm1_styr_L,nm1_endyr_L,nm1_obs_L,nm1_obs_cv_L,
              nm1_styr_D,nm1_endyr_D,nm1_obs_D,nm1_obs_cv_D,
              nm1_styr_cpue,nm1_endyr_cpue,nm1_nyr_cpue,nm1_yrs_cpue,nm1_obs_cpue,nm1_obs_cv_cpue,
              nm1_nyr_lenc, nm1_yrs_lenc ,nm1_nsamp_lenc, nm1_nfish_lenc, nm1_obs_lenc, nm1_pred_lenc, nm1_lenc,
              nm1_nyr_lenc_pool, nm1_yrs_lenc_pool ,nm1_nsamp_lenc_pool, nm1_nfish_lenc_pool, nm1_obs_lenc_pool, nm1_pred_lenc_pool,
              nm1_nyr_agec_pool, nm1_yrs_agec_pool ,nm1_nsamp_agec_pool, nm1_nfish_agec_pool, nm1_obs_agec_pool, nm1_pred_agec_pool,
              nm1_nyr_agec, nm1_yrs_agec ,nm1_nsamp_agec, nm1_nfish_agec, nm1_obs_agec, nm1_pred_agec, nm1_agec,
              nm1_set_log_dm_lenc,nm1_set_log_dm_agec, nm1_set_log_q_cpue,
              nm1_set_log_avg_F_L, nm1_set_log_F_dev_vals_L, nm1_set_log_F_dev_L,
              nm1_set_log_avg_F_D, nm1_set_log_F_dev_vals_D, nm1_set_log_F_dev_D,
              nm1_set_w_cpue,nm1_set_w_lenc,nm1_set_p_lenc,nm1_set_w_agec)

  nm2 <- list(nm2_styr_L,nm2_endyr_L,nm2_obs_L,nm2_obs_cv_L,
              nm2_styr_D,nm2_endyr_D,nm2_obs_D,nm2_obs_cv_D,
              nm2_styr_cpue,nm2_endyr_cpue,nm2_nyr_cpue,nm2_yrs_cpue,nm2_obs_cpue,nm2_obs_cv_cpue,
              nm2_nyr_lenc, nm2_yrs_lenc ,nm2_nsamp_lenc, nm2_nfish_lenc, nm2_obs_lenc, nm2_pred_lenc, nm2_lenc,
              nm2_nyr_lenc_pool, nm2_yrs_lenc_pool ,nm2_nsamp_lenc_pool, nm2_nfish_lenc_pool, nm2_obs_lenc_pool, nm2_pred_lenc_pool,
              nm2_nyr_agec_pool, nm2_yrs_agec_pool ,nm2_nsamp_agec_pool, nm2_nfish_agec_pool, nm2_obs_agec_pool, nm2_pred_agec_pool,
              nm2_nyr_agec, nm2_yrs_agec ,nm2_nsamp_agec, nm2_nfish_agec, nm2_obs_agec, nm2_pred_agec, nm2_agec,
              nm2_set_log_dm_lenc,nm2_set_log_dm_agec, nm2_set_log_q_cpue,
              nm2_set_log_avg_F_L, nm2_set_log_F_dev_vals_L, nm2_set_log_F_dev_L,
              nm2_set_log_avg_F_D, nm2_set_log_F_dev_vals_D, nm2_set_log_F_dev_D,
              nm2_set_w_cpue,nm2_set_w_lenc,nm2_set_p_lenc,nm2_set_w_agec)

  # Identify which elements on nm1 were found in initNames and remove those elements from nm1 and nm2
  ef <- unlist(lapply(nm1,function(x){length(x)!=0}))
  nm1 <- unlist(nm1[ef])
  nm2 <- unlist(nm2[ef])

  # Identify root of set_ parameters, remove set_ and add roots to nm1
  nm1_set <- nm1[substring(nm1,1,4)=="set_"]
  nm1_wo_set <- gsub("set_","",nm1_set)
  nm1 <- c(nm1,nm1_wo_set)


  # Identify root of set_ parameters, remove set_ and add roots to nm2
  nm2_set <- nm2[substring(nm2,1,4)=="set_"]
  nm2_wo_set <- gsub("set_","",nm2_set)
  nm2 <- c(nm2,nm2_wo_set)

  # For all nm1, replace them with nm2
  for(i in seq_along(nm1)){
    nm1_i <- nm1[i]
    nm2_i <- nm2[i]

    dat <- gsub(pattern=nm1_i,replacement=nm2_i,x=dat,fixed=TRUE)
    tpl <- gsub(pattern=nm1_i,replacement=nm2_i,x=tpl,fixed=TRUE)
    cxx <- gsub(pattern=nm1_i,replacement=nm2_i,x=cxx,fixed=TRUE)

  }


  return(list("dat"=dat,"tpl"=tpl,"cxx"=cxx))
}
