#' get_weighted_current
#'
#' Calculate current weight-at-age for landings and discards. This is code from the Beaufort Assessment Model converted from ADMB to R.
#' @param rdat BAM output rdat file object
#' @param selpar_n_yrs_wgted Number years at end of time series over which to average sector Fs, for weighted selectivities. Almost always set to 3 in BAM. numeric vector
#' @param duplicate_sel Specify if selectivities should be duplicated for certain fleets. Otherwise the function assumes that each fleet with F values
#' has a selectivity defined with a fleet abbreviation matching abbreviation used in columns of rdat$parm.cons with "log_avg_F" in the name. Must be specified
#' as a list; e.g. duplicate_sel = data.frame("from"="HB_D","to"="GR_D"), specifies that the selectivity matrix indicated by the abbreviation HB_D should be copied with
#' a name using the abbreviation GR_D.
#' @param wgt_klb Weight-at-age in thousand pounds for fish in the population.
#' @param wgt_klb_L Weight-at-age in thousand pounds for fish in the landings.
#' @param wgt_klb_D Weight-at-age in thousand pounds for fish in the discards.
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' get_weighted_current(rdat = rdat_Tilefish)
#' # In some cases you need to duplicate the selectivity of one fleet to fill in for another, if that was done in the assessment
#' get_weighted_current(rdat = rdat_VermilionSnapper,  duplicate_sel = data.frame("from"="HB_D","to"="GR_D"))
#' }

# NOTES:
# mfexp() in ADMB is basically exp() in R
# elem_prod() in ADMB is "Element-wise multiplication of two vectors" which is what R does normally

# Abbreviated object names in function
# Fa <- log_avg_F
# Fd <- log_F_dev
# ey <- endyr

get_weighted_current <-  function(rdat,
                                  selpar_n_yrs_wgted=3,
                                  duplicate_sel=data.frame("from"=NULL,"to"=NULL),
                                  wgt_klb=NULL,
                                  wgt_klb_L=NULL,
                                  wgt_klb_D=NULL
){
  if(is.null(wgt_klb)){
    wgt_klb <- rdat$a.series$wgt.klb
  }
  
  if(is.null(wgt_klb_L)){
    L_klb_ix <- grep("wgt.wgted.L.klb",names(rdat$a.series))
    if(length(L_klb_ix)==1){
      wgt_klb_L <- rdat$a.series[,L_klb_ix]
    }else{
      wgt_klb_L <- wgt_klb
    }
  }
  
  if(is.null(wgt_klb_D)){
    D_klb_ix <- grep("wgt.wgted.D.klb",names(rdat$a.series))
    if(length(D_klb_ix)==1){
      wgt_klb_D <- rdat$a.series[,D_klb_ix]
    }else{
      wgt_klb_D <- wgt_klb
    }
  }  
  
  
  if(is.null(wholewgt_klb)){wholewgt_klb <- rdat$a.series$wgt.klb}
  
  prm <- rdat$parms
  pc <- rdat$parm.cons
  pts <- rdat$parm.tvec
  sa <- rdat$sel.age
  
  ey <- prm$endyr
  
  # List of selectivity matrices (year,age)
  sm <- local({
    a <- sa[grepl("sel.m",names(sa))]
    names(a) <- gsub("sel.m.","",names(a))
    a
  })
  
  # Duplicate specified selectivity matrices
  if(!is.null(duplicate_sel$from)&!is.null(duplicate_sel$to)){
    for(i in 1:nrow(duplicate_sel)){
      a <- duplicate_sel[i,]
      sm[[a$to]] <- sm[[a$from]]
    }
  }
  
  
  # Data frame of end year selectivities by fleet
  sey <- do.call(cbind,lapply(sm,function(x){x[paste(ey),]}))
  
  
  ## Abbreviate stuff
  ny <- selpar_n_yrs_wgted
 
  # log_avg_F (points)
  Fa <- local({
    a <- names(pc)
    b <- a[grepl("log_avg_F",a)]
    c <- pc[8,b]
    names(c) <- gsub("log_avg_F_","",names(c))
    names(c) <- gsub("L_","",names(c))
    c
  })
  
  # log_F_dev (time series vectors)
  Fd <- local({
    a <- names(pts)
    b <- a[grepl("log.F.dev.",a)]
    c <- pts[,b]
    names(c) <- gsub("log.F.dev.","",names(c))
    names(c) <- gsub(".","_",names(c),fixed=TRUE)
    c    
  })
  
  # Sum of recent log_F_dev (i.e. sum(log_F_dev_fleet((endyr-selpar_n_yrs_wgted+1),endyr))) 
  Fd_ey_sum <- colSums(Fd[paste(as.numeric(ey)-2:0),])

  
  Ffleet <-     local({
    a <- unlist(
      sapply(names(Fa),
             function(x){
               exp((Fa[x]*ny + Fd_ey_sum[x])/ny)
             }
      )
    )
    names(a) <- names(Fa)
    a
  })
  
  # F_temp_sum
  Ffleet_sum <- sum(Ffleet,
    na.rm=TRUE) # Exclude fleets with no landings at the end of the assessment period
  
  # F_prop
  Fprop <-     Ffleet/Ffleet_sum
  
  # log_F_dev_end
  Fd_end <- Fd_ey_sum/ny
  
  F_end_all <- local({
    a <- unlist(
      sapply(names(Fa),
             function(x){
               sey[,x] *  as.numeric(exp(Fa[x] + Fd_end[x]))
             }
      )
    )
    a <- as.data.frame(a)
    names(a) <- names(Fa)
    a
  })
  
  F_end_L <- rowSums(F_end_all[!grepl("_D",names(F_end_all))],na.rm=TRUE)
  F_end_D <- rowSums(F_end_all[grepl("_D",names(F_end_all))],na.rm=TRUE)
  
  F_end <- F_end_L+F_end_D
  F_end_apex <- max(F_end)
  
  sel_wgted_tot <- F_end/F_end_apex
  sel_wgted_L <- sel_wgted_tot * (F_end_L/F_end)
  sel_wgted_D <- sel_wgted_tot * (F_end_D/F_end)
  
  wgt_wgted_L_denom <- sum(Fprop[!grepl("_D",names(Fprop))],na.rm=TRUE) # Sum of Fprop for landings
  wgt_wgted_D_denom <- sum(Fprop[grepl("_D",names(Fprop))],na.rm=TRUE)  # Sum of Fprop for discards
  
  wgt_wgted_L_klb <- local({
    L_abb <- names(Fa)[!grepl("_D",names(Fa))]
    a <- unlist(
      sapply(L_abb,
             function(x){
               Fprop[x] /  wgt_wgted_L_denom*wgt_klb_L
             }
      )
    )
    a <- as.data.frame(a)
    names(a) <- L_abb
    rowSums(a,na.rm=TRUE)
  })
  
  wgt_wgted_D_klb <- local({
    D_abb <- names(Fa)[grepl("_D",names(Fa))]
    a <- unlist(
      sapply(D_abb,
             function(x){
               Fprop[x] /  wgt_wgted_D_denom*wgt_klb_D
             }
      )
    )
    a <- as.data.frame(a)
    names(a) <- D_abb
    rowSums(a,na.rm=TRUE)
  })
  
  return(list(F_end_L=F_end_L,F_end_D=F_end_D, F_end=F_end, F_end_apex=F_end_apex,
              sel_wgted_tot=sel_wgted_tot, sel_wgted_L=sel_wgted_L, sel_wgted_D=sel_wgted_D,
              wgt_wgted_L_klb=wgt_wgted_L_klb,wgt_wgted_D_klb=wgt_wgted_D_klb))
  
}