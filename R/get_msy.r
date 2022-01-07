#' get_msy
#'
#' Calculations of stock dependent quantities over a range of fishing mortalities (F). This is code from the Beaufort Assessment Model converted from ADMB to R.
#' @param rdat BAM output rdat file object
#' @param maxF Maximum fishing mortality rate. numeric vector 
#' @param nFopt length of F vector for grid search optimization
#' @param nFres length of F vector for filtering grid search results to limit output size
#' @param SPR_props proportions of SPR to compute reference points at. numeric vector
#' @param nages Number of ages in population model. numeric vector
#' @param steep Beverton-Holt steepness parameter. numeric vector
#' @param R0 Beverton-Holt R0 parameter (virgin recruitment). Numbers of fish at first age (often age-0 or age-1). numeric vector
#' @param BiasCor bias correction. numeric vector
#' @param SR_switch 1 for Beverton-Holt, 2 for Ricker, 3 for none (average recruitment). integer
#' @param M Natural mortality rate. numeric vector of length \code{nages}
#' @param sel_wgted_L effort-weighted, recent selectivities toward landings. numeric vector of length \code{nages}
#' @param sel_wgted_D effort-weighted, recent selectivities toward discards. numeric vector of length \code{nages}
#' @param wgt_klb Weight-at-age in thousand pounds. numeric vector of length \code{nages}
#' @param reprod Reproductive potential at-age vector. numeric vector of length \code{nages}
#' @param spawn_time_frac time of year of peak spawning, as a fraction of the year. numeric
#' @param wgt_mt Weight-at-age in metric tons for fish in the population. numeric vector of length \code{nages} 
#' @param wgt_klb Weight-at-age in thousand pounds for fish in the population. numeric vector of length \code{nages}
#' @param wgt_klb_L Weight-at-age in thousand pounds for fish in the landings. numeric vector of length \code{nages}
#' @param wgt_klb_D Weight-at-age in thousand pounds for fish in the discards. numeric vector of length \code{nages}
#' @param make_plot Indicate whether or not to plot results. logical
#' @param plot_digits number of significant digits to show in plot. numeric
#' 
#' @keywords bam stock assessment fisheries population dynamics
#' @author Erik Williams, Kyle Shertzer, and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' rdat <- rdat_VermilionSnapper
#' res <- get_msy(rdat,make_plot=TRUE)
#' }

get_msy <-  function(rdat,
                     maxF=NULL,
                     nFopt=NULL,
                     nFres=100,
                     SPR_props = c(0.3,0.4,0.5),
                     nages=NULL,
                     sel_wgted_L=NULL,
                     sel_wgted_D=NULL,
                     M=NULL,
                     reprod=NULL,
                     spawn_time_frac=NULL,
                     R0=NULL,
                     steep=NULL,
                     BiasCor=NULL,
                     SR_switch=NULL,
                     wgt_mt=NULL,
                     wgt_klb=NULL,
                     wgt_klb_L=NULL,
                     wgt_klb_D=NULL,
                     make_plot=FALSE, 
                     plot_digits=3
                     ){
  if(is.null(maxF)){
    maxF <- max(rdat$eq.series$F.eq)
  }
  if(is.null(nFopt)){
    nFopt <- nrow(rdat$eq.series)
  }
  
  if(is.null(nages)){
    nages <- length(rdat$a.series$age)
  }
  if(is.null(sel_wgted_L)){
    sel_wgted_L <- rdat$sel.age$sel.v.wgted.L
  }
  if(is.null(sel_wgted_D)){
    sel_wgted_D <- rdat$sel.age$sel.v.wgted.D
    if(is.null(sel_wgted_D)){
      sel_wgted_D <- sel_wgted_L*0
    }
  }
  if(is.null(M)){
    M <- rdat$a.series$M
  }
  if(is.null(spawn_time_frac)){
    spawn_time_frac <- rdat$parms$spawn.time
  }
  if(is.null(reprod)){
    reprod <- rdat$a.series$reprod
  }
  if(is.null(SR_switch)){  
    SR_switch <- c("BH-steep"="beverton_holt","Ricker-steep"="ricker","Mean"="none")[rdat$info$rec.model]
  }
  if(is.null(R0)){  
    R0 <- exp(rdat$parm.cons$log_R0[8])
  }
  if(is.null(steep)){
  steep <- rdat$parm.cons$steep[8]
  }
  if(is.null(BiasCor)){
    BiasCor <-  rdat$parms[grepl("biascorr",names(rdat$parms))][[1]]
  }
  if(is.null(wgt_mt)){
  wgt_mt <- rdat$a.series$wgt.mt
  }
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
  
  dzero <- 0.00001 
  
  # Abbreviate stuff
  ni <- nFopt
  na <- nages
  selL <- sel_wgted_L
  selD <- sel_wgted_D
  st <- spawn_time_frac
  
  # Stuff initialized early in the BAM tpl (e.g. PARAMETER_SECTION)
  # F_msy <- rep(NA,ni) # not needed
  Fvec <- seq(0,maxF,length=ni)
  N_am <- N_am_sp <- rep(NA, na)
  spr <- SPR <- R_eq <- SSB_eq <- B_eq <- L_eq_klb <- L_eq_knum <- D_eq_klb <- D_eq_knum <- rep(NA,nFopt)
  # iter_inc_msy <- maxF/(ni-1) # not needed
  
  # Usually computed by get_weighted_current()
    
    
  #fill in Fs for msy and per-recruit analyses # not needed
  # Fvec[1]=0.0
  # for (ff in 2:ni) {
  #   Fvec[ff] <- Fvec[ff-1] + iter_inc_msy
  #   }
  
  # compute values as functions of F
# for(ff=1; ff<=ni; ff++)
  # FOR EACH LEVEL OF F..
  for (ff in 1:ni) {
  
  #uses fishery-weighted F's
    Z_am <- 0.0
    F_L_am <- 0.0
    F_D_am <- 0.0
          
    F_L_am <- Fvec[ff]*selL
    F_D_am <- Fvec[ff]*selD    
    Z_am <- M + F_L_am + F_D_am         
    
    N_am[1] <- 1.0
    for (iage in 2:na){
      N_am[iage] <- N_am[iage-1]*exp(-Z_am[iage-1])
      }
    N_am[na] <- N_am[na]/(1.0-exp(-Z_am[na]))
    N_am_sp[1:(na-1)] <- N_am[1:(na-1)]*exp(-Z_am[1:(na-1)]*st)                 
    N_am_sp[na] <- (N_am_sp[na-1]*
                              (exp(-(Z_am[na-1]*(1.0-st) + 
                            Z_am[na]*st) )))/(1.0-exp(-Z_am[na]))
    
    # Spawners per recruit                 
    spr[ff] <- sum(N_am_sp*reprod)
    
    # Spawning potential ratio
    SPR[ff] <- spr[ff]/spr[1]
	        
    R_eq[ff] <- max(dzero,sr_eq(SR_switch=SR_switch, R0=R0, steep=steep, BiasCor=BiasCor, spr_F0=spr[1], spr_F=spr[ff]))
    
    N_am <-    N_am*R_eq[ff]
    N_am_sp <- N_am_sp*R_eq[ff]
    
    L_am <- N_am*(F_L_am/Z_am)*(1-exp(-Z_am))
    D_am <- N_am*(F_D_am/Z_am)*(1-exp(-Z_am))
    
    SSB_eq[ff]    <- sum(N_am_sp*reprod)
  	B_eq[ff]      <- sum(N_am*wgt_mt)
    L_eq_klb[ff]  <- sum(L_am*wgt_klb_L)
    L_eq_knum[ff] <- sum(L_am)/1000  
    D_eq_klb[ff]  <- sum(D_am*wgt_klb_D)   
    D_eq_knum[ff] <- sum(D_am)/1000   
  }  

  # Compute reference points
  msy_klb <- max(L_eq_klb)
  
  ff_msy <- which(L_eq_klb==msy_klb)
  
  ff_SPR_props <- sapply(SPR_props,function(x){which.min(abs(SPR-x))})
  
  # Build data frame of values at F
  eq_data <- data.frame("F"=Fvec,
                        L_klb=L_eq_klb, L_knum=L_eq_knum,
                        D_klb=D_eq_klb, D_knum=D_eq_knum,
                        B_mt=B_eq, SSB=SSB_eq, R=R_eq,
                        spr=spr, SPR=SPR
                        )
  
  # Identify msy reference points
  ref_pts_msy <- eq_data[ff_msy,]
  rownames(ref_pts_msy) <- "Fmsy"
  
  # Identify SPR reference points
  ref_pts_SPR <- eq_data[ff_SPR_props,]
  rownames(ref_pts_SPR) <- paste0("F",SPR_props*100)
  
  # Subset eq_data to return
  eq_data_res <- eq_data[floor(seq(1,nFopt,length=nFres)),]
  
  # Warnings
  if(rdat$parm.cons$steep[4]<0){
    warning(paste0("steepness was fixed at ",rdat$parm.cons$steep[8]," in this assessment"))
  }
  if(ff_msy==nFopt){
    warning("msy_klb was identified at maxF and probably does not represent the true maximum.")
  }
  if(nFopt%in%ff_SPR_props){
    warning("One of the SPR reference points was identified at maxF. Try a higher maxF to correctly identify this value.")
  }
  
  if(make_plot){
    plotNK <- function(...,rp=NA,rp_name=NULL){
      plot(type="l",lwd=2,...)
      if(!is.na(rp)){
        points(ref_pts_msy[,c("F",rp)],type="p",col="blue",pch=16)
        if(is.null(rp_name)){rp_name <- rp}
          rp_labels <- bquote(.(rp_name) == .(signif(ref_pts_msy[,rp],plot_digits)))
                  text(ref_pts_msy[,c("F",rp)],labels=rp_labels,pos=4)  
        }
        }
    par(mar=c(2,2,2,1),mfrow=c(3,2),mgp=c(1.1,0.1,0),tck=0.01)
    with(eq_data_res,{
    plotNK(F,spr)
      points(ref_pts_SPR[,c("F","spr")],type="p",col="blue",pch=16)
      sapply(seq_along(SPR_props),function(x){
        text(ref_pts_SPR[x,c("F","spr")],
             labels=bquote(spr[.(SPR_props[x]*100)] == .(signif(ref_pts_SPR[x,"spr"],plot_digits))),pos=4)})
    plotNK(F,SPR)
      points(ref_pts_SPR[,c("F","SPR")],type="p",col="blue",pch=16)
      sapply(seq_along(SPR_props),function(x){
        text(ref_pts_SPR[x,c("F","SPR")],
             labels=bquote(F[.(SPR_props[x]*100)] == .(signif(ref_pts_SPR[x,"F"],plot_digits))),pos=4)})
    plotNK(F,L_klb, rp="L_klb",rp_name = bquote(MSY))
    plotNK(F,R,rp="R",rp_name = bquote(R[MSY]))
    plotNK(F,SSB,rp="SSB",rp_name = bquote(SSB[MSY]))
    plotNK(F,B_mt,rp="B_mt",rp_name = bquote(B[MSY]))
    })
    }
  
  invisible(list(ref_pts_msy=ref_pts_msy, ref_pts_SPR=ref_pts_SPR, eq_data=eq_data_res))
}