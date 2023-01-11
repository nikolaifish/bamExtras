#' run_proj
#'
#' Deterministic projections
#' @param rdat BAM output rdat (list) object read in with dget().
#' @param run_bam_args list of arguments passed internally to run_bam if any are provided
#' @param bam2r_args list of arguments passed internally to bam2r if any are provided
#' @param pstar Value to use for applying pstar scaling of Fmsy. If set to NULL, it is not used and pstar code will not be run.
#' @param styr_proj start year (i.e. first calendar year) of projection
#' @param nyb number of years in calculations that include recent landings (L) or discards and their cvs, index cvs (U), recruitment (R), or numbers of samples in age or length compositions (comp). Fleets without L or D within the last nyb years of the base model will not be projected forward. named list
#' @param nyc  number of years to maintain current conditions before implementing management
#' @param nyp  number of years of projection. Must be >nyc
#' @param F_cur current fully selected fishing mortality rate
#' @param F_proj fully selected fishing mortality rate
#' @param L_cur current level of landings
#' @param ages   ages. numeric vector
#' @param N_styr_proj abundance at age in styr_proj. numeric vector
#' @param S_styr_proj spawning stock size at age in styr_proj. Often biomass in mt, sometimes eggs in n. numeric vector
#' @param M      natural mortality rate, at age
#' @param wgt    weight at age of population. numeric vector
#' @param wgt_L  weight at age of landings. numeric vector
#' @param wgt_D  weight at age of discards. numeric vector
#' @param wgt_F_flt list of weight vectors (klb by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param len_F_flt list of length vectors (mm by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param reprod reproductive contribution at age to SSB. numeric vector
#' @param sel_L  selectivity at age to compute landings. numeric vector
#' @param sel_D  selectivity at age to compute dead discards. numeric vector
#' @param sel_tot  selectivity at age to compute Z (includes landings and discards). numeric vector
#' @param sel_F_flt list of selectivity vectors (by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param spawn_time time of year for peak spawning
#' @param SR_par list of parameters associated with the stock-recruit (SR) relationship. Currently only works with a Beverton-Holt SR relationship for which the parameters are: h = steepness of spawner-recruit function, R0 = virgin recruitment of spawner-recruit function, Phi0 = virgin spawners per recruit, biascorr = bias correction
#' @param SR_method Spawner-recruit function BH = Beverton-Holt, GM = constant geometric mean of nyb$R recent years of t.series$recruits. By default, run_proj will try to use BH, but if it can't find all the parameters, it will resort to GM (and will return a message letting you know that)
#' @param age_error age error matrix to use for the projection years. By default, the function uses the same age_error matrix from the base model
#' @param project_bam logical. After the projection is run, should the function build a new set of projected bam files?
#' @param plot logical. Produce plots of extended base model including projected data
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer, Erik Williams, and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' }

run_proj <- function(rdat = NULL,
                     run_bam_args = NULL,
                     bam2r_args = NULL,
                     pstar = NULL,
                     styr_proj = NULL,
                     nyb = list(L=3, U=3, R=5, comp=3),
                     nyc = 3,
                     nyp = 5,
                     F_cur = NULL,
                     F_proj = NULL,
                     L_cur = NULL,
                     ages = NULL,
                     N_styr_proj = NULL,
                     S_styr_proj = NULL,
                     sel_L = NULL,
                     sel_D=NULL,
                     sel_tot = NULL,
                     sel_F_flt = NULL,
                     wgt = NULL,
                     wgt_L = NULL,
                     wgt_D = NULL,
                     wgt_F_flt = NULL,
                     len_F_flt = NULL,
                     reprod = NULL,
                     M = NULL,
                     spawn_time = NULL,
                     SR_par = NULL,
                     SR_method = "BH",
                     age_error = NULL,
                     project_bam = FALSE,
                     plot = FALSE

) {
library(ggplot2)
library(tidyr)

if(!is.null(run_bam_args)){
  if(!"return_obj"%in%names(run_bam_args)){
    run_bam_args <- c(run_bam_args,list(return_obj=c("dat","tpl","cxx","rdat")))
  }
  run_bam_out <- do.call(run_bam,run_bam_args)
  rdat <- run_bam_out$rdat
}

  # if(!is.null(rdat)){
    styr <- rdat$parms$styr
    endyr <- rdat$parms$endyr
    yb <- styr:endyr # years of the base model

    a.series <- rdat$a.series
    t.series <- rdat$t.series
    sel_age <- rdat$sel.age

    yrs_L_b <- paste((endyr-(nyb$L-1)):endyr) # years of recent landings (i.e. from the base model)
    yrs_R_b <- paste((endyr-(nyb$R-1)):endyr)
    R_b <- t.series[yrs_R_b,"recruits"]
    R_b_gm <- bamExtras::geomean2(R_b)

    # Identify F-at-age for each fleet in endyr
    # L = (F/Z)*N*(1-exp(-Z))
    # L/(N*(1-exp(-Z))) = F/Z
    # L*Z/(N*(1-exp(-Z))) = F
    # F = L*Z/(N*(1-exp(-Z)))
    Cn <- rdat$CLD.est.mats[names(rdat$CLD.est.mats)[grepl("^(Ln|Dn)((?!total).)*$",names(rdat$CLD.est.mats),perl=TRUE)]]
    names(Cn) <- gsub("^([LD])(n)(.*)","Cn.\\1\\3",names(Cn))
    Cw <- rdat$CLD.est.mats[names(rdat$CLD.est.mats)[grepl("^(Lw|Dw)((?!total).)*$",names(rdat$CLD.est.mats),perl=TRUE)]]
    names(Cw) <- gsub("^([LD])(w)(.*)","Cw.\\1\\3",names(Cw))
    C_ts_cv <- t.series[paste(yb),grepl("^cv.[DL]",names(t.series))] # cvs associated with components of the catch (removals) during the base years
    Z <- rdat$Z.age[paste(yb),]
    N <- rdat$N.age[paste(yb),]
    Nmdyr <- rdat$N.age.mdyr[paste(yb),]
    R <- N[,1] # Recruits
    # Like F_fleet in bam tpl (e.g. F_cHL). Same computation as bam
    F_flt <- lapply(Cn,function(x){x*Z/(N*(1-exp(-Z)))}) # Computes a list of F (year,age) for each fleet, which are not in the rdat
    # names(F_flt) <- paste("F",gsub("([LD])(n)(.*)","\\1\\3",names(F_flt)),sep=".")
    names(F_flt) <- gsub("^Cn","F",names(F_flt))

    if(is.null(ages)){ages <- a.series$age}
    nages <- length(ages)
    len <- rdat$a.series$length
    if(is.null(age_error)){age_error <- rdat$age.error$error.mat
    if(is.null(age_error)){ # If there is no age error matrix in the rdat, just use an identity matrix
      age_error <- diag(length(ages))
      dimnames(age_error) <- list("age"=ages,"age"=ages)
      age_error
    }
    }
    if(is.null(spawn_time)){spawn_time <- rdat$parms$spawn.time}
    if(is.null(reprod)){reprod <- a.series$reprod}
    if(is.null(styr_proj)){styr_proj <- endyr+1}
    yp <- styr_proj:(styr_proj+nyp-1) # years of the projection period
    if(is.null(F_cur)){F_cur <- bamExtras::geomean2(t.series[yrs_L_b,"F.full"],na.rm=TRUE)}
    if(is.null(F_proj)){F_proj <- rdat$parms$Fmsy}
    if(is.null(L_cur)){
      is.total.L.klb <- "total.L.klb"%in%names(t.series)
      if(!is.total.L.klb){
        total.L.name <- names(t.series)[grepl("total.L",names(t.series))][1]
        message(paste("message: total.L.klb not found in t.series. L_cur is computed from",total.L.name,"instead\n"))
        total.L <- t.series[,total.L.name,drop=FALSE]
      }else{
        total.L <- t.series[,"total.L.klb",drop=FALSE]
        message(paste("message: L_cur is computed from t.series$total.L.klb\n"))

      }
      L_cur <- mean(total.L[yrs_L_b,])
      }


    if(is.null(sel_L)){sel_L <- sel_age$sel.v.wgted.L}

    if(is.null(sel_D)){
      if(!is.null(sel_age$sel.v.wgted.D)){
        sel_D <- sel_age$sel.v.wgted.D
      }else{
        message("message: sel_age$sel.v.wgted.D not found. Assessment may not model discards?\n")
      }
    }

    if(is.null(sel_tot)){sel_tot <- sel_age$sel.v.wgted.tot}

    if(is.null(sel_F_flt)){
      LD_abb <- gsub(".pr$","",names(t.series)[grepl("^[LD].*.pr$",names(t.series),perl=TRUE)]) # landings abbreviations
      LD_abb_ts <- paste0("F.",gsub("^(D.)(.*)","\\2.D",gsub("^L.","",LD_abb)))
      Fsum_flt <- t.series[paste(yb),LD_abb_ts] # F time series for each fleet
      names(Fsum_flt) <- paste0("F.",LD_abb)
      Fsum <- rowSums(Fsum_flt) # Should be equal to t.series$Fsum

      # Compute selectivity for each fleet associated with an F value.
      # NOTE: Computing selectivities is somewhat more reliable than trying to find
      # then in the rdat, in instances when selectivities from one fleet are used
      # for multiple fleets (e.g. when the headboat selectivity is used for the MRIP landings)
      sel_F_flt <- lapply(names(F_flt),function(x){
        a <- pmax(pmin(F_flt[[x]]/Fsum_flt[,x],1),0)
        a[is.na(a)] <- 0 # Replace NA with zero (NaN will occur when Fsum_flt values are zero, as in years when fleets do not have any landings)
        a
      })
      names(sel_F_flt) <- gsub("^F.","sel.",names(F_flt))

      F_ybgm_flt <- apply(tail(Fsum_flt,nyb$L),2,bamExtras::geomean2) # geomean F by fleet during the last nyb years of the base model
      F_prop_flt <- F_ybgm_flt/sum(F_ybgm_flt) # Proportion of F attributed to each fleet at the end of the assessment
      # This is just F_ybgm_flt multiplied by the selectivity in the endyr of the base model for each fleet
      # In the bam tpl, these vectors are named F_end although they are summed by landings (F_end_L), discards (F_end_D), or both (F_end)
      F_ybgm_flt_a <- as.data.frame(
        lapply(1:length(sel_F_flt),function(i){
          nm_x <- names(sel_F_flt)[i]
          sel_endyr_x <- sel_F_flt[[nm_x]][paste(endyr),]
          F_ybgm_flt_x <- F_ybgm_flt[[gsub("^sel.","F.",nm_x)]]
          F_ybgm_flt_x*sel_endyr_x
        }
        )
      )
      names(F_ybgm_flt_a) <- names(F_ybgm_flt)

      # Sum across fleets (these objects are named as in bam tpl)
      F_end_L <- rowSums(F_ybgm_flt_a[grepl("^F.L",names(F_ybgm_flt_a))])
      F_end_D <- rowSums(F_ybgm_flt_a[grepl("^F.D",names(F_ybgm_flt_a))])
      F_end <- rowSums(as.data.frame(F_ybgm_flt_a))
      F_end_apex <- max(F_end)
      sel_wgted_tot <- F_end/F_end_apex # Should equal to rdat$sel.age$sel.v.wgted.tot
      sel_wgted_L <-   sel_wgted_tot*(F_end_L/F_end) # Should equal rdat$sel.age$sel.v.wgted.L # F_end_L/F_end_apex # as in tpl
      sel_wgted_D <-   sel_wgted_tot*(F_end_D/F_end) # Should equal rdat$sel.age$sel.v.wgted.D # F_end_D/F_end_apex # as in tpl
      sel_wgted_F_flt <- F_ybgm_flt_a/F_end_apex     # by extension, scale selectivity at age for each fleet
      # Note that: round((sel_wgted_L+sel_wgted_D) == rowSums(sel_wgted_F_flt)
    }

    ## Weights of fish (year, age)
    if(is.null(wgt)){wgt <- a.series$wgt.mt}
    if(is.null(wgt_L)){
      if(!is.null(a.series$wholewgt.wgted.L.klb)){
        wgt_L <- a.series$wholewgt.wgted.L.klb
        message("message: wgt_L set equal to a.series$wholewgt.wgted.L.klb\n")
      }else if(!is.null(a.series$gutwgt.wgted.L.klb)){
        wgt_L <- a.series$gutwgt.wgted.L.klb
        message("message: wgt_L set equal to a.series$gutwgt.wgted.L.klb\n")
      }
    }
    if(is.null(wgt_D)){
      if(!is.null(a.series$wholewgt.wgted.D.klb)){
        wgt_D <- a.series$wholewgt.wgted.D.klb
      }else{
        message("message: a.series$wholewgt.wgted.D.klb not found. Assessment may not model discards.\n")
      }
    }

    ## Weights of fish (year, age) by fleet.
    if(is.null(wgt_F_flt)){
      wgt_F_flt <- rdat$size.age.fishery[grepl("^[a-zA_Z].*wgt",names(rdat$size.age.fishery))]
      # Reformat names to wgt.fleetType.fleet
      # (e.g. wgt.L.cHL for weight of landings in the commercial hook and line fleet)
      wgt_F_flt_abb <- local({
        a <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\3",names(wgt_F_flt),perl=TRUE)
        gsub("^([^D.])([a-zA-Z]+)","L.\\1\\2",gsub("(.*?)(?<=.)(.D)$","D.\\1",a,perl=TRUE))
      })

      wgt_F_flt_units <- local({
        a <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\1.\\5",names(wgt_F_flt),perl=TRUE)
        # a <- gsub("wgt","W",a)
        # a <- gsub("whole","W",a)
        # a <- gsub("gut","G",a)
        return(a)
      })
      names(wgt_F_flt) <- names(wgt_F_flt_units) <- paste("wgt", #wgt_F_flt_units,
                                                          wgt_F_flt_abb,sep=".")

      # Convert weights in lb to klb
      for(nm_i in names(wgt_F_flt)){
        xi <- wgt_F_flt[[nm_i]]
        if(grepl(".lb$",wgt_F_flt_units[[nm_i]])){
          wgt_F_flt[[nm_i]] <- xi/1000
          wgt_F_flt_units[[nm_i]] <- gsub(".lb$",".klb",wgt_F_flt_units[[nm_i]])
        }
      }
    }

    # Lengths of fish (year, age) by fleet.
    if(is.null(len_F_flt)){
      len_F_flt <- rdat$size.age.fishery[grepl("^[a-zA_Z]*.*len",names(rdat$size.age.fishery))]
      # Reformat names to len.fleetType.fleet
      # (e.g. len.L.cHL for length of landings in the commercial hook and line fleet)
      len_F_flt_abb <- local({
        a <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\3",names(len_F_flt),perl=TRUE)
        gsub("^([^D.])([a-zA-Z]+)","L.\\1\\2",gsub("(.*?)(?<=.)(.D)$","D.\\1",a,perl=TRUE))
      })

      len_F_flt_units <- local({
        a <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\1.\\5",names(len_F_flt),perl=TRUE)
        # a <- gsub("len","W",a)
        # a <- gsub("whole","W",a)
        # a <- gsub("gut","G",a)
        return(a)
      })
      names(len_F_flt) <- names(len_F_flt_units) <- paste("len", #len_F_flt_units,
                                                          len_F_flt_abb,sep=".")
    }


    # Recompute landings by fleet to compare with calculations used in projections
    Cn2 <- lapply(Cn,function(x){x*NA})
    Cw2 <- lapply(Cw,function(x){x*NA})
    for(i in 1:length(yb)){
    for(j in 1:ncol(Fsum_flt)){
      Cn2[[j]][i,] <- L_calc(F_flt[[j]][i,], Z[i,], N[i,])
      Cw2[[j]][i,] <- L_calc(F_flt[[j]][i,], Z[i,], N[i,], wgt_F_flt[[j]][i])
    }
    }


    # CPUE
    #   From bam tpl for SEDAR53:
    #     N_cHL(iyear)=elem_prod(elem_prod(Nmdyr(iyear),sel_cHL(iyear)),wholewgt_cHL_klb(iyear));
    #     pred_cHL_cpue(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
    U_a <- list()
    NU <- list()
    wgt_U <- list()
    U_ts <- t.series[paste(yb),grepl("^U.*pr$",names(t.series)),drop=FALSE]
    U_ts_cv <- t.series[paste(yb),grepl("^cv.U",names(t.series)),drop=FALSE] # cvs associated with U during the base years
    # names(U_ts) <- gsub("^(U.)(.*)(.)(pr)$","\\1\\4.\\2",names(U_ts))
    U_abb <- gsub("^U.|.pr$","",names(U_ts))
    names(U_ts) <- U_abb
    q_mn <- t.series[paste(yb),paste0("q.",U_abb),drop=FALSE]  # This is really a kind of mean q, since q_rate and q_DD_mult might scale it to compute U
    q_DD_mult <- t.series[paste(yb),"q.DD.mult",drop=FALSE] # density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
    q_rate <- q_mn*0+1
    q <- q_mn*NA # Initialize values for final q value (q_mn * q_rate * q_DD_mult) multiplied by N to compute CPUE

    sel_U <- sel_age[which(gsub("sel.[vm].|(wgted.)","",names(sel_age))%in%U_abb)]
    names(sel_U) <- U_abb

    for(i in 1:length(sel_U)){
      sel_U_i <- sel_U[[i]]
      sel_nm_i <- names(sel_U)[i]
      sel_abb_i <- gsub("sel.[vm].","",sel_nm_i)
      q_nm_i <- paste0("q.",sel_abb_i)
      U_type_i <- gsub("^wgt.(L.)?","",sel_abb_i) # 2023-1-10 I'm not sure this gsub is necessary anymore. I think it was coded when the names were different.
      wgt_U_i <- local({ # Set default wgt = 1
        a <- sel_U_i
        a*0+1
      })
      # If the U type is commercial, compute N_fleet in weight
      if(gsub("^(.{1})(.*)","\\1",U_type_i)=="c"){
        wgt_U_nm_i <- paste0("wgt.L.",sel_abb_i)
        if(wgt_U_nm_i%in%names(wgt_F_flt)){
          wgt_U_i <- wgt_F_flt[[paste0("wgt.L.",sel_abb_i)]]
        }else{
          message(paste("message: I could not find",wgt_U_nm_i, "to multiply by", sel_nm_i,"when computing", paste0("N_",sel_abb_i,"."), paste0("U_pr_",sel_abb_i), "will be computed in numbers instead of weight.\n"))
        }
      }
      NU_i <- Nmdyr*sel_U_i*wgt_U_i
      NU[[i]] <- NU_i
      q_mn_i <- q_mn[paste(yb),q_nm_i,drop=FALSE]
      q_rate_nm_i <- paste0("q.",sel_abb_i,".rate.mult")
      if(q_rate_nm_i%in%names(t.series)){
        q_rate_i <- t.series[paste(yb),q_rate_nm_i,drop=FALSE]
      }else{
        q_rate_i <- q_mn_i*0+1
        message(paste("message:",q_rate_nm_i,"not found in t.series. Setting q rate multiplier values to 1.\n"))
      }
      q_rate[,q_nm_i] <- q_rate_i

      q_i <- q_mn_i*q_rate_i*q_DD_mult
      q[,q_nm_i] <- q_i
      U_a_i <- unlist(q_i)*NU_i
      U_a[[i]] <- U_a_i
      wgt_U[[i]] <- wgt_U_i
      # U_pr_i <- rowSums(U_a_i)
    }
    names(U_a) <- names(wgt_U) <- names(NU) <- U_abb

    U <- as.data.frame(lapply(U_a,rowSums))

    ## age compositions from the base years
    # observed
    acomp_b_ob <- rdat$comp.mats[grepl("^acomp.*.ob$",names(rdat$comp.mats))]
    names(acomp_b_ob) <- local({
      a <- gsub("^(acomp.)(.*)(.ob)$","\\2",names(acomp_b_ob))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    # predicted
    acomp_b_pr <- rdat$comp.mats[grepl("^acomp.*.pr$",names(rdat$comp.mats))]
    names(acomp_b_pr) <- local({
      a <- gsub("^(acomp.)(.*)(.pr)$","\\2",names(acomp_b_pr))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })

    ## length compositions from the base years
    # observed
    lcomp_b_ob <- rdat$comp.mats[grepl("^lcomp.*.ob$",names(rdat$comp.mats))]
    names(lcomp_b_ob) <- local({
      a <- gsub("^(lcomp.)(.*)(.ob)$","\\2",names(lcomp_b_ob))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    # predicted
    lcomp_b_pr <- rdat$comp.mats[grepl("^lcomp.*.pr$",names(rdat$comp.mats))]
    names(lcomp_b_pr) <- local({
      a <- gsub("^(lcomp.)(.*)(.pr)$","\\2",names(lcomp_b_pr))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })

    lenprob <- list()

    ## Build age-length conversion matrix associated with each set of length comps
    ## (will often be the same for all comps)
    for(i in 1:length(lcomp_b_pr)){
      avail_len_cv_val <- names(rdat$parm.cons)[grepl("^len.cv",names(rdat$parm.cons))]
      if(length(avail_len_cv_val)>1){
        message(paste("message: found multiple len.cv.val in rdat:",paste(avail_len_cv_val,collapse=", "),"Currently only using len.cv.val"))
      }

      len_cv_val_i <- rdat$parm.cons$len.cv.val[8]
      len_sd_i <- len*len_cv_val_i

      lc_i <- lcomp_b_pr[[i]]
      lc_bm_i <- as.numeric(colnames(lc_i))  # bin mid point
      lc_bw_i <- median(diff(lc_bm_i)) # bin width
      lc_bl_i <- lc_bm_i-lc_bw_i/2  # bin lo (minimum)
      lenprob_i <- local({
        mat_i <- matrix(NA,nrow=length(ages),ncol=length(lc_bm_i),dimnames = list("age"=ages,"lenbin"=lc_bm_i))
        mat_i[,1] <- pnorm((lc_bl_i[2]-len)/len_sd_i)
        for(i in 2:(ncol(mat_i)-1)){
          mat_i[,i] <- pnorm((lc_bl_i[i+1]-len)/len_sd_i)-pnorm((lc_bl_i[i]-len)/len_sd_i)
        }
        mat_i[,ncol(mat_i)] <- 1-rowSums(mat_i[,1:(ncol(mat_i)-1)])
        mat_i
      })
      lenprob[[i]] <- lenprob_i
    }

    names(lenprob) <- names(lcomp_b_pr)


    ## Recompute predicted age and length comps
    # Note: These should be the same or very similar to age and length comps computed by BAM.
    #  But these are recomputed for the base years as a check since the computation
    #  will be used in the projection years
    # acomp_b_pr_2 <-  list()
    # for(nm_i in names(acomp_b_pr)){
    #   yrs_ac_i <- rownames(acomp_b_pr[[nm_i]])
    #   if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
    #     n_i <- Cn[[paste0("Cn.",nm_i)]][yrs_ac_i,]
    #   }else{
    #     n_i <- NU[[nm_i]][yrs_ac_i,] # If not (i.e. it's a fishery independent survey) use numbers associated with the index
    #   }
    #   P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    #   acomp_b_pr_2_i <- t(apply(P_n_i,1,function(x){colSums(x*age_error)}))
    #   acomp_b_pr_2[[nm_i]] <- acomp_b_pr_2_i
    # }

    # lcomp_b_pr_2 <-  list()
    # for(nm_i in names(lcomp_b_pr)){
    #   yrs_lc_i <- rownames(lcomp_b_pr[[nm_i]])
    #   lenprob_i <- lenprob[[nm_i]]
    #   if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
    #     n_i <- Cn[[paste0("Cn.",nm_i)]][yrs_lc_i,]
    #   }else{
    #     n_i <- NU[[nm_i]][yrs_lc_i,] # If not (i.e. it's a fishery independent survey) use numbers associated with the index
    #   }
    #   P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    #   lcomp_b_pr_2_i <- t(apply(P_n_i,1,function(x){colSums(x*lenprob_i)}))
    #   lcomp_b_pr_2[[nm_i]] <- lcomp_2_i
    # }

    ##
    # number of trips in compositions
    ntrip <- t.series[paste(yb),grepl("^([la]comp).*(.n)$",names(t.series))]
    ntrip[ntrip==-99999] <- NA
    names(ntrip) <- local({
      a <- gsub("(.n)$","",names(ntrip))
      b <- gsub("^([al]comp.)(.*)(.)(D)$","\\1D.\\2",a) # Convert .D suffix to D. prefix
      c <- gsub("^([al]comp.)([cr])(.*)","\\1L.\\2\\3",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    # number of fish in compositions
    nfish <- t.series[paste(yb),grepl("^([la]comp).*(.nfish)$",names(t.series))]
    nfish[nfish==-99999] <- NA
    names(nfish) <- local({
      a <- gsub("(.nfish)$","",names(nfish))
      b <- gsub("^([al]comp.)(.*)(.)(D)$","\\1D.\\2",a) # Convert .D suffix to D. prefix
      c <- gsub("^([al]comp.)([cr])(.*)","\\1L.\\2\\3",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })

    #Initial conditions
    #Recruits in yr 1 constrained to S-R curve, or else varies in stochastic runs
    if(is.null(N_styr_proj)){N_styr_proj <- tail(rdat$N.age,1)}
    N_endyr <- rdat$N.age[paste(endyr),]
    if(is.null(S_styr_proj)){
      S_styr_proj <- local({
        Z_styr_proj <- rdat$Z.age[paste(endyr),]
        N_endyr_spn <- N_endyr*exp(-1.0*Z_styr_proj*spawn_time)
        sum(N_endyr_spn*reprod)
      })
    }

    if(SR_method=="BH"){
      if(is.null(SR_par)){
        h <-        rdat$parms$BH.steep
        R0 <-       rdat$parms$BH.R0
        Phi0 <-     rdat$parms$BH.Phi0
        biascorr <- rdat$parms$BH.biascorr
      }else{
        for(i in 1:length(SR_par)){xi <- SR_par[i]; assign(x=names(xi),value=xi)}
      }
      SR_par_null <- c("h","R0","Phi0","biascorr")[c(is.null(h),is.null(R0),is.null(Phi0),is.null(biascorr))]
      # If not all of the BH SR parameters are found, resort to using the GM method
      if(length(SR_par_null)>0){
        message(paste("message: BH SR parameters not found:",paste(SR_par_null,collapse=", ")," Using SR_method = GM.\n"))
        SR_method <- "GM"
      }else{
        message("message: all BH SR parameters found. Using SR_method = BH.\n")
      }
    }else{
      message(paste("message: Using SR_method = GM.\n"))
    }

    if(is.null(M)){
      M <- a.series$M
    }

  # } # end if(!is.null(rdat))


  ###### Setup projections ######
    if(nyp>0){ # Only run this if there is actually a projection
  ## Build empty objects
  # generic empty objects
    ema <-  setNames(rep(NA,nages),ages) # age vector; a = ages
    emyp <- setNames(rep(NA,nyp),yp)     # year vector; yp = years of the projection period
    emypa <- matrix(NA,nrow=nyp,ncol=nages,dimnames=list("year"=yp,"age"=ages)) # matrix (year,age) during the projection period
    emuyp <- matrix(NA,nrow=nyp,ncol=length(U_a),dimnames=list("year"=yp,"U"=names(U_a))) # matrix (year, U) during the projection period

    emfypa <- local({
      a <- lapply(1:length(Fsum_flt),function(x){emypa})
      names(a) <- names(Fsum_flt)
      a
    })
    emuypa <- local({
      a <- lapply(1:length(U_a),function(x){emypa})
      names(a) <- names(U_a)
      a
    })

  ## Old objects
  # Note: these _a objects will be phase out once I build the new ones
  Z_p <-      emypa  # total mortality by age
  Nspwn_p <- emypa  # numbers at age at spawn_time
  Nmdyr_p <- emypa  # numbers at age at midyear
  N_p <-      emypa  # numbers at age by year
  NU_p <-    emuypa # numbers (or biomass) at age for each index of abundance, used in U_a calculations

  S_p <-   emyp # spawning stock (often biomass in mt, sometimes eggs in n)
  B_p <-   emyp # population biomass (mt)
  R_p <-   emyp # recruits (n)
  Fsum_p <-   emyp # sum F across all fleets (/yr)

  Ln_p <- emyp	# total landings (n)
  Lw_p <- emyp	# total landings (weight)
  Dn_p <- emyp	# total dead discards (n)
  Dw_p <- emyp	# total dead discards (weight)

  ## New objects
  # sel_F_flt_p <- lapply(1:length(Fsum_flt),function(x){emypa})
  # names(sel_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  U_a_p <-  emuypa    # predicted cpue at age by fleet
  wgt_U_p <- emuypa # weights associated with cpue at age by fleet
  q_p <- emuyp      # catchability (q)

  F_flt_p <- emfypa

  wgt_F_flt_p <- emfypa
  names(wgt_F_flt) <- gsub("^F","wgt",names(emfypa))

  len_F_flt_p <- emfypa
  names(len_F_flt) <- gsub("^F","len",names(emfypa))

  Cn_p <- emfypa # Catch in numbers by fleet (landings or discards)
  names(Cn_p) <- gsub("^F","Cn",names(emfypa))

  Cw_p <- emfypa # Catch in weight by fleet (landings or discards)
  names(Cw_p) <- gsub("^F","Cw",names(emfypa))

  ## Initialization
  N_p[1,] <- N_styr_proj
  S_p[1] <- S_styr_proj
  B_p[1] <- sum(N_p[1,]*wgt)

  Fsum_p[1:nyc] <- F_cur
  Fsum_p[(nyc+1):nyp] <- F_proj

  # Compute total Z and F at-age for each year of the projection
  # Z_p <-   outer(Fsum_p,M+sel_tot,FUN="*") # Z for population
  Z_p <-   t(M+t(outer(Fsum_p,sel_tot,FUN="*")))   # Z for population
  F_L_p <-       outer(Fsum_p,sel_L,  FUN="*")     # F for landings
  F_D_p <-       outer(Fsum_p,sel_D,  FUN="*")     # F for discards

  ## sel during projection years (by py, age, fleet)
  # for any source of F (landings or discards)
  # Fill with weighted selectivity from last year of base (endyr)
  sel_F_flt_p <- lapply(1:length(sel_F_flt),function(i){
    nm_i <- names(sel_F_flt)[i]
    sel_endyr_i <- sel_wgted_F_flt[,gsub("^sel.","F.",nm_i)]
    # sel_endyr_i <- sel_F_flt[[nm_i]][paste(endyr),] # Gets unweighted selectivities
    matrix(sel_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  # sel during projection years (by py, age, fleet) for any source of indices of abundance (cpue)
  # Fill with selectivity from last year of base (endyr)
  sel_U_p <- lapply(1:length(sel_U),function(i){
    nm_i <- names(sel_U)[i]
    sel_endyr_i <- sel_U[[nm_i]][paste(endyr),]
    matrix(sel_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_U_p) <- names(sel_U)

  # F during projection years (by py, age, fleet) for any source of F (landings or discards)
  # (analogous to F_L_p or F_D_p but by fleet)
  # Note that:
  #      F.L.tmp <- Reduce("+",F_flt_p[grepl("^F.L.",names(F_flt_p))])
  #      F.L.tmp == F_L_p
  F_flt_p <- lapply(1:length(F_flt_p),function(i){
    sel_F_flt_p[[i]]*Fsum_p
  })
  names(F_flt_p) <- names(Fsum_flt)

  # wgt during projection years (by py, age, fleet) for any source of F (landings or discards)
  # Fill with fish weights from last year of base (endyr)
  wgt_F_flt_p <- lapply(1:length(F_flt_p),function(i){
    nm_i <- names(F_flt_p)[i]
    wgt_i <- wgt_F_flt[[gsub("^F.","wgt.",nm_i)]]
    wgt_endyr_i <- wgt_i[paste(endyr),]
    matrix(wgt_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(wgt_F_flt_p) <- gsub("^F.","wgt.",names(F_flt_p))

  # len during projection years (by py, age, fleet) for any source of F (landings or discards)
  # Fill with fish weights from last year of base (endyr)
  len_F_flt_p <- lapply(1:length(F_flt_p),function(i){
    nm_i <- names(F_flt_p)[i]
    len_i <- len_F_flt[[gsub("^F.","len.",nm_i)]]
    len_endyr_i <- len_i[paste(endyr),]
    matrix(len_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(len_F_flt_p) <- gsub("^F.","len.",names(F_flt_p))

  # weight associated with cpue during the projection years
  # (mostly used to convert commercial cpue to weight. Otherwise set to 1 which leaves cpue in numbers)
  wgt_U_p <- lapply(1:length(wgt_U),function(i){
    nm_i <- names(wgt_U)[i]
    x_i <- wgt_U[[nm_i]]
    x_endyr_i <- x_i[paste(endyr),]
    matrix(x_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(wgt_U_p) <- names(wgt_U)

  # catchability associated with cpue during the projection years
  q_p[paste(yp),] <- rep(unlist(q[paste(endyr),]),each=nyp)

  # compositions
  acomp_p <- list()
  for(i in names(acomp_b_pr)){
  ac_i <- acomp_b_pr[[i]]
  acomp_p[[i]] <- matrix(NA,nrow=nyp,ncol=ncol(ac_i),dimnames=list("year"=yp,"age"=colnames(ac_i)))
  }

  lcomp_p <- list()
  for(i in names(lcomp_b_pr)){
    lc_i <- lcomp_b_pr[[i]]
    lcomp_p[[i]] <- matrix(NA,nrow=nyp,ncol=ncol(lc_i),dimnames=list("year"=yp,"lenbin"=colnames(lc_i)))
  }

  ntrip_p <- local({
    a <- ceiling(apply(tail(ntrip,nyb$comp),2,function(x){bamExtras::geomean2(x)}))
    matrix(a,nrow=nyp,ncol=length(a),dimnames=list(year=paste(yp),fleet=names(a)),byrow=TRUE)
  })

  nfish_p <- local({
    a <- ceiling(apply(tail(nfish,nyb$comp),2,function(x){bamExtras::geomean2(x)}))
    matrix(a,nrow=nyp,ncol=length(a),dimnames=list(year=paste(yp),fleet=names(a)),byrow=TRUE)
  })

  #### Projection loop
  ### years 1 to nyp-1
  if(nyp>1){
    for (i in 1:(nyp-1)) {
      yp_i <- paste(yp[i])

      ## Population (year i)
      B_p[i] <- sum(N_p[i,]*wgt)
      Nspwn_p[i,] <- N_p[i,]*(exp(-1.0*Z_p[i,]*spawn_time))
      Nmdyr_p[i,] <- N_p[i,]*(exp(-1.0*Z_p[i,]*0.5))
      S_p[i] <- sum(Nspwn_p[i,]*reprod)

      ## Population (year i+1)
      if(SR_method=="BH"){ # Beverton-Holt stock-recruit relationship
        N_p[i+1,1] <- biascorr*(0.8*R0*h*S_p[i])/(0.2*R0*Phi0*(1.0-h)+ (h-0.2)*S_p[i])
      }
      if(SR_method=="GM"){ # Geometric mean recruitment
        N_p[i+1,1] <- R_b_gm
      }

      N_p[(i+1),2:nages] <- N_p[i,1:(nages-1)]*(exp(-1.0*Z_p[i,1:(nages-1)]))
      N_p[(i+1),nages]   <- N_p[(i+1),nages] +
        N_p[i,nages]*(exp(-1.0*Z_p[i,nages])) #plus group
      ## cpue
      # By fleet
      for(j in 1:length(U_a_p)){
        nm_j <- names(U_a_p)[j]
        NU_p_ji <- Nmdyr_p[yp_i,]*sel_U_p[[nm_j]][yp_i,]*wgt_U_p[[nm_j]][yp_i,]
        NU_p[[nm_j]][yp_i,] <- NU_p_ji
        U_a_p[[nm_j]][yp_i,] <- NU_p_ji*q_p[yp_i,nm_j]

      }

      ## Fishery (year i)
      # Aggregate (Standard calculations)
      Ln_p[i] <- sum(L_calc(F_L_p[i,], Z_p[i,], N_p[i,]))
      Lw_p[i] <- sum(L_calc(F_L_p[i,], Z_p[i,], N_p[i,], wgt_L))
      if(any(!is.null(c(sel_D,wgt_D)))){
        Dn_p[i] <- sum(L_calc(F_D_p[i,], Z_p[i,], N_p[i,]))
        Dw_p[i] <- sum(L_calc(F_D_p[i,], Z_p[i,], N_p[i,], wgt_D))
      }
      # By fleet
      for(j in 1:ncol(Fsum_flt)){
        Cn_p[[j]][i,] <- L_calc(F_flt_p[[j]][i,], Z_p[i,], N_p[i,])
        Cw_p[[j]][i,] <- L_calc(F_flt_p[[j]][i,], Z_p[i,], N_p[i,], wgt_F_flt_p[[j]][i,])
      }

    }  # end i
  }  # end if(nyp>1)

  ### year nyp
  ## Population (year nyp)
  B_p[nyp] <- sum(N_p[nyp,]*wgt)
  Nspwn_p[nyp,] <- N_p[nyp,]*(exp(-1.0*Z_p[nyp,]*spawn_time))
  Nmdyr_p[nyp,] <- N_p[nyp,]*(exp(-1.0*Z_p[nyp,]*0.5))
  S_p[nyp] <- sum(Nspwn_p[nyp,]*reprod)
  R_p <- N_p[,1]

  ## cpue
  # By fleet
  for(j in 1:length(U_a_p)){
    nm_j <- names(U_a_p)[j]
    NU_p_ji <- Nmdyr_p[nyp,]*sel_U_p[[nm_j]][nyp,]*wgt_U_p[[nm_j]][nyp,]
    NU_p[[nm_j]][nyp,] <- NU_p_ji
    U_a_p[[nm_j]][nyp,] <- NU_p_ji*q_p[nyp,nm_j]

  }

  ## Fishery (year nyp)
  # Aggregate (Standard calculations)
  Ln_p[nyp] <- sum(L_calc(F_L_p[nyp,], Z_p[nyp,], N_p[nyp,]))
  Lw_p[nyp] <- sum(L_calc(F_L_p[nyp,], Z_p[nyp,], N_p[nyp,], wgt_L))
  if(any(!is.null(c(sel_D,wgt_D)))){
    Dn_p[nyp] <- sum(L_calc(F_D_p[nyp,], Z_p[nyp,], N_p[nyp,]))
    Dw_p[nyp] <- sum(L_calc(F_D_p[nyp,], Z_p[nyp,], N_p[nyp,], wgt_D))
  }

  # By fleet
  for(j in 1:ncol(Fsum_flt)){
    Cn_p[[j]][nyp,] <- L_calc(F_flt_p[[j]][nyp,], Z_p[nyp,], N_p[nyp,])
    Cw_p[[j]][nyp,] <- L_calc(F_flt_p[[j]][nyp,], Z_p[nyp,], N_p[nyp,], wgt_F_flt_p[[j]][nyp,])
  }

  U_p <- as.data.frame(lapply(U_a_p,rowSums))

  ## Compute predicted age and length comps during projection years
  for(nm_i in names(acomp_p)){
    yrs_acomp_i <- rownames(acomp_b_ob[[nm_i]])
    yrs_acomp_p_i <- rownames(acomp_p[[nm_i]])
    # If the last year of comps is within the last nyb$comp of the assessment, then project them.
    if(max(as.numeric(yrs_acomp_i))%in%tail(as.numeric(yb),nyb$comp)){
    if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
      n_i <- Cn_p[[paste0("Cn.",nm_i)]][yrs_acomp_p_i,,drop=FALSE]
    }else{
      n_i <- NU_p[[nm_i]][yrs_acomp_p_i,,drop=FALSE] # If not (i.e. it's a fishery independent survey) use numbers associated with the index
    }
    P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    acomp_p_i <- t(apply(P_n_i,1,function(x){colSums(x*age_error)}))
    acomp_p[[nm_i]] <- acomp_p_i
    }
  }

  for(nm_i in names(lcomp_p)){
    yrs_lcomp_i <- rownames(lcomp_b_ob[[nm_i]])
    yrs_lcomp_p_i <- rownames(lcomp_p[[nm_i]])
    # If the last year of comps is within the last nyb$comp of the assessment, then project them.
    if(max(as.numeric(yrs_lcomp_i))%in%tail(as.numeric(yb),nyb$comp)){
      lenprob_i <- lenprob[[nm_i]]
      if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
        n_i <- Cn_p[[paste0("Cn.",nm_i)]][yrs_lcomp_p_i,,drop=FALSE]
      }else{
        n_i <- NU_p[[nm_i]][yrs_lcomp_p_i,,drop=FALSE] # If not (i.e. it's a fishery independent survey) use numbers associated with the index
      }
      P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
      lcomp_p_i <- t(apply(P_n_i,1,function(x){colSums(x*lenprob_i)}))
      lcomp_p[[nm_i]] <- lcomp_p_i
    }
  }
  # Project cvs for catch and cpue
  C_ts_cv_p <- local({
    a <- matrix(apply(C_ts_cv,2,function(x){rep(bamExtras::geomean2(tail(x,nyb$L)),nyp)}),
                nrow=nyp,dimnames=list(paste(yp),names(C_ts_cv)))
    # rownames(a) <- paste(yp)
    a
  })

  U_ts_cv_p <- local({
    a <- matrix(apply(U_ts_cv,2,function(x){rep(bamExtras::geomean2(tail(x,nyb$U)),nyp)}),
                nrow=nyp,dimnames=list(paste(yp),names(U_ts_cv)))
    # rownames(a) <- paste(yp)
    a
  })

# Add projected values to base values
  N <- rbind(N,N_p)
  R <- c(R,R_p)
  U <- rbind(U,U_p)
  C_ts_cv <- rbind(C_ts_cv,C_ts_cv_p)
  U_ts_cv <- rbind(U_ts_cv,U_ts_cv_p)
  nfish <- rbind(nfish,nfish_p)
  ntrip <- rbind(ntrip,ntrip_p)

  lcomp <- local({
    a <- lapply(seq_along(names(lcomp_b_ob)),function(x){
      nm_x <- names(lcomp_b_ob)[x]
      ax <- rbind(lcomp_b_ob[[x]],lcomp_p[[x]])
      ax[complete.cases(ax),,drop=FALSE]
    })
    names(a) <- names(lcomp_b_ob)
    a
  })

  acomp <- local({
    a <- lapply(seq_along(names(acomp_b_ob)),function(x){
      nm_x <- names(acomp_b_ob)[x]
      ax <- rbind(acomp_b_ob[[x]],acomp_p[[x]])
      ax[complete.cases(ax),,drop=FALSE]
    })
    names(a) <- names(acomp_b_ob)
    a
  })

  Fsum <- c(Fsum,Fsum_p)

  for(nm_i in names(Cn)){
    Cn[[nm_i]] <- rbind(Cn[[nm_i]],Cn_p[[nm_i]])
  }
  for(nm_i in names(Cw)){
    Cw[[nm_i]] <- rbind(Cw[[nm_i]],Cw_p[[nm_i]])
  }

  } # end if nyp>0

    B <- wgt*N
    Nsum <- rowSums(N)
    Bsum <- rowSums(B)

  # Values by fleet  in long format (for ggplot)
  # Landings (n)
  Cn.L <- local({
    a <- as.data.frame(lapply(Cn[grepl("^Cn.L.",names(Cn))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cn.L.","",colnames(b))
    gather(data=b,key=fleet,value=Ln,names(b)[names(b)!="year"],factor_key = TRUE)

  })
  # Discards (n)
  Cn.D <- local({
    a <- as.data.frame(lapply(Cn[grepl("^Cn.D.",names(Cn))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cn.D.","",colnames(b))
    gather(data=b,key=fleet,value=Dn,names(b)[names(b)!="year"],factor_key = TRUE)

  })

  # Landings (w)
  Cw.L <- local({
    a <- as.data.frame(lapply(Cw[grepl("^Cw.L.",names(Cw))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cw.L.","",colnames(b))
    gather(data=b,key=fleet,value=Lw,names(b)[names(b)!="year"],factor_key = TRUE)

  })
  # Discards (w)
  Cw.D <- local({
    a <- as.data.frame(lapply(Cw[grepl("^Cw.D.",names(Cw))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cw.D.","",colnames(b))
    gather(data=b,key=fleet,value=Dw,names(b)[names(b)!="year"],factor_key = TRUE)

  })

  Cn.L.tot <- tapply(Cn.L[,"Ln"],Cn.L[,"year"],sum)
  Cn.D.tot <- tapply(Cn.D[,"Dn"],Cn.D[,"year"],sum)

  Cw.L.tot <- tapply(Cw.L[,"Lw"],Cw.L[,"year"],sum)
  Cw.D.tot <- tapply(Cw.D[,"Dw"],Cw.D[,"year"],sum)


## Extend data inputs and build projected dat file
  # Identify the parts of init that need to be changed
  # Identify the new parts and match them with the parts that need to change
  if(!is.null(bam2r_args)){
    bam <- do.call(bam2r,bam2r_args)
    init_p <- init_b <- bam$init

    ### Identify temporal objects in init that need to change
    # styr
    # init_styr <- init_b[grepl("^styr",names(init_b))]
    # endyr
    init_endyr <- init_b[grepl("^endyr",names(init_b))]
    # init_obs <- init_b[grepl("^obs_(L|released|cpue|cv|lenc|agec|maturity)",names(init_b))]
    # init_obs_maturity <- init_b[grepl("^obs_maturity",names(init_b))]
    #
    # # tv (time varying objects currently in Atlantic Menhaden model)
    # init_tv <- init_b[grepl("_tv$",names(init_b))]

    #### Update temporal values
    ## _endyr
    # Any _endyr for select values and for data sets that extend to the terminal
    # year of the assessment should be extended by nyp
    yr_nm_add_p <- c("endyr","endyr_rec_dev","endyr_rec_phase2","endyr_rec_spr","styr_regs","endyr_proj",
                     names(init_endyr[grepl("^endyr_([DL]|cpue)",names(init_endyr))&as.numeric(init_endyr)>=endyr])
    )
    init_p[yr_nm_add_p] <- lapply(init_p[yr_nm_add_p],function(x){paste(as.numeric(x)+nyp)})

    yrs_rec_dev <- init_p$styr_rec_dev:init_p$endyr_rec_dev
    init_p[["set_log_dev_vals_rec"]] <- setNames(rep("0.0",length(yrs_rec_dev)),yrs_rec_dev)

    ## landings
    # Note: need to make sure you get the right units (n or w)
    init_obs_L_nm <- names(init_b)[grepl("^obs_L",names(init_b))]
    for(nm_i in init_obs_L_nm){
      abb_i <- gsub("^(obs_L)(_)(.*)","\\3",nm_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]]))
      yrsbi <- names(xbi)
      xni <- rowSums(Cn[[paste0("Cn.L.",abb_i)]])
      xwi <- rowSums(Cw[[paste0("Cw.L.",abb_i)]])
      xkni <- xni/1000
      xkwi <- xwi/1000
      xcvi <- init_b[[paste0("obs_cv_L_",abb_i)]]
      # Compute the sums of squared deviations to determine whether landings are
      # in numbers or in weight. Landings in the appropriate units are then added
      # to the observed landings to add to the new dat file.
      SSxbi <- unlist(lapply(list(xni=xni,xwi=xwi,xkni=xkni,xkwi=xkwi),function(x){sum((x[yrsbi]-xbi)^2)}))
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvi))))
      xpi <- round(get(names(SSxbi)[which.min(SSxbi)])[paste(yp)],ndigi)
      xi <- c(xbi,xpi[xpi>0]) # Only include projected values greater than zero
      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      init_p[[paste0("obs_cv_L_",abb_i)]] <- setNames(paste(round(C_ts_cv[,paste0("cv.L.",abb_i)],ndigcvi)),rownames(C_ts_cv))[names(xi)]
      ## set_log_dev_vals_F_L
      init_p[[paste0("set_log_dev_vals_F_L_",abb_i)]] <- setNames(rep("0.0",length(xi)),names(xi))
    }

    ## discards
    # Note: Should be in numbers, but the function still checks both n and w
    # need to convert predicted discards
    init_obs_released_nm <- names(init_b)[grepl("^obs_released",names(init_b))]
    for(nm_i in init_obs_released_nm){
      abb_i <- gsub("^(obs_released)(_)(.*)","\\3",nm_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]]))
      yrsbi <- names(xbi)
      DMi <- as.numeric(init_b[[paste0("set_Dmort_",abb_i)]])
      xni <- rowSums(Cn[[paste0("Cn.D.",abb_i)]])*1/DMi # Convert to released (live discards)
      xwi <- rowSums(Cw[[paste0("Cw.D.",abb_i)]])*1/DMi # Convert to released (live discards)
      xkni <- xni/1000
      xkwi <- xwi/1000
      xcvi <- init_b[[paste0("obs_cv_D_",abb_i)]]
      SSxbi <- unlist(lapply(list(xni=xni,xwi=xwi,xkni=xkni,xkwi=xkwi),function(x){sum((x[yrsbi]-xbi)^2)}))
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvi))))
      xpi <- round(get(names(SSxbi)[which.min(SSxbi)])[paste(yp)],ndigi)
      xi <- c(xbi,xpi[xpi>0]) # Only include projected values greater than zero
      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      init_p[[paste0("obs_cv_D_",abb_i)]] <- setNames(paste(round(C_ts_cv[,paste0("cv.D.",abb_i)],ndigcvi)),rownames(C_ts_cv))[names(xi)]
      ## set_log_dev_vals_F_D
      init_p[[paste0("set_log_dev_vals_F_D_",abb_i)]] <- setNames(rep("0.0",length(xi)),names(xi))
    }

    ## cpue
    init_obs_cpue_nm <- local({
      nm <- names(init_b)[grepl("^obs_cpue",names(init_b))]
      nm_abb <- gsub("^(obs_cpue)(_)(.*)","\\3",nm)
      nm_yes <- nm[which(nm_abb%in%names(U))]
      nm_no <- nm[which(!nm_abb%in%names(U))]
      if(length(nm_no)>0){
        message(paste("message: ",paste(nm_no,collapse=", "), "were found in the bam tpl but not in the rdat. They may not be included in the likelihood."))
      }
      # Only include names found in U (indices reported in the rdat)
      return(nm_yes)
    })

    for(nm_i in init_obs_cpue_nm){
      abb_i <- gsub("^(obs_cpue)(_)(.*)","\\3",nm_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]]))
      yrsbi <- names(xbi)
      xproji <- setNames(U[,abb_i],rownames(U))
      xcvi <- init_b[[paste0("obs_cv_cpue_",abb_i)]]
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvi))))

      xpi <- round(xproji[paste(yp)],ndigi)
      xi <- c(xbi,xpi[!is.na(xpi)])
      yrsi <- names(xi)
      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      init_p[[paste0("obs_cv_cpue_",abb_i)]] <- setNames(paste(round(U_ts_cv[,paste0("cv.U.",abb_i)],ndigcvi)),rownames(U_ts_cv))[names(xi)]

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]
    }

    ## agec
    init_obs_agec_nm <- names(init_b)[grepl("^obs_agec",names(init_b))]
    init_obs_agec_nm_key <- setNames(names(acomp),
                                     paste0("obs_agec_",gsub("^(D.)([A-Za-z]+)","\\2_D",gsub("^L.","",names(acomp)))))
    for(nm_i in init_obs_agec_nm){
      nm2_i <- init_obs_agec_nm_key[[nm_i]]
      xbi <- init_b[[nm_i]]
      x2bi <- acomp[[nm2_i]]
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      obsi <- apply(x2bi,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
      dimnames(obsi) <- dimnames(x2bi)
      yrsi <- rownames(obsi)
      init_p[[nm_i]] <- obsi

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]

      # update nfish and nsamp as necessary
      nfish_nm_i <- gsub("obs","nfish",nm_i)
      nsamp_nm_i <- gsub("obs","nsamp",nm_i)
      nfishbi <- init_b[[nfish_nm_i]]
      nsampbi <- init_b[[nsamp_nm_i]]
      nfishpi <- setNames(rep(bamExtras::geomean2(nfishbi[paste(tail(yb,nyb$comp))]),nyp),paste(yp))
      nsamppi <- setNames(rep(bamExtras::geomean2(nsampbi[paste(tail(yb,nyb$comp))]),nyp),paste(yp))
      nfishdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nfishbi))))
      nsampdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nsampbi))))
      nfishi <- setNames(sprintf(paste0("%.",nfishdigi,"f"),round(as.numeric(c(nfishbi,nfishpi)))),names(c(nfishbi,nfishpi)))
      nsampi <- setNames(sprintf(paste0("%.",nsampdigi,"f"),round(as.numeric(c(nsampbi,nsamppi)))),names(c(nsampbi,nsamppi)))

      init_p[[nfish_nm_i]] <- nfishi[yrsi]
      init_p[[nsamp_nm_i]] <- nsampi[yrsi]
    }

    ## lenc
    init_obs_lenc_nm <- names(init_b)[grepl("^obs_lenc",names(init_b))]
    init_obs_lenc_nm_key <- setNames(names(lcomp),
                                     paste0("obs_lenc_",gsub("^(D.)([A-Za-z]+)","\\2_D",gsub("^L.","",names(lcomp)))))
    # for(nm_i in init_obs_lenc_nm){
    #   nm2_i <- init_obs_lenc_nm_key[[nm_i]]
    #   xbi <- init_b[[nm_i]]
    #   x2bi <- lcomp[[nm2_i]]
    #   ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
    #   obsi <- apply(x2bi,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
    #   dimnames(obsi) <- dimnames(x2bi)
    #   yrsi <- rownames(obsi)
    #   init_p[[nm_i]] <- obsi
    #   # update styr, endyr, yrs, and nyr as necessary
    #   yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
    #                       paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
    #   )
    #   yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
    #   init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]
    #
    #   # update nfish and nsamp as necessary
    #   nfish_nm_i <- gsub("obs","nfish",nm_i)
    #   nsamp_nm_i <- gsub("obs","nsamp",nm_i)
    #   nfishbi <- init_b[[nfish_nm_i]]
    #   nsampbi <- init_b[[nsamp_nm_i]]
    #   nfishdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nfishbi))))
    #   nsampdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nsampbi))))
    #   nfish2i <- setNames(sprintf(paste0("%.",nfishdigi,"f"), nfish[,paste0("lcomp.",nm2_i)]),rownames(nfish))
    #   nsamp2i <- setNames(sprintf(paste0("%.",nsampdigi,"f"), ntrip[,paste0("lcomp.",nm2_i)]),rownames(ntrip))
    #
    #
    #   init_p[[nfish_nm_i]] <- nfish2i[yrsi]
    #   init_p[[nsamp_nm_i]] <- nsamp2i[yrsi]
    # }
    for(nm_i in init_obs_lenc_nm){
      nm2_i <- init_obs_lenc_nm_key[[nm_i]]
      xbi <- init_b[[nm_i]]
      x2bi <- lcomp[[nm2_i]]
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      obsi <- apply(x2bi,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
      dimnames(obsi) <- dimnames(x2bi)
      yrsi <- rownames(obsi)
      init_p[[nm_i]] <- obsi

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]

      # update nfish and nsamp as necessary
      nfish_nm_i <- gsub("obs","nfish",nm_i)
      nsamp_nm_i <- gsub("obs","nsamp",nm_i)
      nfishbi <- init_b[[nfish_nm_i]]
      nsampbi <- init_b[[nsamp_nm_i]]
      nfishpi <- setNames(rep(bamExtras::geomean2(nfishbi[paste(tail(yb,nyb$comp))]),nyp),paste(yp))
      nsamppi <- setNames(rep(bamExtras::geomean2(nsampbi[paste(tail(yb,nyb$comp))]),nyp),paste(yp))
      nfishdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nfishbi))))
      nsampdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nsampbi))))
      nfishi <- setNames(sprintf(paste0("%.",nfishdigi,"f"),round(as.numeric(c(nfishbi,nfishpi)))),names(c(nfishbi,nfishpi)))
      nsampi <- setNames(sprintf(paste0("%.",nsampdigi,"f"),round(as.numeric(c(nsampbi,nsamppi)))),names(c(nsampbi,nsamppi)))

      init_p[[nfish_nm_i]] <- nfishi[yrsi]
      init_p[[nsamp_nm_i]] <- nsampi[yrsi]
    }

    bam_p <- bam2r(dat_obj=bam$dat,tpl_obj=bam$tpl,cxx_obj=bam$cxx,init=init_p)
  }else{
    bam_p <- NULL
  } # end if(!is.null(bam2r_args))

  # plot stuff
  if(plot){
  par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
  # N
  plot(names(Nsum),Nsum,type="o")
  abline(v=endyr)
  # R
  plot(names(R),R,type="o")
  abline(v=endyr)
  # B
  plot(names(Bsum),Bsum,type="o")
  abline(v=endyr)
  # Fsum
  plot(names(Fsum),Fsum,type="o")
  abline(v=endyr)

  # Landings and discards
  # Cn.L
  p <- ggplot(Cn.L,mapping=aes(x=year,y=Ln))+
    geom_area(aes(fill=fleet))+
    theme_bw()+
    scale_fill_brewer(palette="Spectral")+
    stat_summary(fun = sum, geom = "line", size = 1)+
    stat_summary(fun = sum, geom = "point", size = 2)+
    geom_vline(xintercept = endyr, linetype="dashed", size = 0.3)
  p2 <- p + geom_text(aes(x=endyr, label="endyr\n",y=max(Ln)), angle=90)
  print(p2)


  # Cn.D
  p2 <- p %+% Cn.D + aes(y=Dn) +
    geom_text(aes(x=endyr, label="endyr\n",y=max(Dn)), angle=90)
  print(p2)


  # Cw.L
  p2 <- p %+% Cw.L + aes(y=Lw) +
    geom_text(aes(x=endyr, label="endyr\n",y=max(Lw)), angle=90)
  print(p2)


  # Cw.D
  p2 <- p %+% Cw.D + aes(y=Dw) +
    geom_text(aes(x=endyr, label="endyr\n",y=max(Dw)), angle=90)
  print(p2)


  # cpue
    matplot(as.numeric(rownames(U)),U,type="o",xlab="",xlim=c(styr,endyr+nyp),pch=1)
    # matpoints(as.numeric(rownames(U_p)),U_p,type="o",pch=1)
    legend("topleft",legend=colnames(U),col=1:ncol(U),lty=1:ncol(U),pch=1)
    abline(v=endyr)
  }

  # Return results

  invisible(list(
              # years_p=yp,
              # F=Fsum_p,
              # L_wgt=Lw_p,
              # L_num=Ln_p,
              # D_wgt=Dw_p,
              # D_num=Dn_p,
              # S=S_p,
              # B=B_p,
              # R=R_p,
              # N=N_p,
              bam_p = bam_p
              )
         )
}
