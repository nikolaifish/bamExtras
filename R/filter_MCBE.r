#' Identify which runs should be filtered out of results of MCBE uncertainty analysis
#'
#' @param sim_summary Output object from summarize_MCBE
#' @param filter_info list of ranges used to filter out MCBE results. set to NULL to retain all results.
#' @returns This function returns a list including: 1. \code{sim_pass} is a
#' numeric vector indicating which simulation runs should be retained after filtering.
#' 2. \code{test_vals} is a data frame of the values used to conduct filtering for all runs.
#' 3. \code{test_logical} is a logical matrix indicating the results of each
#' filtering test.
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer, Erik Williams, and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Run MCBE, writing files to dir_bam_sim
#' run_MCBE("GrayTriggerfish", dir_bam_sim="sim_GrTr")
#' # Summarize results of MCBE and assign to object
#' ss_GrTr <- summarize_MCBE(dir_bam_sim="sim_GrTr")
#' sim2keep <- filter_MCBE(ss_GrTr,nm_rp = list("spr.brps"=c("Fref"="F30","F.Fref"="Fend.F30.mean")))
#' # Plot MCBE results
#' }

filter_MCBE <- function(sim_summary,
                        filter_info = list("gradient.max"  = c(0,0.01),
                                           "Fref"          = c(0,5),
                                           "F.Fref"        = c(0,5),
                                           "R.sigma.par"   = c(0.2,1.0),
                                           "R0_prob"       = c(0.005,0.995), # Not limits but probabilities used to compute limits
                                           "avoid_bounds"  = TRUE # logical. If TRUE, filter out runs with parameters near the bounds
                        ),
                        nm_rp = list("parms"=c("Fref"="Fmsy","F.Fref"="Fend.Fmsy.mean"))
                        ){
  ss <- sim_summary
  fi <- filter_info

  parms <- ss$parms
  parm.cons <- ss$parm.cons
  a.series <- ss$a.series

  nsim <- length(parms$modRunName)
  simID2 <- simID <- 1:nsim

  nearbound <- parms$nearbound
  Nnearbound <- length(which(nearbound))
  nm_R0 <- local({
    a <- c("BH.R0","R0")
    b <- a[a%in%names(ss$parms)]
    if(length(b)==0){
      message("neither BH.R0 nor R0 found in names(parms). R0_prob was removed from filter_info.")
      fi <- fi[[names(fi)!="R0_prob"]]
    }else{
    b[1]
    }
  })
  R0 <- parms[[nm_R0]]
  R.sigma.par <- parms$R.sigma.par
  gradient.max <- ss$like$gradient.max
  for(i in seq_along(nm_rp)){
    nm_i <- names(nm_rp)[i]
    nm_rp_i <- nm_rp[[nm_i]]
    xi <- ss[[nm_i]]
    for(j in seq_along(nm_rp_i)){
      nm_rp_ij <- nm_rp_i[j]
      xij <- xi[,nm_rp_ij]
      assign(names(nm_rp_ij),xij)
    }
  }

  nm_fi <- names(fi)

  # Compute R0 limits from R0_prob
  R0_lim <- c(0,1)
  if("R0_prob"%in%nm_fi){
    R0_lim <- as.numeric(quantile(R0, probs=c(fi$R0_prob[1],fi$R0_prob[2])))
    fi$R0 <- R0_lim
    fi <- fi[which(names(fi)!="R0_prob")]
  }

  nm_fi <- names(fi)

  nm_filter <- c("gradient.max","Fref","F.Fref","R.sigma.par","R0")
  mt_pass <- matrix(TRUE,nrow=nsim,ncol=length(nm_filter),dimnames=list("sim"=simID,"par"=nm_filter))
  mt_vals <- as.data.frame(mt_pass*NA)
  for(nm_i in nm_filter){
    if(nm_i%in%nm_fi){
       xi_lim <- fi[[nm_i]]
       xi_obs <- get(nm_i)
      if(length(xi_lim)!=2){warning(paste("filtering criteria for",nm_i,"must be a numeric vector of length 2."))}
       mt_pass[,nm_i] <- (xi_obs>xi_lim[1] & xi_obs<xi_lim[2])
       mt_vals[,nm_i] <- xi_obs
    }
  }

  if("avoid_bounds"%in%nm_fi){
    mt_vals <- cbind(mt_vals,"nearbound"=nearbound)
    mt_pass <- cbind(mt_pass,"not_nearbound"=!nearbound)
  }

    sim_pass <- which(apply(mt_pass,1,all))
    simID2 <- simID[sim_pass]
  return(list("sim_pass"=simID2,"test_vals"=mt_vals,"test_logical"=mt_pass))
}
