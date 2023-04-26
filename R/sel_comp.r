#' sel_comp
#'
#' Estimate selectivity from age or length compositions using optim()
#' @param par_init named vector of initial guesses for parameters to estimate
#' @param par_lo lower bounds for par
#' @param par_up upper bounds for par
#' @param fn_sel name of selectivity function which can be called with two arguments,
#' "x" and "par", where "x" is a vector of ages and "par" is a named numeric vector
#' where the names correspond the names of the parameters in fn_sel.
#' @param optim_control control argument to pass to \code{\link[stats]{optim}}
#' @param optim_method method argument to pass to \code{\link[stats]{optim}}
#' @param age_pop vector of ages in population (often includes 0)
#' @param age_agec vector of ages in age comp data (often does not include 0)
#' @param P_agec vector of proportion-at-age in age comp data
#' @param P_lenc vector of proportion-at-age in length comp data
#' @param lenprob age-length conversion matrix. Only used if P_lenc is specified
#' @param F_full fixed value of full F. numeric vector of length 1.
#' @param M_age vector of natural mortality at age, corresponding to age_agec
#' @param N0 initial population size
#' @param penalty_out_of_bounds large value to add as penalty to objective function value (NLL or SSQ) if parameter estimates go outside of bounds
#' @param dbeta_args list of arguments to pass to dbeta() function for defining penalty function.
#' The shape parameters should be equal for symmetrical U-shaped penalty profile, which approach
#' infinity as the parameter nears the bounds. Smaller values of the shape parameters make the U-shape
#' increasingly flat-bottomed while larger values result in more curvature.
#' @param value_type "SSQ" to return sum of squared deviations (default). "NLL" to return negative log-likelihood.
#' @param output_type "value" (default) returns objective function value (NLL or SSQ). Use "value" with optim(). "predict" to return predicted age comp values.
#' @param pass_dat_par_trace logical. Should dat_par_trace be passed to the global environment? Allows the user to follow the optimization of par.
#' @param penalty_scale value to multiply by penalties added to objective function value. Defaults to 1 which leaves penalties as is. A value of zero effectively removes all penalties.
#' @param fit_to_zeros logical. Should observed comps with values of zero be included in the objective function?
#' @param plot_format list of plot formatting defaults
#' @details As of 2023-04-20, this function seems to work pretty well based on testing with red snapper data (see example).
#' However fitting to age comps works better than length comps. Overall, fitting seems to be sensitive
#' to things like dbeta_args, penalty_out_of_bounds, and the variability of length at age in lenprob.
#' So, at this point I'm going to use it with caution and continue to improve it.
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' rdat <- rdat_RedSnapper
#' agec <- tail(rdat$comp.mats$acomp.sCT.ob,1)
#' lenc <- tail(rdat$comp.mats$lcomp.cHL.ob,1)
#'
#' # Estimate selectivity
#' sela1 <- sel_comp(age_pop=rdat$a.series$age,
#'                     age_agec=as.numeric(dimnames(agec)[[2]]),
#'                     P_agec=as.numeric(agec),
#'                     M_age=rdat$a.series$M
#' )
#' # Run the code again with the estimates from sela1 for par_init
#' sela2 <- sel_comp(par_init = sela1$fit$par,
#'                     age_pop=rdat$a.series$age,
#'                     age_agec=as.numeric(dimnames(agec)[[2]]),
#'                     P_agec=as.numeric(agec),
#'                     M_age=rdat$a.series$M
#' )
#'
#' # Build length-at-age probability matrix
#' lenprob <- agelenprob(
#'   age=rdat$a.series$age,
#'   len=rdat$a.series$length,
#'   binmid_lenc=as.numeric(colnames(rdat$comp.mats$lcomp.cHL.ob)),
#'   cv_lenc=rdat$parm.cons$len.cv.val.L[8],
#'   plot=TRUE)
#'
#' sell1 <- sel_comp(
#'   age_pop=rdat$a.series$age,
#'   age_agec=as.numeric(dimnames(agec)[[2]]),
#'   P_lenc=as.numeric(lenc),
#'   M_age=rdat$a.series$M,
#'   lenprob=lenprob
#' )
#'
#' sell2 <- sel_comp(par_init = sell1$fit$par,
#'                     age_pop=rdat$a.series$age,
#'                     age_agec=as.numeric(dimnames(agec)[[2]]),
#'                     P_lenc=as.numeric(lenc),
#'                     M_age=rdat$a.series$M,
#'                     lenprob=lenprob
#' )
#' }

sel_comp <- function(
  par_init = c("a1"=0.5, "b1"=3, "a2"=0.5,"b2"=6),
  par_lo=par_init*0.001,
  par_up=par_init*1000,
  fn_sel="dbllgs",
  optim_method = "Nelder-Mead",
  optim_control=list(reltol=sqrt(.Machine$double.eps),  maxit=1000),
  age_pop,
  age_agec,
  P_agec=NULL,
  P_lenc=NULL,
  lenprob,
  F_full=0.2,
  M_age,
  N0 = 1,
  penalty_out_of_bounds=100,
  dbeta_args=list(shape1=10e-6, shape2=10e-6),
  value_type="SSQ",
  output_type="value",
  pass_dat_par_trace=TRUE,
  penalty_scale=1,
  fit_to_zeros=FALSE,
  plot_format =  list(
    col=list("init"="green2","obs"="purple","pr"="black"),
    lty=list("init"=2,"obs"=0,"pr"=1),
    lwd=list("init"=2,"obs"=2,"pr"=2),
    pch=list("init"=NA,"obs"=1,"pr"=NA),
    type=list("init"="l","obs"="p","pr"="l"))
){

  # replace function name with actual function
  fn_sel <- get(fn_sel)

  if(!is.null(P_agec)&!is.null(P_lenc)){
    warning("P_agec and P_lenc are both specified. Please specify only one of these. This function won't produce results for both agec and lenc at the same time.")
  }

  if(!exists("env_par_trace",envir=.GlobalEnv,mode="environment")){
    # Define environment for saving optim() trace information
    assign("env_par_trace",new.env(),envir=.GlobalEnv)
  }else{
    # If it exists, clear it
    rm(list=ls(envir=env_par_trace),envir=env_par_trace)
  }

  objfn <- function(par,
                    output_type="value"){

    sel <-  fn_sel(x=age_pop,par=par)
    sel_scaled <- sel/max(sel)

    # Mortality
    F_age <- sel_scaled*F_full
    Z_age <- M_age + F_age

    # Numbers-at-age
    N_age <- local({
      N_age <- age_pop*NA
      N_age[1] <- N0
      for(i in 2:(length(N_age))){
        N_age[i] <- N_age[i-1]*exp(-Z_age[i-1])
      }
      N_age[length(N_age)] <- N_age[length(N_age)]/(1-exp(-Z_age[length(N_age)]))
      N_age
    })

    if(!is.null(P_agec)){
      # Predicted age comp
      P_agec_pr <- local({
        a <- N_age*sel/sum(N_age*sel)
        names(a) <- age_pop
        b <- a[paste(age_agec)]
        c <- b/sum(b)
        c
      })

      P_agec_pr <- pmin(pmax(P_agec_pr,0),1) # keep P_agec_pr between zero and 1
    }

    if(!is.null(P_lenc)){
      # Number-at-length
      N_len <- colSums(N_age*lenprob)

      # Predicted length comp
      P_lenc_pr <-  local({
        a <- N_age*sel/sum(N_age*sel)
        names(a) <- age_pop
        b <- colSums(a*lenprob)
        c <- b/sum(b)
        c
      })

      P_lenc_pr <- pmin(pmax(P_lenc_pr,0),1) # keep P_lenc_pr between zero and 1
    }

    # Add penalties to NLL for each parameter to keep each away from bounds
    pnl_nearbound <- unlist(
      lapply(names(par),function(nm_p){
        p <- par[[nm_p]]
        lo_p <- par_lo[[nm_p]]
        up_p <- par_up[[nm_p]]
        do.call(dbeta,c(list(x=as.numeric((p-lo_p)/(up_p-lo_p))),dbeta_args))
      })
    )

    pnl_outofbound <- unlist(
      lapply(names(par),function(nm_p){
        p <- par[[nm_p]]
        lo_p <- par_lo[[nm_p]]
        up_p <- par_up[[nm_p]]
        ifelse(p<lo_p|p>up_p,penalty_out_of_bounds,0)
      })
    )

    penalties <- sum(c(pnl_nearbound,pnl_outofbound))*penalty_scale

    # Compute objective functions
    if(!is.null(P_agec)){
      if(fit_to_zeros){
        P_agec_o <- P_agec
        P_agec_pr_o <- P_agec_pr
      }else{
        P_agec_o <-    P_agec[which(P_agec>0)]
        P_agec_pr_o <- P_agec_pr[which(P_agec>0)]
      }

      agec_sd_resid <- mean(sqrt((P_agec_o-P_agec_pr_o)^2))
      NLL <- -sum(log(pmax(dnorm(x=P_agec_o, mean = P_agec_pr_o, sd=agec_sd_resid, log = FALSE),1e-15))) + penalties

      SSQ <-  sum((P_agec_pr_o-P_agec_o)^2) + penalties
    }
    if(!is.null(P_lenc)){
      if(fit_to_zeros){
        P_lenc_o <- P_lenc
        P_lenc_pr_o <- P_lenc_pr
      }else{
        P_lenc_o <-    P_lenc[which(P_lenc>0)]
        P_lenc_pr_o <- P_lenc_pr[which(P_lenc>0)]
      }

      lenc.sd.resid <- mean(sqrt((P_lenc_o-P_lenc_pr_o)^2))
      NLL <- -sum(log(pmax(dnorm(x=P_lenc_o, mean = P_lenc_pr_o, sd=lenc.sd.resid, log = FALSE),1e-15))) + penalties

      SSQ <-  sum((P_lenc_pr_o-P_lenc_o)^2) + penalties
    }

    if(value_type=="NLL"){value <- NLL}
    if(value_type=="SSQ"){value <- SSQ}

    # Get parameter values from each iteration of the optimizer
    if(!exists("dat_par_trace",envir = env_par_trace,inherits=FALSE)){
      env_par_trace$dat_par_trace <- as.data.frame(matrix(numeric(0),ncol=length(par)+1,dimnames=list(NULL,c(names(par),"value")))) # initialize empty list
    }
    trace_step_i <- as.data.frame(t(c(par,"value"=value)))
    env_par_trace$dat_par_trace <- rbind(env_par_trace$dat_par_trace,
                                         trace_step_i
    )

    # Return value
    if(output_type=="predict"){
      if(!is.null(P_agec)){
        out <- P_agec_pr
      }
      if(!is.null(P_lenc)){
        out <- P_lenc_pr
      }
    }

    if(output_type=="value"){
      out <- value
    }
    out
  }  # end objfn

  # Compute selectivity values based on initial parameter (par_init) values
  sel_scaled_init <- local({
    sel <- fn_sel(x=age_pop,par=par_init)
    sel_scaled <- sel/max(sel)
  })

  if(!is.null(P_agec)){
    # Calculate predicted age comps based on initial parameter (par_init) values
    agec_init <- objfn(
      par=par_init,
      output_type="predict"
    )


    # Fit selectivity parameters
    fit <- optim(
      fn=objfn,
      par=par_init,
      control=optim_control,
      method=optim_method
    )

    # Calculate predicted age comps based on fitted selectivity
    agec_pr <- objfn(
      par=fit$par,
      output_type="predict")
  }

  if(!is.null(P_lenc)){
    # Calculate predicted age comps based on initial parameter (par_init) values
    lenc_init <- objfn(
      par=par_init,
      output_type="predict"
    )


    # Fit selectivity parameters
    fit <- optim(
      fn=objfn,
      par=par_init,
      control=optim_control,
      method=optim_method
    )

    # Calculate predicted age comps based on fitted selectivity
    lenc_pr <- objfn(
      par=fit$par,
      output_type="predict")
  }

  # Compute predicted selectivity scaled so that the maximum value is equal to one
  sel_scaled_pr <- local({
    sel <- fn_sel(x=age_pop,par=fit$par)
    sel_scaled <- sel/max(sel)
  })

  # Plot stuff
  env_here <- environment()
  lapply(seq_along(plot_format),function(i){
    x <- plot_format[i]
    assign(names(x),x[[1]],envir=env_here)
  })

  par(mfrow=c(2,1),mar=c(2,2,2,1),oma=c(0,0,3,0),mgp=c(1,0.2,0),tck=-0.01,lend="butt")

  # Plot initial selectivity
  plot(x=age_pop,  y=sel_scaled_init,
       col=col$init, lwd=lwd$init, lty=lty$init, pch=pch$init, type=type$init,
       main="selectivity", xlab="age",ylab="selectivity",ylim=c(0,1))
  # add predicted selectivity
  points(x=age_pop, y=sel_scaled_pr,
         col=col$pr, lwd=lwd$pr, lty=lty$pr, pch=pch$pr, type=type$pr
  )
  grid()

  # Plot comps

  if(!is.null(P_agec)){
    x <- age_agec
    y_obs <- P_agec
    y_init <- agec_init
    y_pr <- agec_pr

    xlab <- "age"
    main <- "age comp"
  }
  if(!is.null(P_lenc)){
    x <- as.numeric(colnames(lenprob))
    y_obs <- P_lenc
    y_init <- lenc_init
    y_pr <- lenc_pr

    xlab <- "length"
    main <- "length comp"
  }

  # Plot observed comps
  plot(x,y_obs,ylim=range(c(y_obs,y_init,y_pr)),
       col=col$obs, lwd=lwd$obs, lty=lty$obs, pch=pch$obs, type=type$obs,
       main=main, xlab=xlab, ylab="proportion")

  # Plot predicted comps based on par_init
  points(x,y_init,
         col=col$init, lwd=lwd$init, lty=lty$init, pch=pch$init, type=type$init
  )

  # Plot predicted comps based on fitted selectivity
  points(x,y_pr,
         col=col$pr, lwd=lwd$pr, lty=lty$pr, pch=pch$pr, type=type$pr
  )
  grid()

  # if(!is.null(title_text)){
  #   title(main=title_text,outer=TRUE)
  # }

  legend("right",bty="n",
         legend=c("initial","observed","predicted"),
         col=unlist(col), lwd=unlist(lwd), lty=unlist(lty), pch=unlist(pch)
  )

  # Return values
  dat_sel <- data.frame("age"=age_pop,"sel_scaled_init"=sel_scaled_init,"sel_scaled_pr"=sel_scaled_pr)
  dat_comp <- data.frame("x"=x,"comp_init"=y_init,"comp_obs"=y_obs,"comp_pr"=y_pr)
  names(dat_comp) <- gsub("x",xlab,names(dat_comp))

  invisible(list("par_init"=par_init,
                 "sel"=dat_sel,
                 "comp"= dat_comp,
                 "fit"=fit,
                 "dat_par_trace"=env_par_trace$dat_par_trace
  )
  )
}
