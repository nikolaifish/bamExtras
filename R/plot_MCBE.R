#' Plot results from MCBE uncertainty analysis
#'
#' @param sim_summary Output object from summarize_MCBE
#' @param bamout
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Run MCBE, writing files to dir_bam_sim
#' run_MCBE("GrayTriggerfish", dir_bam_sim="sim_GrTr")
#' # Summarize results of MCBE and assign to object
#' MCBE_GrTr <- summarize_MCBE(dir_bam_sim="sim_GrTr")
#' # Plot MCBE results
#' }



plot_MCBE <- function(sim_summary,
                      base_out=NULL,
                      dir_figs = "figs",
                      filter_info = list(gradient.max       = c(0,0.01),
                                        Fmsy           = c(0,5),
                                        Fend.Fmsy.mean = c(0,5),
                                        R.sigma.par    = c(0.2,1.0),
                                        R0_prob        = c(0.005,0.995) # Not limits but probabilities used to compute limits
                      )
                      ){
ss <- sim_summary
# test args

### Do stuff with base run
# spp <- dget(file.path(dir_bam,paste(filename,'rdat', sep=".")))


# nm_sim <- sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim)
#
#
# ## Get base model stuff
#
# # Identify working directory
#   wd <- getwd()
#   message(paste("working directory:",wd))
#

  nsim <- length(ss$parms$modRunName)
  simID2 <- simID <- 1:nsim

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ####   plotsConvergence   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # pdf(file.path(dir_figs,"Fig-MCB-converge.pdf"))
    with(ss$parms,{

      sd.eda1 <- sd.eda2 <- sd.eda3 <- rep(NA,nsim)
      for(i in 1:nsim)
      {
        sd.eda1[i]=sd(msy.klb[1:i])
        sd.eda2[i]=sd(Fmsy[1:i])
        sd.eda3[i]=sd(SSBmsy[1:i])
      }
      par(mfrow=c(3,1),mar=c(2,2,1,1),mgp=c(1.2,.3,0))

      plot(1:nsim, sd.eda1,xlab="", ylab="SE(msy)", type="l", lwd=2)
      plot(1:nsim, sd.eda2,xlab="", ylab="SE(Fmsy)", type="l", lwd=2)
      plot(1:nsim, sd.eda3,xlab="Number bootstrap replicates", ylab="SE(SSBmsy)", type="l", lwd=2)
    })
    # dev.off()

  #   pdf(file.path(dir_figs,"lk.total.v.gradient.max.pdf"))
    par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=c(1.2,.3,0))
    x <- ss$like$lk.total
    y <- ss$like$gradient.max
    plot(x=x,
         y=y,
         pch=16,col=rgb(0,0,0,0.2),
         xlab="lk.total",ylab="gradient.max",
         xlim=range(x)+c(-1,1),ylim=range(y)+diff(range(y))*c(0,.1))
    # points(x=spp$like[["lk.total"]],
    #        y=spp$like[["gradient.max"]],col="red",pch=3)
    legend("topright",legend=paste("n=",length(x)),bty="n")

    plot_boot_density(ss$like,"lk.total",plot_args = list(xlab = "lk.total", ylab="density",main = "", type = "l"))
    # abline(v=spp$like[["lk.total"]],col="red",lwd=2)
    plot_boot_density(ss$like,"gradient.max",plot_args = list(xlab = "gradient.max", ylab="density",main = "", type = "l"))

    # abline(v=spp$like[["gradient.max"]],col="red",lwd=2)
    # dev.off()

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ####       trimRuns       ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # trim some runs that didn't converge or params hit upper bounds
    nearbound <- ss$parms$nearbound
    Nnearbound <- length(which(nearbound))
    R0 <- ss$parms$BH.R0
    R.sigma.par <- ss$parms$R.sigma.par
    gradient.max <- ss$like$gradient.max
    Fend.Fmsy.mean <- ss$parms$Fend.Fmsy.mean
    Fmsy <- ss$parms$Fmsy
    # q.cL <- ss$parms$q.cL
    # q.sM <- ss$parms$q.sM

    R0_lim <- quantile(R0, probs=c(filter_info$R0_prob[1],filter_info$R0_prob[2]))

    # simID2 <- simID[which(gradient.max<0.01 & R0>=R0.trim[1] & R0<=R0.trim[2] & !nearbound)] # Include runs that meet certain criteria
    # Most of these criteria for trimming MCB runs are from the SEDAR 25 2016 update, except for the !nearbound which I added in SEDAR 66 (NK 2021-04-02)
    # Several commented out were used in the SEDAR 25 2017 update
    filter_tests <- data.frame(gradient.max = gradient.max>filter_info$gradient.max[1] & gradient.max<filter_info$gradient.max[2],
                               Fend.Fmsy.mean = Fend.Fmsy.mean > filter_info$Fend.Fmsy.mean[1] & Fend.Fmsy.mean < filter_info$Fend.Fmsy.mean[2],
                               Fmsy = Fmsy > filter_info$Fmsy[1] & Fmsy < filter_info$Fmsy[2],
                               R.sigma.par = R.sigma.par > filter_info$R.sigma.par[1] & R.sigma.par < filter_info$R.sigma.par[2],
                               R0 = R0 > R0_lim[1] & R.sigma.par < R0_lim[2],
                               notNearBound = !nearbound
                               )
    sim_pass <- which(apply(filter_tests,1,all))
    simID2 <- simID[sim_pass]

    # pdf(file.path(dir_figs,"gradient.max.density.pdf"))
    par(mfrow=c(1,1))
    plot_boot_density(data=ss$like[simID2,],x_name="gradient.max")
    # dev.off()

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ####     plotMCBData      ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # if(diagnostics.i=="MCB"){
    base.t.series <- spp$t.series[spp$t.series$year%in%yr.plot,] # Remove last year (projection) from t.series
    # Plot input data
    # Landings
    # pdf(file.path(dir_figs,"Fig-MCB-landings.ob.pdf"))

    names_L_ob <- names(ss$t.series)[grepl("^L.*ob$",names(ss$t.series))]

    par(mfrow=c(length(names_L_ob),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    for(i in names_L_ob){
      obj.i <- ss$t.series[[i]][simID2,]
      plot_boot_vec(obj.i,xlab="",ylab=i)
    }
    # dev.off()

    # Discards
    # pdf(file.path(dir_figs,"Fig-MCB-discards.ob.pdf"))

    names_D_ob <- names(ss$t.series)[grepl("^D.*ob$",names(ss$t.series))]

    par(mfrow=c(length(names_D_ob),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    for(i in names_D_ob){
      obj.i <- ss$t.series[[i]][simID2,]
      plot_boot_vec(obj.i,xlab="",ylab=i)
    }
    # dev.off()

    # Indices of abundance
    # pdf(file.path(dir_figs,"Fig-MCB-indices.ob.pdf"))

    names_U_ob <- names(ss$t.series)[grepl("^U.*ob$",names(ss$t.series))]

    par(mfrow=c(length(names_U_ob),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    for(i in names_U_ob){
      obj.i <- ss$t.series[[i]][simID2,]
      plot_boot_vec(obj.i,xlab="",ylab=i)
    }
    # dev.off()

  #   # Indices of abundance
  #   pdf(file.path(dir_figs,"Fig-MCB-indices.ob.pdf"))
  #   D.boot.U.names <- ls()[grepl("D.boot.U",ls())&grepl(".ob",ls())]
  #   par(mfrow=c(length(D.boot.U.names),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
  #   for(i in D.boot.U.names){
  #     obj.i <- get(i)
  #     obj.i[obj.i==-99999] <- NA
  #
  #     obj.varname.i <- gsub("D.boot.","",i)
  #
  #     y.ref.i <- base.t.series[,obj.varname.i]
  #     y.ref.i[which(y.ref.i==-99999)] <- NA
  #
  #     plotBootSeries(dataBoot=obj.i[,paste(yr.plot)],runsToKeep = simID2,
  #                    x.ref=base.t.series$year,y.ref=y.ref.i,
  #                    xlab="",ylab=obj.varname.i)
  #   }
  #   #
  #   # plotBootSeries(dataBoot=D.boot.U.sCT.ob[,names(which(!is.na(colSums(D.boot.U.sCT.ob[,paste(yr.plot)]))))],
  #   #                runsToKeep = simID2,
  #   #                x.ref=base.t.series$year,y.ref=base.t.series$U.sCT.ob,
  #   #                xlab="",ylab="Chevron trap/video index",xlim=range(yr.plot))
  #   # plotBootSeries(dataBoot=D.boot.U.rHb.ob[,names(which(!is.na(colSums(D.boot.U.rHb.ob[,paste(yr.plot)]))))],
  #   #                runsToKeep = simID2,
  #   #                x.ref=base.t.series$year,y.ref=base.t.series$U.rHb.ob,
  #   #                xlab="",ylab="headboat index",xlim=range(yr.plot))
  #   dev.off()
  #
  #   # Natural mortality
  #   pdf(file.path(dir_figs,"Fig-MCB-naturalMortality.ob.pdf"))
  #   par(mfrow=c(2,1),mar=c(3,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="M.constant",xlab="M constant")
  #   abline(v=spp$parms$M.constant,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"M.constant"]),lty=2,lwd=2)
  #
  #   plotBootSeries(dataBoot=D.boot.M,runsToKeep = simID2,
  #                  x.ref=spp$a.series$age,y.ref=spp$a.series$M,
  #                  xlab="age",ylab="natural mortality")
  #   dev.off()
  #
  #   # Discard mortality
  #   # pdf(file.path(dir_figs,"Fig-MCB-discardMortality.ob.pdf"))
  #   # par(mfrow=c(2,1),mar=c(3,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
  #   #
  #   # plotBootDensity(data=D.boot.parms[simID2,],par.name="D.mort.cHl",xlab="commercial handline")
  #   # abline(v=spp$parms$D.mort.cHl,lty=1,lwd=2)
  #   # abline(v=median(D.boot.parms[simID2,"D.mort.cHl"]),lty=2,lwd=2)
  #   #
  #   # plotBootDensity(data=D.boot.parms[simID2,],par.name="D.mort.rHb",xlab="recreational")
  #   # abline(v=spp$parms$D.mort.rHb,lty=1,lwd=2)
  #   # abline(v=median(D.boot.parms[simID2,"D.mort.rHb"]),lty=2,lwd=2)
  #   #
  #   # dev.off()
  #
  # }
  #
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # ####     plotTSeries      ####
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # if(diagnostics.i=="MCB"){
  #   par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-apicalF.pdf"))
  #   plotBootSeries(dataBoot=D.boot.F.full[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$F.full,
  #                  xlab="",ylab="Apical F")
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-FdFmsy.pdf"))
  #   plotBootSeries(dataBoot=D.boot.F.Fmsy[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$F.Fmsy,
  #                  xlab="",ylab="F/Fmsy")
  #   abline(h=1,lwd=2)
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-N.pdf"))
  #   plotBootSeries(dataBoot=D.boot.N[,paste(yr.plot)], runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$N,
  #                  xlab="",ylab="Abundance (N)")
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-SSB.pdf"))
  #   plotBootSeries(dataBoot=D.boot.SSB[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$SSB,
  #                  xlab="",ylab="SSB")
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-SSBdMSST.pdf"))
  #   plotBootSeries(dataBoot=D.boot.SSB.msst[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$SSB.msst,
  #                  xlab="",ylab="SSB/MSST")
  #   abline(h=1,lwd=2)
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-SSBdSSBmsy.pdf"))
  #   plotBootSeries(dataBoot=D.boot.SSB.SSBmsy[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$SSB.SSBmsy,
  #                  xlab="",ylab="SSB/SSBmsy")
  #   abline(h=1,lwd=2)
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-recruits.pdf"))
  #   plotBootSeries(dataBoot=D.boot.recruits[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$recruits,
  #                  xlab="",ylab="Recruits (numbers)")
  #   dev.off()
  #
  #   # Benchmarks time series, three panel plot (SSBdMSST, SSBdSSBmsy, and FdFmsy)
  #   pdf(file.path(dir_figs,"Fig-MCB-status-ts.pdf"), width=8,height=10)
  #   par(mfrow=c(2,1),mar=c(1.1,2,1,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
  #   plotBootSeries(dataBoot=D.boot.SSB.msst[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$SSB.msst,
  #                  xlab="",ylab="SSB/MSST")
  #   abline(h=1,lwd=2)
  #   # plotBootSeries(dataBoot=D.boot.SSB.SSBmsy[,paste(yr.plot)],runsToKeep = simID2,
  #   #                x.ref=base.t.series$year,y.ref=base.t.series$SSB.SSBmsy,
  #   #                xlab="",ylab="SSB/SSBmsy")
  #   # abline(h=1,lwd=2)
  #   plotBootSeries(dataBoot=D.boot.F.Fmsy[,paste(yr.plot)],runsToKeep = simID2,
  #                  x.ref=base.t.series$year,y.ref=base.t.series$F.Fmsy,
  #                  xlab="",ylab="F/Fmsy")
  #   abline(h=1,lwd=2)
  #   dev.off()
  #
  #   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   ####     plotASeries      ####
  #   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   pdf(file.path(dir_figs,"Fig-MCB-Mage.pdf"))
  #   plotBootSeries(dataBoot=D.boot.M,runsToKeep = simID2,
  #                  x.ref=spp$a.series$age,y.ref=spp$a.series$M,
  #                  xlab="age",ylab="M")
  #   dev.off()
  #
  #   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   ####     plotDensity      ####
  #   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   pdf(file.path(dir_figs,"Fig-MCB-status.pdf"),width=8,height=10)
  #   #par(mfrow=c(3,1),mar=c(2,2,1,1),mgp=c(1.2,.3,0))
  #   par(mfrow=c(2,1),mar=c(2,2,1,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
  #   # plotBootDensity(data=D.boot.parms[simID2,],par.name="SSBend.SSBmsy",xlab=paste("SSB(",max(yr.plot),")/SSBmsy",sep=""))
  #   # abline(v=spp$parms$SSBend.SSBmsy,lty=1,lwd=2)
  #   # abline(v=median(D.boot.parms[simID2,"SSBend.SSBmsy"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="SSBend.MSST",xlab=paste("SSB(",max(yr.plot),")/MSST",sep=""),
  #                   xlim=c(0.01,6),ylim=c(0,0.9))
  #   abline(v=spp$parms$SSBend.MSST,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"SSBend.MSST"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="Fend.Fmsy.mean",xlab=paste("F(",max(yr.plot)-2,"-",max(yr.plot),")/Fmsy",sep=""),
  #                   xlim=c(0.01,6),ylim=c(0,0.9))
  #   abline(v=spp$parms$Fend.Fmsy.mean,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"Fend.Fmsy.mean"]),lty=2,lwd=2)
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-benchmarks.pdf"))
  #   par(mfrow=c(2,2),mar=c(3,2,1,1),mgp=c(1.5,.3,0),tck=-0.01,lend="butt")
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="Fmsy",
  #                   xlab=expression(Fmsy~(yr^-1)))
  #   abline(v=spp$parms$Fmsy,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"Fmsy"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="SSBmsy",
  #                   xlab="SSBmsy (mt)")
  #   abline(v=spp$parms$SSBmsy,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"SSBmsy"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="msy.klb",
  #                   xlab="MSY (1000 lb)")
  #   abline(v=spp$parms$msy.klb,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"msy.klb"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="Bmsy",
  #                   xlab="Bmsy (mt)")
  #   abline(v=spp$parms$Bmsy,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"Bmsy"]),lty=2,lwd=2)
  #   dev.off()
  #
  #   pdf(file.path(dir_figs,"Fig-MCB-SR-parms.pdf"))
  #   par(mfrow=c(2,2),mar=c(3,2,1,1),mgp=c(1.2,.3,0),tck=-0.01,lend="butt")
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="R0")
  #   abline(v=spp$parms$R0,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"R0"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="BH.steep",
  #                   xlab="steepness")
  #   abline(v=spp$parms$BH.steep,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"BH.steep"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="BH.Phi0",
  #                   xlab="Unfished spawning biomass per recruit") # AKA sprF0
  #   abline(v=spp$parms$BH.Phi0,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"BH.Phi0"]),lty=2,lwd=2)
  #
  #   plotBootDensity(data=D.boot.parms[simID2,],par.name="R.sigma.par",
  #                   xlab="SD of log recruitment residuals")
  #   abline(v=spp$parms$R.sigma.par,lty=1,lwd=2)
  #   abline(v=median(D.boot.parms[simID2,"R.sigma.par"]),lty=2,lwd=2)
  #   dev.off()
  #
  #   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   ####      plotPhase      ####
  #   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   pdf(file.path(dir_figs,"Fig-MCB-phase.pdf"))
  #   par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
  #
  #   # Fig-MCB-phaseMSST
  #   plotBootPhase(x.Boot=D.boot.parms[simID2,"Fend.Fmsy.mean"],y.Boot = D.boot.parms[simID2,"SSBend.MSST"],
  #                 x.ref = spp$parms$Fend.Fmsy.mean, y.ref = spp$parms$SSBend.MSST,
  #                 xlab=paste("F(",max(yr.plot)-2,"-",max(yr.plot),")/Fmsy",sep=""),
  #                 ylab=paste("SSB(",max(yr.plot),")/MSST",sep=""),
  #                 text.y.adj=c(0.2,0,0.2,0))
  #
  #   # Fig-MCB-phaseSSB
  #   # plotBootPhase(x.Boot=D.boot.parms[simID2,"Fend.Fmsy.mean"],y.Boot = D.boot.parms[simID2,"SSBend.SSBmsy"],
  #   #               x.ref = spp$parms$Fend.Fmsy.mean, y.ref = spp$parms$SSBend.SSBmsy,
  #   #               xlab=paste("F(",max(yr.plot)-2,"-",max(yr.plot),")/Fmsy",sep=""),
  #   #               ylab=paste("SSB(",max(yr.plot),")/SSBmsy",sep=""),
  #   #               text.y.adj=c(0.2,0,0.2,0))
  #   dev.off()
  #
  # }



}
