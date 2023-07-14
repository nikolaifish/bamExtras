#' F_calc
#'
#' Calculations of stock dependent quantities over a range of fishing mortalities (F). This code calculates spawning potential ratio (SPR), and then maximum sustainable yield (MSY) using equilibrium methods incorporating a stock-recruit relationship. Many of the calculations, particularly the MSY calculations, are based on code
#' written by Kyle Shertzer in a very similar function.
#' @param ac age classes. numeric vector
#' @param h Beverton-Holt steepness parameter
#' @param R0 Beverton-Holt R0 parameter. Numbers of fish at age-a (often age-0 or age-1).
#' @param sigma lognormal bias correction -- exp(sigma^2/2)
#' @param M Natural mortality rate
#' @param Fmax Maximum fishing mortality rate
#' @param F_n Number of fishing mortality rates to try
#' @param sel Selectivity at age
#' @param PS Proportion of fish to include in stock calculation, at age (e.g. mature females at age)
#' @param W weight at age
#' @param ep egg production proxy at age (e.g. weight, fecundity)
#' @param Px reference proportion of unfished spawning potential ratio (e.g. Px=.40)
#' @param plot_digits number of significant digits to show in plot
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer and Nikolai Klibansky
#' @export
#' @examples
#' rdat <- rdat_VermilionSnapper
#' pr <- rdat$parms
#' as <- rdat$a.series
#' F_calc(ac=as$age, h=pr$BH.steep, R0=pr$R0, M=as$M, sel=rdat$sel.age$sel.v.wgted.tot,
#'        PS=as$mat.female, W=as$weight, ep=as$wgt.klb,
#'        plots=TRUE)
#'
#' # When there is no SRR, the function will still compute SPR without computing MSY-based equibrium values
#' rdat <- rdat_RedSnapper
#' pr <- rdat$parms
#' as <- rdat$a.series
#' F_calc(ac=as$age, M=as$M, sel=rdat$sel.age$sel.v.wgted.tot,
#'        PS=as$mat.female, W=as$weight, ep=as$wgt.klb,
#'        plots=TRUE)

F_calc <-  function(ac, h=NULL, R0=NULL, sigma=0, M, Fmax=1, F_n=101, sel, PS, W, ep, Px=.40,
                 plots=FALSE, plus_group=TRUE, plot_digits=3){

F <- seq(0,Fmax,length=F_n)
ac_n=length(ac) # Number of age classes
a_steps <- ac_n/(max(ac)-min(ac)+1)  # Number of age steps per age

spr <- SPR <- R_eq <- S_eq <- B_eq <- Y_eq <- E_eq <- rep(0,length(F)) # Initialize storage vectors

is_SRR <- !is.null(h)&!is.null(R0) # Are stock recruit relationship parameters provided?
if(!is_SRR){
  warning("Either h or R0 is null. Equilibrium calculations will not be computed.")
}

# FOR EACH LEVEL OF F..
for (fi in 1:F_n) {
  # NON-EQUILIBRIUM CALCULATIONS
  # ALL FISH
  N_a <- exp_decay(age=ac,Z=M+sel*F[fi],N0=1,plus_group=plus_group)

  # Numbers at age in the unfished population (when F=0)
  if(F[fi]==0){
    N_a_F0 <- data.frame("age_class"=ac, "N_F0"=N_a)
  }

  # Number of fish in the spawning stock (e.g. number of mature females) at age
  N_S_a <- N_a*PS
  # Spawners per recruit
  spr[fi] <- sum(N_S_a*ep)
  # Unfished spawning biomass per recruit
  phi0 <- spr[1]
  # Spawning potential ratio
  SPR[fi] <- spr[fi]/spr[1]

  # EQUILIBRIUM CALCULATIONS
  if(is_SRR){
    # Incorporate the SRR into the following calculations
    # Recruitment (equilibrium recruitment based on SSB/R at F, and R_eq)
    R_eq[fi] <- local({BC<-exp(sigma^2/2.0)              #multiplicative bias correction
    R_eq <- (R0*(spr[fi]*4*h*BC-phi0*(1-h)))/(spr[fi]*(5*h-1))
    ifelse(R_eq<1e-7, 1e-7, R_eq)   # I think this just keeps R_eq from getting too close to zero.
    })
    # Numbers at age
    N_a_eq <- R_eq[fi]*N_a
    # Number of fish in the spawning stock (e.g. number of mature females) at age
    N_S_a_eq <- N_a_eq*PS
    # Spawning stock at F
    S_eq[fi] <- sum(N_S_a_eq*ep)
    # Total biomass at F
    B_eq[fi] <- sum(N_a_eq*W)
    # Catch (numbers) at age (Gabriel et al. 1989)
    C_a_eq <- ((sel*F[fi])/(M+sel*F[fi])) * (1-exp(-(M+sel*F[fi])/a_steps)) * N_a_eq
    # Yield (weight)
    Y_eq[fi] <- sum(C_a_eq*W)
    # Exploitation rate (total catch/number of fish)
    E_eq[fi] <- sum(C_a_eq)/sum(N_a_eq)
  }
}

  # Reference points
    SPR_Px<-Px  # SPR_Px
    #SPR_Px<<-which(abs(SPR-Px)==min(abs(SPR-Px)))  # Index value corresponding to SPR_Px
    # Calculate close approximation of F at SPR_Px, by interpolating between the nearest values
        x1_ix <- which(SPR==max(SPR[SPR<=SPR_Px]))
        x2_ix <- which(SPR==min(SPR[SPR>SPR_Px]))
        x1 <- SPR[x1_ix]  # Value just below SPR_Px
        x2 <- SPR[x2_ix]   # Value just above SPR_Px

        # F_Px
          F_Px <- local({
          y1 <- F[x1_ix]      # F at x1
          y2 <- F[x2_ix]      # F at x2
          B1 <- (y2-y1)/(x2-x1)            # Slope of the line connecting the adjacent points
          B0 <- y1-B1*x1                   # Intercept of the line connecting the adjacent points
          return(B0+B1*SPR_Px)          # F_Px calculated through linear interpolation
          })


    #F_Px <-    F[SPR_Px_ix]    # F at SPR reference value
    #SPR_Px <-  SPR[SPR_Px_ix]  # SPR at SPR reference value
          if(is_SRR){
            MSY <-      max(Y_eq)          # maximum sustainable yield
            MSY_ix <-   which(Y_eq==MSY) # MSY F-index
            FMSY <-    F[MSY_ix]           # fishing rate at MSY
            sprMSY <-  spr[MSY_ix]         # spawners per recruit at MSY
            SPRMSY <-   sprMSY/phi0        # spawning potential ratio at MSY
            RMSY <-    R_eq[MSY_ix]        # equilibrium recruitment at MSY
            SMSY <-    S_eq[MSY_ix]        # spawning stock at MSY (e.g. SSB, fecundity)
            BMSY <-    B_eq[MSY_ix]        # total biomass (male and female) at MSY
            EMSY <-    E_eq[MSY_ix]        # exploitation rate at MSY
          }else{
            MSY <- FMSY <- sprMSY <- SPRMSY <- RMSY <- SMSY <- BMSY <- EMSY <- NA
          }

ref_points <- if(sum(ep)>0){data.frame(F_Px, SPR_Px, MSY, FMSY, sprMSY, SPRMSY, RMSY, SMSY, BMSY, EMSY)
           }else{data.frame("F_Px"=NA,"SPR_Px"=SPR_Px,"MSY"=NA,"FMSY"=NA,"sprMSY"=NA,"SPRMSY"=NA,
                            "RMSY"=NA,"SMSY"=NA,"BMSY"=NA,"EMSY"=NA)}
data <- data.frame("F"=F, "spr"=spr, "SPR"=SPR,
                "R_eq"=R_eq, "S_eq"=S_eq, "B_eq"=B_eq, "Y_eq"=Y_eq, "E_eq"=E_eq)
if(plots){
  par(mar=c(2,2,2,1),mfrow=c(3,3),mgp=c(1.1,0.1,0),tck=0.01)
  plotNK <- function(...,rp=NA){
    plot(type="l",lwd=2,...)
    if(!is.na(rp)){points(FMSY,rp,type="p",col="blue",pch=16)}
  }

  plotNK(F,spr,  main="Non-Equilibrium", rp=sprMSY)
  plotNK(F,SPR,  main="Non-Equilibrium", rp=SPRMSY)
  usr <- par("usr")
  #text(FMSY,SPRMSY, labels=bquote(SPR[MSY] == .(signif(SPRMSY,plot_digits))),pos=4)
  text(FMSY,SPRMSY, labels=bquote(list(F[MSY] == .(signif(FMSY,plot_digits)),SPR[MSY] == .(signif(SPRMSY,plot_digits)))),pos=4)
  points(F_Px,SPR_Px,type="p",col="red",pch=16)
  text(F_Px,SPR_Px, labels=bquote(list(F[.(Px)] == .(signif(F_Px,plot_digits)),SPR[.(Px)] == .(signif(SPR_Px,plot_digits)))),pos=4)
  #text(F_Px,usr[3]+diff(usr[3:4])*.05, labels=bquote(F[.(Px)] == .(signif(F_Px,plot_digits))),pos=4)
  if(is_SRR){
    plotNK(F,Y_eq, main="Equilibrium",     rp=MSY)
    text(FMSY,MSY, labels=bquote(MSY == .(signif(MSY,plot_digits))),pos=4)
    plotNK(F,R_eq, main="Equilibrium",     rp=RMSY)
    text(FMSY,RMSY, labels=bquote(R[MSY] == .(signif(RMSY,plot_digits))),pos=4)
    plotNK(F,S_eq, main="Equilibrium",     rp=SMSY)
    text(FMSY,SMSY, labels=bquote(S[MSY] == .(signif(SMSY,plot_digits))),pos=4)
    plotNK(F,B_eq, main="Equilibrium",     rp=BMSY)
    text(FMSY,BMSY, labels=bquote(B[MSY] == .(signif(BMSY,plot_digits))),pos=4)
    plotNK(F,E_eq, main="Equilibrium",     rp=EMSY)
    text(FMSY,EMSY, labels=bquote(E[MSY] == .(signif(EMSY,plot_digits))),pos=4)
  }
}

invisible(list("N_a_F0"=N_a_F0, "RefPts"=ref_points, "Data"=data))
}
