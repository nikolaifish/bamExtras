#' Generating bootstrap vectors of with lognormal error
#'
#' This function was designed to be used with time series of indices of abundance or removals. Based on code from previous BAM
#' Monte-Carlo Bootstrap code
#' @param x observed vector
#' @param cv observed cv for x
#' @param bootN number of bootstrap replicates. Allows multiple bootstrapped vectors to be calculated without looping
#' @param standardize divide resulting vectors by their own mean so that they have a mean of 1
#' @param digits Number of significant digits to round output to. If NULL, output is not rounded.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky, Kyle Shertzer, Erik Williams
#' @export
#' @examples
#' \dontrun{
#' rdat <- rdat_VermilionSnapper
#' year <- rdat$t.series$year
#' U <- setNames(rdat$t.series$U.cHL.ob,year)
#' cv_U <- setNames(rdat$t.series$cv.U.cHL,year)
#' nsim <- 100
#' U_boot <- lnorm_vector_boot(U,cv_U,nsim)
#' dimnames(U_boot) <- list("year"=year,"sim"=sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim))
#' plot_boot_vec(t(U_boot))
#' points(as.numeric(year),U,type="o",col="blue")
#' }

lnorm_vector_boot <- function(x,
                              cv,
                              bootN=1,
                              standardize=FALSE,
                              digits=NULL){
  cc <- which(complete.cases(x,cv))
  xboot_obs <- rep(x[cc],bootN)
  M_xBoot <- matrix(NA,nrow=length(x),ncol=bootN) # Create matrix from xboot

  n <- length(xboot_obs)
  sd <- (log(1.0+cv[cc]^2))^0.5 # Calculate vector of lognormal sd from lognormal CV
  lnormResid <- rlnorm(n=n, mean=-(sd^2)/2, sd=sd) # Generate vector of lognormal residuals (proportions)
  xBoot <- xboot_obs*lnormResid # Multiply observed index vector by bootstrap lognormal residuals
  xBoot[which(xBoot<=0)]=xboot_obs[which(xBoot<=0)] # If any of the bootstrap values are <=0, replace them with the observed values
  M_xBoot_cc <- matrix(xBoot,nrow=length(cc),ncol=bootN) # Create matrix from xboot
  M_xBoot[cc,] <- M_xBoot_cc # Add values to matrix including NA rows

  if(standardize){
    #xBoot <- xBoot/mean(xBoot) # Standardize bootstrap vector
    M_xBoot <- apply(M_xBoot,MARGIN=2,function(y){y/mean(y,na.rm=TRUE)}) # Standardize each column
  }
  if(is.numeric(digits)){
    M_xBoot <- signif(M_xBoot,digits=digits) # Round bootstrap vector
  }
  if(!is.null(names(x))){rownames(M_xBoot) <- names(x)}
  return(M_xBoot)
}
