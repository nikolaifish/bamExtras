#' Compute mean squared errors (MSE) for indices
#'
#' Compute mean squared errors (MSE) for indices
#' @param x a BAM rdat object, i.e. x <- dget(rdat)
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' mse_calc(rdat_BlackSeaBass)
#' }
#'

mse_calc <- function(x){

  t.ser <- x$t.series
  t.ser[t.ser==-99999.00] <- NA
  t.ser.nm <- names(t.ser)

  comp.mats <- x$comp.mats
  comp.mats.nm <- names(comp.mats)

  # Indices
  U.nm.ob <- t.ser.nm[substr(t.ser.nm,1,2)=="U." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]
  U.nm.pr <- t.ser.nm[substr(t.ser.nm,1,2)=="U." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="pr"]
  U.nm <- gsub(".ob","",U.nm.ob)

  U.ob <- t.ser[,sort(U.nm.ob),drop=FALSE]
  U.pr <- t.ser[,sort(U.nm.pr),drop=FALSE]

  sq.errors <- (U.ob-U.pr)^2

  U.mse <- apply(sq.errors,2,function(x){mean(x,na.rm=TRUE)})
  names(U.mse) <- gsub(".ob","",names(U.mse))

  # Landings
  L.nm.ob <- t.ser.nm[substr(t.ser.nm,1,2)=="L." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]
  L.nm.pr <- t.ser.nm[substr(t.ser.nm,1,2)=="L." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="pr"]

  L.ob <- t.ser[,sort(L.nm.ob),drop=FALSE]
  L.pr <- t.ser[,sort(L.nm.pr),drop=FALSE]

  sq.errors <- (L.ob-L.pr)^2

  L.mse <- apply(sq.errors,2,function(x){mean(x,na.rm=TRUE)})
  names(L.mse) <- gsub(".ob","",names(L.mse))

  # comps
  comp.nm.ob <- comp.mats.nm[substr(comp.mats.nm,1,6)%in%c("acomp.","lcomp.") & substr(comp.mats.nm,nchar(comp.mats.nm)-1,nchar(comp.mats.nm))=="ob"]
  comp.nm.pr <- comp.mats.nm[substr(comp.mats.nm,1,6)%in%c("acomp.","lcomp.") & substr(comp.mats.nm,nchar(comp.mats.nm)-1,nchar(comp.mats.nm))=="pr"]

  comp.ob <- comp.mats[sort(comp.nm.ob)]
  comp.pr <- comp.mats[sort(comp.nm.pr)]

  comp.mse <- rep(NA,length(comp.nm.ob))
  names(comp.mse) <- gsub(".ob","",comp.nm.ob)

  for(i in 1:length(comp.nm.ob)){
    comp.mse[i] <- mean((comp.ob[[i]]-comp.pr[[i]])^2)
  }

  out <- c(U.mse,L.mse,comp.mse)
  names(out) <- paste("mse",names(out),sep=".")

  return(out)
}
