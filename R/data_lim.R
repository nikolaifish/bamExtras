#' Truncate variables in a numeric matrix to specified limits
#'
#' Truncate variables in a numeric matrix `data` to specified limits. This function was designed to be used after running `bamExtras::data_polate`, to limit life history variables to realistic values after linear extrapolation.
#'
#' @param data a numeric matrix or data frame that can be coerced to a numeric matrix
#' @param xname column names of variables to truncate
#' @param xlim numeric vector of length 2 with upper and lower limits to trim `xname` variables to.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' rdat <- standardize_rdat(rdat_RedGrouper)
#' mydata <- rdat$a.series
#'  mydata2 <- data_polate(mydata,xout=0:max(mydata$age))
#'  mydata3 <- data_lim(mydata2,xlim=c(0,Inf))
#'  out <-     data_lim(mydata3,xname=c("prop.female","prop.male","mat.female","mat.male"),xlim=c(0,1))
#'
#'  par(mfrow=c(1,1))
#' for(i in colnames(out[,-1])){
#'  x <- mydata[,1]
#'  xout <- out[,1]
#'  yout_i <- out[,i]
#'  plot(xout,yout_i,ylab=i)
#'  points(x,mydata[,i],col="red",pch=16)
#'  readline(prompt="Press [enter] to continue")
#' }
#' }

data_lim <- function(data,xname=NULL,xlim=c(-Inf,Inf)){
  data <- as.matrix(data)
  if(is.null(xname)){xname <- colnames(data)}
  xname <- xname[xname%in%colnames(data)] # Limit xname to variables that are actually in data
  data2 <- data[,xname,drop=FALSE]
  data3 <- apply(data2,2,function(x){pmin(pmax(x,xlim[1]),xlim[2])})
  data[,xname] <- data3
  return(data)
}
