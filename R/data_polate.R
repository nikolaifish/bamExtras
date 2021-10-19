#' Linearly extrapolate/interpolate variables in a numeric matrix
#'
#' Linearly extrapolate/interpolate variables in a numeric matrix `data`. This function uses `Hmisc::approxExtrap()` and `base::apply()` to extrapolate variables y_i to y_n in a data frame based on univariate relationships with a variable x (`x`; which is another variables in `data`) to a new variable `xout`.
#'
#' @param data a numeric matrix or data frame that can be coerced to a numeric matrix with column names containing the x variable and y_i variables to extrapolate
#' @param xname the name of the column containing the x variable. If `is.null(xname)` the first column will be used.
#' @param ynames the names of the columns containing the y variables. If `is.null(ynames)` all column names other than `xname` are used.
#' @param xout a numeric vector of values to x values to extrapolate y_i to. If `if(is.null(xout)){xout <- x}`, and the function doesn't extrapolate anything.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @note Users should be cautious about conducting linear extrapolation for relationships that are known to be non-linear over a large range, because they will undoubtedly deviate from the true non-linear relationships. Linear interpolation will also deviate, though the risks will tend to be less severe especially if the resolution of the original data is relatively fine.
#' @export
#' @examples
#' \dontrun{
#' mydata <- rdat_RedGrouper$a.series
#'  out <- data_polate(mydata,xout=seq(min(mydata$age),max(mydata$age),by=0.5))
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

data_polate <- function(data,xname=NULL,ynames=NULL,xout=NULL){
  data <- as.matrix(data)
  if(is.null(xname)){xname <- colnames(data)[1]}
  if(is.null(ynames)){ynames <- colnames(data)[colnames(data)!=xname]}
  x <- data[,xname,drop=FALSE]
  if(is.null(xout)){xout <- x}
  xout <- as.matrix(xout)
  colnames(xout) <- colnames(x)

  yout <- apply(data[,c(ynames)],2,function(y){Hmisc::approxExtrap(x=x,y=y,xout=xout)$y})
  out <- cbind(xout,yout)

  return(out)
}

