% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_lim.R
\name{data_lim}
\alias{data_lim}
\title{Truncate variables in a numeric matrix to specified limits}
\usage{
data_lim(data, xname = NULL, xlim = c(-Inf, Inf))
}
\arguments{
\item{data}{a numeric matrix or data frame that can be coerced to a numeric matrix}

\item{xname}{column names of variables to truncate}

\item{xlim}{numeric vector of length 2 with upper and lower limits to trim \code{xname} variables to.}
}
\description{
Truncate variables in a numeric matrix \code{data} to specified limits. This function was designed to be used after running \code{bamExtras::data_polate}, to limit life history variables to realistic values after linear extrapolation.
}
\examples{
\dontrun{
rdat <- standardize_rdat(rdat_RedGrouper)
mydata <- rdat$a.series
 mydata2 <- data_polate(mydata,xout=0:max(mydata$age))
 mydata3 <- data_lim(mydata2,xlim=c(0,Inf))
 out <-     data_lim(mydata3,xname=c("prop.female","prop.male","mat.female","mat.male"),xlim=c(0,1))

 par(mfrow=c(1,1))
for(i in colnames(out[,-1])){
 x <- mydata[,1]
 xout <- out[,1]
 yout_i <- out[,i]
 plot(xout,yout_i,ylab=i)
 points(x,mydata[,i],col="red",pch=16)
 readline(prompt="Press [enter] to continue")
}
}
}
\author{
Nikolai Klibansky
}
\keyword{assessment}
\keyword{bam}
\keyword{fisheries}
\keyword{stock}
