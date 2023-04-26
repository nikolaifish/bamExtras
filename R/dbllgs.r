#' dbllgs
#'
#' A double logistic function, which is the product of two logistic functions lgs1 and lgs2.
#' When all parameters are positive, lgs1 is a function of parameters a1 and b1,
#' and increases with x, while lgs2 is a function of parameters a2 and b2, and decreases
#' with x.
#' @param x numeric vector
#' @param a1 slope parameter for lgs1
#' @param b1 point of inflection (e.g. a50) for lgs1
#' @param a2 slope parameter for lgs2
#' @param b2 point of inflection (e.g. a50) for lgs2
#' @param type output type. Return the double logistic ("dbl"), "lgs1", "lgs2",
#' or "all" to return a data frame of "x", "lgs1", "lgs2", and "dbl"
#' @param par named vector of parameters a1, b1, a2, b2. If this is specified it will override specifications of individual parameters
#' @author Kyle Shertzer and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Nice domed shape
#' x <- 1:40
#' y <- dbllgs(x=x, a1=1, b1=5, a2=1, b2=10)
#' plot(x,y)
#'
#' # Alternative way to specify parameters
#' y <- dbllgs(x=x, par=c(x=x, a1=1, b1=5, a2=1, b2=10))
#' plot(x,y)
#'
#' # Return lgs1
#' y <- dbllgs(x=x, a1=1, b1=5, a2=1, b2=10, type="lgs1")
#' plot(x,y)
#'
#' # Return lgs2
#' y <- dbllgs(x=x, a1=1, b1=5, a2=1, b2=10, type="lgs2")
#' plot(x,y)
#'
#' # Return all
#' dat <- dbllgs(x=x, a1=1, b1=5, a2=1, b2=10, type="all")
#' matplot(dat[,1],dat[,2:4],type="l",col=c("red","blue","purple"), lwd=2, lty=c(2,2,1))
#' }

dbllgs <- function(x, a1, b1, a2, b2,par=NULL,type="dbl"){
  if(!is.null(par)){
    a1 <- par["a1"]
    b1 <- par["b1"]
    a2 <- par["a2"]
    b2 <- par["b2"]
  }
  lgs1 <- 1/(1+exp(-a1*(x-b1)))
  lgs2 <- 1-1/(1+exp(-a2*(x-b2)))
  dbl <- lgs1*lgs2

  switch(type,
         lgs1=lgs1,
         lgs2=lgs2,
         dbl=dbl,
         all=data.frame(x=x,lgs1=lgs1,lgs2=lgs2,dbl=dbl)
         )
}
