#' Simulate a random walk
#' 
#' Simple function to simulate a random walk
#' @param x independent variable values (e.g. vector of time steps)
#' @param y0 Initial value of y
#' @param sd Standard deviation of normal error distribution
#' @param mu multiplicative drift parameter
#' 
#' @export
#' @examples
#' \dontrun{
#' plot(random_walk(),type="l")
#' plot(random_walk(),random_walk(),type="l")
#' # Random walks can be pretty
#' x <- 1:1000
#' y1 <- random_walk(x)
#' y2 <- random_walk(x)
#' plot(y1,y2,type="l")
#' points(y1,y2,type="p",pch=16,col=rainbow(length(x)))
#' }

random_walk <- function(x=1:100, y0=100, sd=1, mu=0) {
  n <- length(x)
  err<-cumsum(rnorm(n=n, mean=0, 
                  sd=sd))
  y<-y0+x*mu+err
  return(y)
}