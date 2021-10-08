#' Standard deviation of lognormal distribution
#'
#' Standard deviation of lognormal distribution x, given the mean (mu) and standard deviation (sigma) of log(x). Bolker 2008 page 137
#' @param mu mean of log(x)
#' @param sigma standard deviation of log(x)
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @references Bolker, B. M. 2008. Ecological Models and Data in R. Princeton University Press, Princeton, NJ. Page 137
#' @export
#' @examples
#' set.seed(1234)
#' x <- rnorm(1000,mean=100,sd=3)
#' m <- mean(log(x))
#' s <- sd(log(x))
#' lnorm_sd(mu=m,sigma=s)

lnorm_sd <- function(mu,sigma){
  sqrt(exp(2*mu+sigma^2)*(exp(sigma^2)-1))
}
