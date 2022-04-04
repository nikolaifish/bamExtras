#' Coefficient of variation (cv) of lognormal distribution
#'
#' Coefficient of variation (cv) of lognormal distribution x, given the standard deviation (sigma) of log(x). Bolker 2008 page 137
#' @param sigma standard deviation of log(x)
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @references Bolker, B. M. 2008. Ecological Models and Data in R. Princeton University Press, Princeton, NJ. Page 137
#' @export
#' @examples
#' \dontrun{
#' set.seed(1234)
#' x <- rnorm(1000,mean=100,sd=3)
#' s <- sd(log(x))
#' lnorm_cv(sigma=s)
#' }
#'

lnorm_cv <- function(sigma){
  function(sigma){sqrt(exp(sigma^2)-1)}
}
