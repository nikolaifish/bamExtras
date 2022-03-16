#' Convert lognormal cv to lognormal sd
#'
#' @param cv lognormal cv
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' set.seed(1234)
#' x <- rnorm(1000,mean=100,sd=3)
#' m <- mean(log(x))
#' s <- sd(log(x))
#' lnorm_sd(mu=m,sigma=s)

lnorm_cv_to_sd <- function(cv){
  (log(1.0+cv^2))^0.5
}
