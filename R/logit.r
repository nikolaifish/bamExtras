#' logit
#'
#' simple logit function.
#' @param x numeric vector (e.g. probabilities, 0 < x < 1)
#' @note This is really the same as \code{qlogis()} but it's nice to be able to just call \code{logit()}.
#' Returns Inf when x = 0 or x = 1 and NaN for x < 0 or x > 1.
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' x <- 1:20
#' y <- bamExtras::lgs(x=x)
#' z <- bamExtras::logit(x=y)
#' z2 <- stats::qlogis(y)
#' plot(x,z)
#' points(x,z2,type="l",col="blue")
#' }

logit <- function(x){
  log(x/(1-x))
}
