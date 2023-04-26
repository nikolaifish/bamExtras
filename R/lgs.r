#' lgs
#'
#' logistic function
#' @param x numeric vector
#' @param a slope parameter
#' @param b point of inflection (e.g. a50)
#' @param par named vector of parameters a and b. If this is specified it will
#' @param type which form of the logistic function to use: "slope_intercept" for
#' the slope-intercept form, "slope_inflection" for the slope-point-of-inflection form.
#' For both forms, a is the slope parameter.
#' override specifications of individual parameters
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' x <- 1:40
#' y <- lgs(x=x, a=1, b=5)
#' plot(x,y)
#'
#' # Alternative way to specify parameters
#' y <- lgs(x=x, par=c(x=x, a=1, b=5))
#' plot(x,y)
#' }

lgs <- function(x, a=1, b=0,par=NULL,type="slope_inflection"){
  if(!is.null(par)){
    a <- par["a"]
    b <- par["b"]
  }
  switch(type,
         slope_inflection = 1/(1+exp(-a*(x-b))),
         slope_intercept = exp(a*x+b)/(1+exp(a*x+b))
         )
}
