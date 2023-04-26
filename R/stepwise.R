#' stepwise
#'
#' very simple step function used for estimating selectivity (i.e. independent proportions for each age or length). The function
#' structure allows it to be passed to \code{fn_sel} in \code{\link[bamExtras]{sel_comp}}.
#' @param x numeric vector (e.g. )
#' @param par named vector of parameters
#' @author Nikolai Klibansky
#' @export
#' @note This is really a dummy function to use as \code{fn_sel} in \code{\link[bamExtras]{sel_comp}}
#' @examples
#' \dontrun{
#' x <- 1:10
#' parms <- setNames(round(lgs(x,1,median(x)),2),paste0("a",x))
#' stepwise(x=x,par=parms)
#' }

stepwise <- function(x, par=NULL){
  par[1:length(x)]
}
