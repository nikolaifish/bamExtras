#' Convert Bayesian Information Criterion (BIC) to likelihood (lk) value
#'
#' Converts Bayesian Information Criterion (BIC) to likelihood (lk) value
#' @param bic Bayesian Information Criterion (BIC) value
#' @param k number of parameters estimated
#' @param n number of observations
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' bic2lk()

bic2lk <- function(bic,k,n){exp((k*log(n)-bic)/2)}

#' @rdname bic2lk
#' @export
bic_to_lk <- bic2lk
