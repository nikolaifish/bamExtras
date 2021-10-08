#' Computes mean excluding NA values
#'
#' Version of mean() where default na.rm=TRUE
#' @param x numeric vector
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' x <- c(1,2,10,NA)
#' mean(x)
#' mean(x,na.rm=TRUE)
#' mean_noNA(x)

mean_noNA <- function(x){mean(x,na.rm=TRUE)}
