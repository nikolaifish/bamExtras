#' Compute mean lengths by year above a certain minimum size (e.g. length at full selection).
#' @param CAL catch-at-length matrix, in numbers (rows=years, columns=length bins)
#' @param minL minimum length. Allows you to truncate the matrix to values above this value
#' @keywords bam stock assessment fisheries DLMtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#'
#' }
ML_LFS <- function(CAL,minL=0){
  CAL <-CAL[,as.numeric(colnames(CAL))>minL]
  year <- as.numeric(rownames(CAL))
  mlen <- apply(t(t(CAL)*as.numeric(colnames(CAL))),1,sum)/rowSums(CAL)
  ss <- rowSums(CAL)
  return(data.frame("year"=year,"mlen"=mlen,"ss"=ss))
}
