#' agelenprob
#'
#' Build age-length probability matrix associated with a set of length comps
#' @param age vector of ages
#' @param len vector of lengths corresponding to age
#' @param binmid_lenc vector of lengths corresponding to midpoints of length composition bins. These do not need to match len.
#' @param cv_lenc cv associated with variation in length at age associated with lenc
#' @param plot logical. If TRUE, plot an image of the resulting probability matrix
#' @param label logical. If TRUE, add probabilities as text to image plot
#' @param label_min numeric. Set minimum probability to label in plot
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' rdat <- rdat_RedSnapper
#' agelenprob(
#'   age=rdat$a.series$age,
#'   len=rdat$a.series$length,
#'   binmid_lenc=as.numeric(colnames(rdat$comp.mats$lcomp.cHL.ob)),
#'   cv_lenc=rdat$parm.cons$len.cv.val.L[8],
#'   plot=TRUE)
#' }

agelenprob <- function(age,
                       len,
                       binmid_lenc,
                       cv_lenc,
                       plot=FALSE,
                       label=TRUE,
                       label_min=0.01){

  sd_lenc <- len*cv_lenc
  binmid_lenc <- as.numeric(colnames(lenc))  # bin mid point
  binwid <- local({
    a <- diff(binmid_lenc)
    if(length(unique(a))>1){
      warning("not all length bins in binmid_lenc have the same width")
    }
    median(a) # bin width
  })
  binmin_lenc <- binmid_lenc-binwid/2  # bin lo (minimum)
  lenprob <- matrix(NA,nrow=length(age),ncol=length(binmid_lenc),dimnames = list("age"=age,"lenbin"=binmid_lenc))
  lenprob[,1] <- pnorm((binmin_lenc[2]-len)/sd_lenc)
  for(i in 2:(ncol(lenprob)-1)){
    lenprob[,i] <- pnorm((binmin_lenc[i+1]-len)/sd_lenc)-pnorm((binmin_lenc[i]-len)/sd_lenc)
  }
  lenprob[,ncol(lenprob)] <- 1-rowSums(lenprob[,1:(ncol(lenprob)-1)])

  if(plot){
    lenbin <- binmid_lenc
    image(age,lenbin,lenprob)
    if(label){
      txt_x <- rep(age,length(lenbin))
      txt_y <- rep(lenbin,each=length(age))
      txt_z <- as.character(round(as.numeric(lenprob),2))
      txt_z[which(as.numeric(txt_z)<label_min)] <- ""
      text(txt_x,txt_y,labels=txt_z)
    }
  }
  return(lenprob)
}

