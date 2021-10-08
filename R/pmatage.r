#' Calculate proportion mature population from rdat. Either proportion female-at-age (for gonochoristic) or proportion combined proportion male and female maturity at age (for protogynous)
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param Mat_age1_max Limit maximum value of proportion mature of first age class (usually age-1). Models sometimes fail when maturity of first age class is too high (e.g. >0.5)
#' @param herm Is the species hermaphroditic? If "gonochoristic", use female maturity. If "protogynous", use a function of male and female maturity.
#' @param age Vector of ages to include
#' @keywords bam stock assessment fisheries DLMtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Identify proportion mature-at-age
#' pmat <- pmatage(rdat_RedPorgy)
#' plot(names(pmat),pmat)
#' }

pmatage <- function(
  rdat, Mat_age1_max = 0.49, herm, age=NULL
){
rdat <- standardize_rdat(rdat)
a.series <- rdat$a.series
if(is.null(age)){
age <- a.series$age
}
a.series <- a.series[paste(age),]

mat_female <- setNames(a.series$mat.female,age)

# If male maturity is provided use it, otherwise assume all males are mature.
if("mat.male"%in%names(a.series)){
  mat_male <- a.series$mat.male
}else{
  mat_male <- rep(1,nrow(a.series))
}
names(mat_male) <- age

prop_female <- setNames(a.series$prop.female,age)

# If gonochoristic use female maturity
if(herm=="gonochoristic"){
  pmat <- mat_female
}

# If protogynous, use a function of male and female maturity
if(herm=="protogynous"){
  pmat <- mat_female*prop_female+mat_male*(1-prop_female)
  }

names(pmat) <- age
pmat[1] <- min(pmat[1],Mat_age1_max)

return(list(prop_female=prop_female,mat_female=mat_female,mat_male=mat_male,pmat=pmat))
}
