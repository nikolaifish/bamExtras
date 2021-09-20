#' Exponential decay
#'
#' Compute numbers of individuals alive at time in a population experiencing exponential decay, allowing age-varying Z. Can be used to compute numbers at age at the beginning of the year.
#' @param age ages
#' @param Z instantaneous mortality rate
#' @param N0 number of individuals at age zero at the beginning of the year
#' @param plus.group Should the function include a plus group? logical
#' @keywords bam stock assessment fisheries population dynamics
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' expDecay(age=1:10,Z=0.2,N0=100)
#' expDecay(age=1:10,Z=0.2,N0=100,plus.group = TRUE)
#' expDecay(age=c(2,5,10),Z=0.2,N0=100,plus.group = TRUE)
#' }
expDecay <- function(age,Z,N0=1,plus.group=FALSE){
  ac <- 1:(max(age)+1) # Vector of all age classes
  ages <- ac-1
  ac.n <- length(ac)
  N.a <-  setNames(rep(0,ac.n),ages)    # Numbers at age
  if(length(Z)==1){
    Z <- rep(Z,ac.n)
  }
  a.steps <-  ac.n/(max(ac)-min(ac)+1) # Number of age steps per age class
  # Numbers at age (based on Gabriel et al. 1989)
  for (a in ac){
    N.a[a] <- local({
      if(a==1)   {return(N0)}
      if(a>1)    {
        if(plus.group&a==ac.n){return(N.a[ac.n-1]*exp(-Z[ac.n-1]/a.steps)/   # Plus group
                                        (1-exp(-Z[ac.n-1]/a.steps)))
        }else{return(N.a[a-1]* exp(-Z[a-1]/a.steps))}
      }
    })
  }
  return(N.a[paste(age)])
}
