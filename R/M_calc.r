#' Compute constant natural mortality estimates from life history parameters
#' 
#' @param type Type of M calculation to use. If type="all" the function tries all calculations
#'  and returns a vector of all results. NA values indicate that a particular calculation couldn't
#'  be computed, probably because a necessary life history parameter was not provided.
#'  Possible values include four calculations specified in Hoenig (1983; Table 1;
#'  "Hng_Mollusks", "Hng_Fish", "Hng_Cetaceans", "Hng_All") and six calculations
#'   from Then et al. (2014; Table 3; "Thn_1partmax", "Thn_Hnglm", "Thn_Hngnls", 
#'   "Thn_1parK", "Thn_2parK", "Thn_Plynls").
#' @param Linf length infinity in cm (von Bertalanffy growth equation parameter)
#' @param K growth coefficient (von Bertalanffy growth equation parameter)
#' @param t0 age (time) at length zero (von Bertalanffy growth equation parameter)
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @references Hoenig, J. M. 1983. Empirical use of longevity data to estimate mortality-rates. Fishery Bulletin 81:898-903.
#' @references Then, A. Y., J. M. Hoenig, N. G. Hall, and D. A. Hewitt. 2014. Evaluating the predictive performance of empirical estimators of natural mortality rate using information on over 200 fish species. ICES Journal of Marine Science.
#' @export
#' @examples
#' \dontrun{
#' # Return all possible estimates based on tmax
#' M_calc(tmax=11)
#' # Return all possible estimates based on tmax or K
#' M_calc(tmax=11,K=0.173)
#' # Return all possible estimates based on tmax, K, or Linf (in centimeters)
#' M_calc(tmax=11,Linf=50.2,K=0.173)
#' # Return just the Hoenig_nls estimate published in Then et al. (2014) based on tmax
#' M_calc(type="Thn_Hngnls",tmax=11)
#' }

M_calc <- function(type="all",tmax,Linf,K) {
  
  errorFcn <-  function(msg){return(NA)}
  
  # Try to run all of these values. Return NA if error.
  # From Table 1 of Hoenig, J. M. 1983. Empirical use of longevity data to estimate mortality-rates. Fishery Bulletin 81:898-903.
  M <- list(
    Hng_Mollusks =  tryCatch({3.42123*tmax^-0.832}, error=errorFcn),   # originally written as exp(1.23  - 0.832*log(tmax))
    Hng_Fish = tryCatch({4.30596*tmax^-1.01}, error=errorFcn),        # originally written as exp(1.46  - 1.01* log(tmax))
    Hng_Cetaceans = tryCatch({2.562543*tmax^-0.873}, error=errorFcn), # originally written as exp(0.941 - 0.873*log(tmax))
    Hng_All = tryCatch({4.220696*tmax^-0.982}, error=errorFcn),       # originally written as exp(1.44  - 0.982*log(tmax))
    # From Table 3 of  Then, A. Y., J. M. Hoenig, N. G. Hall, and D. A. Hewitt. 2014. Evaluating the predictive performance of empirical estimators of natural mortality rate using information on over 200 fish species. ICES Journal of Marine Science.
    Thn_1partmax = tryCatch({5.109/tmax}, error=errorFcn),
    Thn_Hnglm = tryCatch({5.5678*tmax^-1.01}, error=errorFcn), # originally written as exp(1.717  - 1.01* log(tmax))
    Thn_Hngnls = tryCatch({4.899*tmax^-0.916}, error=errorFcn),
    Thn_1parK = tryCatch({1.692*K}, error=errorFcn),
    Thn_2parK = tryCatch({0.098 + 1.55*K}, error=errorFcn),
    Thn_Plynls = tryCatch({4.118 * (K^0.73) * (Linf^-0.33)}, error=errorFcn)
  )
  
  if(type=="all"){
    M
  }else{
    M[[type]]
  }
}