#' Complete a length or age composition data frame so that all have desired dimensions
#'
#' Complete a length or age composition data frame so that all have desired dimensions. This is a modified version of comp_complete. The two functions should probably be combined with options to choose between sub-functions.
#' @param comp_mat: comp matrices, where column names are sample size columns (e.g. n.fish, n.trips) as well as bin values (e.g. length or age), and rows are unique observations (e.g. years)
#' @param n_colnames: column names for sample size columns
#' @param val_rownames:  row names you want to have included in each comp matrix
#' @param val_colnames: column names you want to have included in each comp matrix (excluding sample size columns)
#' @param minusGroup: logical. When truncating bin range, should values below the smallest bin be summed into the smallest bin? The default is to simply truncate.
#' @param plusGroup: logical. When truncating bin range, should values above the largest bin be summed into the largest bin? The default is to simply truncate.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export

comp_complete_bam <- function(comp_mat,
                              n_colnames=NULL,
                              val_rownames=NULL,
                              val_colnames=NULL,
                              minusGroup=FALSE,
                              plusGroup=FALSE){
  colnames(comp_mat) <- tolower(colnames(comp_mat))
  if(is.null(n_colnames)){
    n_colnames <- colnames(comp_mat)[grepl(pattern = "^n",x=colnames(comp_mat))]
  }
  val_colnames_obs <- colnames(comp_mat)[!colnames(comp_mat)%in%c(n_colnames,"year","yr")]

  if(is.null(val_rownames)){
    val_rownames <- paste(min(as.numeric(rownames(comp_mat))):max(as.numeric(rownames((comp_mat)))))
  }

  if(is.null(val_colnames)){
    val_colnames <- paste(min(as.numeric(val_colnames_obs)):max(as.numeric(val_colnames_obs)))
  }

  D_n <- comp_mat[,n_colnames,drop=FALSE]

  M_cmp_obs <- as.matrix(comp_mat[,val_colnames_obs])

  # Add minus group if specified, and if data are going to be truncated
  # (by default, function simply excludes values below minimum val_colnames)
  if(minusGroup&min(as.numeric(val_colnames))>min(as.numeric(val_colnames_obs))){
    M_cmp_obs_minus <- M_cmp_obs[,as.numeric(val_colnames_obs) <= min(as.numeric(val_colnames)),drop=FALSE]
    M_cmp_obs[,paste(min(as.numeric(val_colnames)))] <-  rowSums(M_cmp_obs_minus) # Add plus group values
    M_cmp_obs[,which(as.numeric(val_colnames_obs) < min(as.numeric(val_colnames)))] <- 0 # Replace values in bins above with zeros
  }

  # Add plus group if specified, and if data are going to be truncated
  # (by default, function simply excludes values above maximum val_colnames)
  if(plusGroup&max(as.numeric(val_colnames))<max(as.numeric(val_colnames_obs))){
    M_cmp_obs_plus <- M_cmp_obs[,as.numeric(val_colnames_obs) >= max(as.numeric(val_colnames)),drop=FALSE]
    M_cmp_obs[,paste(max(as.numeric(val_colnames)))] <-  rowSums(M_cmp_obs_plus) # Add plus group values
    M_cmp_obs[,which(as.numeric(val_colnames_obs) > max(as.numeric(val_colnames)))] <- 0 # Replace values in bins above with zeros
  }

  M_cmp_out <- matrix(NA,nrow=length(val_rownames),ncol=length(val_colnames),dimnames=list(val_rownames,val_colnames))

  # Add observed values to the appropriate cells in the output matrix
  for(obs_rowName_i in rownames(M_cmp_obs)[which(rownames(M_cmp_obs)%in%rownames(M_cmp_out))]){
    obs_colName_i <- colnames(M_cmp_obs)[which(colnames(M_cmp_obs)%in%colnames(M_cmp_out))]
    M_cmp_out[obs_rowName_i,] <- rep(0,ncol(M_cmp_out)) # First fill with zeros
    M_cmp_out[obs_rowName_i,obs_colName_i] <- M_cmp_obs[obs_rowName_i,obs_colName_i] # Now fill with observed data
  }

  M_cmp_out <- t(apply(M_cmp_out,MARGIN=1,function(x){x/sum(x)})) # rescale so that each row sums to 1

  # Add sample size columns
  D_n2 <- D_n[rownames(M_cmp_out),,drop=FALSE]
  rownames(D_n2) <- rownames(M_cmp_out)
  cbind(D_n2,M_cmp_out)

  comp_mat_out <- as.data.frame(cbind(D_n2,M_cmp_out))
  return(comp_mat_out)
}
