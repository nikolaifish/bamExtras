#' Complete each data frame in  a list of comp data frames so that all have the same dimensions
#'
#' @param comp_list list of composition matrices, where column names are bin values (e.g. length or age), and rows are unique observations (e.g. years)
#' @param comp_data_n data frame with annual sample sizes associated with each comp matrix
#' @param val_rownames Names of data rows for each comp matrix (e.g. years) desired in output.
#' @param val_colnames Names of data columns for each comp matrix (e.g. age bins, length bins) desired in output.
#' @param xBinByDataframe If true, x values are determined separately for each data frame (e.g. length comps may be binned separately by fleet)
#' @param n_tag tag(s) added to end of columns in comp_data_n which indicate sample sizes (e.g. ".n", ".nfish")
#' @param n_colname_new column name values to associate with n_tag in output
#' @param output_type Type of comps to output: "input"= format used in spreadsheets when submitting data, including ntrip and and nfish columns, "prop" = matrix of proportions, "nfish" = matrix of numbers of fish
#' @param valsToNA Values that should be replaced with NA (e.g. -99999)
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Combine age comp matrices for Red Porgy
#' cm <- rdat_RedPorgy$comp.mats
#' cma <- cm[grepl("^acomp.*ob$",names(cm))]
#'
#' cman <- rdat_RedPorgy$t.series[,gsub("ob$","n",names(cma))]
#' cmanfish <- rdat_RedPorgy$t.series[,gsub("ob$","nfish",names(cma))]
#'
#' # Output proportions at age
#' comp_complete(cma,output_type="prop")
#'
#' # Output numbers of fish at age
#' names(cma) <- gsub(".ob$","",names(cma))
#' comp_complete(cma,cmanfish,output_type="nfish",n_tag=".nfish")
#'
#' # Output age comps in data input format
#' comp_complete(cma,cbind(cman,cmanfish),output_type="input")
#' }
#'

comp_complete <- function(comp_list,comp_data_n,xBinByDataframe=FALSE,
                          n_tag= c(".n",".nfish"),
                          n_colname_new= c("ntrip","nfish"),
                          val_rownames=NULL,
                          val_colnames=NULL,
                          output_type="prop",
                          minusGroup=FALSE,
                          plusGroup=FALSE,
                          valsToNA=-99999){
  rownames_all <- sort(unique(unlist(lapply(comp_list,rownames))))
  colnames_all <- paste(sort(as.numeric(unique(unlist(lapply(comp_list,colnames))))))

  if(is.null(val_rownames)){
    val_rownames <- paste(min(as.numeric(rownames_all)):max(as.numeric(rownames_all)))
  }

  if(is.null(val_colnames)){
    val_colnames <- colnames_all
  }

  for (compName_i in names(comp_list)){
    M_i <- comp_list[[compName_i]]

    colnames_obs_i <- colnames(M_i)

    if(xBinByDataframe){
      val_colnames_i <- colnames(M_i)
    }else{
      val_colnames_i <- val_colnames
    }

    # Identify the specific row and column names missing comp matrix i
    row_missing <- val_rownames[!val_rownames%in%rownames(M_i)]
    col_missing <- val_colnames_i[!val_colnames_i%in%colnames(M_i)]

    ## Add missing rows
    # Create empty matrix of missing rows
    M_row_missing <- matrix(NA,nrow=length(row_missing),ncol=ncol(M_i),dimnames=list(row_missing,colnames(M_i)))
    # Actually add missing rows and columns to comp matrix i
    M_i2 <- rbind(M_i,M_row_missing)
    M_i2 <- local({M_i2 %>% rownames %>% order -> o; M_i2[o,]}) # Sort rows of M_i2

    ## Add missing columns
    # Create empty matrix of missing columns
    M_col_missing <- matrix(NA,nrow=nrow(M_i2),ncol=length(col_missing),dimnames=list(row.names(M_i2),col_missing))
    # Actually add missing rows and columns to comp matrix i
    M_i3 <- cbind(M_i2,M_col_missing)
    M_i3 <- local({M_i3 %>% colnames %>% as.numeric %>% order -> o; M_i3[,o]}) # Sort columns of M_i3

    # Add minus group if specified, and if data are going to be truncated
    # (by default, function simply excludes values below minimum val_colnames_i)
    if(minusGroup&min(as.numeric(val_colnames_i))>min(as.numeric(colnames_obs_i))){
      M_i3_minus <- M_i3[,as.numeric(colnames_obs_i) <= min(as.numeric(val_colnames_i)),drop=FALSE]
      M_i3[,paste(min(as.numeric(val_colnames_i)))] <-  rowSums(M_i3_minus) # Add plus group values
      M_i3[which(!is.na(rowSums(M_i3))),which(as.numeric(colnames_obs_i) < min(as.numeric(val_colnames_i)))] <- 0 # Replace values in bins above with zeros
    }

    # Add plus group if specified, and if data are going to be truncated
    # (by default, function simply excludes values above maximum val_colnames_i)
    if(plusGroup&max(as.numeric(val_colnames_i))<max(as.numeric(colnames_obs_i))){
      M_i3_plus <- M_i3[,as.numeric(colnames_obs_i) >= max(as.numeric(val_colnames_i)),drop=FALSE]
      M_i3[,paste(max(as.numeric(val_colnames_i)))] <-  rowSums(M_i3_plus) # Add plus group values
      M_i3[which(!is.na(rowSums(M_i3))),which(as.numeric(colnames_obs_i) > max(as.numeric(val_colnames_i)))] <- 0 # Replace values in bins above with zeros
    }

    # Replace NA in rows (years) with data with zeros
    M_i3 <- t(apply(M_i3, 1, function(x){
      if(all(is.na(x))){
        x}else{
          x[is.na(x)] <- 0
          x
        }
      }))

    M_out_i <- matrix(NA,nrow=length(val_rownames),ncol=length(val_colnames_i),dimnames=list(val_rownames,val_colnames_i))

    # Add observed values to the appropriate cells in the output matrix
    for(obs_rowName_j in rownames(M_i3)[which(rownames(M_i3)%in%rownames(M_out_i))]){
      obs_colname_i <- colnames(M_i3)[which(colnames(M_i3)%in%colnames(M_out_i))]
      #M_out_i[obs_rowName_j,] <- rep(0,ncol(M_out_i)) # First fill with zeros
      # Linear interpolation to match val_colnames_i (won't change anything if existing column names are desired)
      x_ij <- M_i3[obs_rowName_j,]
      y_ij <- rep(NA,length(val_colnames_i))
      if(any(!is.na(x_ij))){
      y_ij <- approx(as.numeric(names(x_ij)),x_ij,xout = as.numeric(val_colnames_i))$y
      }
      M_out_i[obs_rowName_j,] <- y_ij # Now fill with observed data
    }

    M_out_i <- t(apply(M_out_i,MARGIN=1,function(x){x/sum(x)})) # rescale so that each row sums to 1

    # Replace original matrix with new completed matrix
    if(output_type=="prop"){
      comp_i <- M_out_i
    }else{
      # Add sample size column
      n_tag_cols <- matrix(NA,ncol=length(n_colname_new),nrow=nrow(M_out_i),dimnames=list(rownames(M_out_i),n_colname_new))
      comp_i <- as.data.frame(cbind(n_tag_cols,M_out_i))
      RootComp_i <- compName_i
      RootComp_n <- gsub(pattern=n_tag[1],replacement="",x=names(comp_data_n))
      Comp_n_i <- comp_data_n[,paste(RootComp_i,n_tag,sep=""),drop=FALSE]
      rows_Comp_n_i <- match(rownames(comp_i),rownames(Comp_n_i))
      comp_i[,n_colname_new] <- Comp_n_i[rows_Comp_n_i,paste(RootComp_i,n_tag,sep="")]
    }

    if(output_type=="input"){
      comp_i <- comp_i
    }
    if(output_type=="nfish"){
      comp_i <- round(comp_i[,!names(comp_i)%in%n_colname_new]*comp_i$nfish)
    }

    comp_i[comp_i==valsToNA] <- NA # Change specified values to NA (e.g. -99999)

    comp_list[[compName_i]] <- comp_i
  }
  return(comp_list)
}
