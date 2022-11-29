#' Complete and combine age or length comps
#'
#' @param comp_list list of composition matrices, where column names are bin values (e.g. length or age), and rows are unique observations (e.g. years)
#' @param FUN function (unquoted name) to apply to analogous cells of comp matrices. Default is to sum values.
#' @param scale_rows logical. Should rows of resulting matrix be scaled to sum to 1? Defaults to TRUE.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Combine age comp matrices for Red Porgy
#' cm <- rdat_RedPorgy$comp.mats
#' cma <- cm[grepl("^acomp.*ob$",names(cm))]
#' comp_combine(cma)
#' }
#'

comp_combine <- function(comp_list, FUN=function(x){sum(x,na.rm=TRUE)},scale_rows=TRUE){
  # Identify the complete set of row and column names for all comp matrices
  all.rownames <- sort(unique(unlist(lapply(comp_list,rownames))))
  all.colnames <- as.character(sort(as.numeric(unique(unlist(lapply(comp_list,colnames))))))

  for (i in names(comp_list)){
    M.i <- comp_list[[i]]

    # Identify the specific row and column names missing comp matrix i
    row.missing <- all.rownames[!all.rownames%in%rownames(M.i)]
    col.missing <- all.colnames[!all.colnames%in%colnames(M.i)]

    ## Add missing rows
    # Create empty matrix of missing rows
    M.row.missing <- matrix(NA,nrow=length(row.missing),ncol=ncol(M.i),dimnames=list(row.missing,colnames(M.i)))
    # Actually add missing rows and columns to comp matrix i
    M.i2 <- rbind(M.i,M.row.missing)
    M.i2 <- local({M.i2 %>% rownames %>% order -> o; M.i2[o,]}) # Sort rows of M.i2

    ## Add missing columns
    # Create empty matrix of missing columns
    M.col.missing <- matrix(NA,nrow=nrow(M.i2),ncol=length(col.missing),dimnames=list(row.names(M.i2),col.missing))
    # Actually add missing rows and columns to comp matrix i
    M.i3 <- cbind(M.i2,M.col.missing)
    M.i3 <- local({M.i3 %>% colnames %>% as.numeric %>% order -> o; M.i3[,o]}) # Sort columns of M.i3

    # Replace original matrix with new completed matrix
    comp_list[[i]] <- M.i3
  }
  #  return(comp_list)
  # Convert comp_list to array
  A.cmp <- array(data=unlist(comp_list), dim = c(length(all.rownames), length(all.colnames), length(comp_list)))

  # Apply function to analogous cells across all matrices
  M.cmp.combined <- apply(X=A.cmp,MARGIN=c(1,2),FUN=FUN)

  # Scale rows to sum to 1
  if(scale_rows){
    M.cmp.combined <- M.cmp.combined/rowSums(M.cmp.combined)
  }

  # Rename rows and columns
  rownames(M.cmp.combined) <- all.rownames
  colnames(M.cmp.combined) <- all.colnames

  # For rows with no data, fill with NA
  M.cmp.combined[rowSums(M.cmp.combined)==0] <- NA

  return(M.cmp.combined)
}
