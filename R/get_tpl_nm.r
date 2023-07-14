#' Get the names of all objects named in the BAM tpl file, with object classes in \code{tpl_object_classes}.
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param tpl_file tpl file path
#' @param tpl_obj_classes Character vector of names of tpl object classes. I think I listed all of the classes used in BAM, but I might have missed something.
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' get_tpl_nm(tpl_AtlanticMenhaden)
#' get_tpl_nm(tpl_BlackSeaBass)
#' 

get_tpl_nm <- function(tpl_obj=NULL,
                       tpl_file=NULL,
                       tpl_object_classes = c(
                         "init_int",
                         "number", "init_number", "init_bounded_number",
                         "const double",
                         "vector", "init_vector", "init_bounded_vector",
                         "init_ivector",
                         "init_bounded_dev_vector",
                         "matrix", "init_matrix"
                       )
                       
){
  
  if(!is.null(tpl_obj)){
    tpl <- tpl_obj
  }

  # Read in dat, tpl, and cxx files
  if(!is.null(tpl_file)){
    tpl <- readLines(con=tpl_file)
  }

  regex_A <- "(?<=^" #
  regex_Z <- " )(.*?)(?=[=(;])" #" )(.*?)(?=[^A-Za-z1-9_\\s])"

  tpl_names <- list()
  tpl2 <- gsub("(\\/\\/.*)","",tpl) # Remove comments
  tpl3 <- unlist(strsplit(tpl2,split=";"))
  tpl4 <- trimws(tpl3)
  tpl5 <- paste0(tpl4,";") # Add semicolons to end each line

  for(class_i in tpl_object_classes){
    tpl_names[[class_i]] <- unlist(str_extract_all(tpl5,paste0(regex_A,class_i, regex_Z)))
  }
  
  return(tpl_names)
}
