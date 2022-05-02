#' Add elements to rdat
#'
#' @param rdat BAM output rdat (list) object read in with dget()
#' @param fleet_sel_abb_key Data frame for specifying which selectivities to use with time series of fishing mortality (F) for individual fleets
#' when the abbreviations in the rdat do not exactly match with a selectivity matrix. Column F lists the abbreviations F time series. Column sel
#' lists the abbreviations associated with the selectivity you want to use with the F time series. This is used to compute fleet-specific F-age
#' matrices (year,age).
#' @param F_names_exclude Character vector of names of columns in t.series that you want to exclude. The function identifies any columns
#'  that begin with "F." The default value of \code{F_names_exclude} is a character vector of column names found in some BAM rdat files that
#'  that begin with "F." but are not fleet-specific F time series.
#' to a single fleet.
#' @details
#'
#' @keywords bam stock assessment fisheries MSEtool
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Add to any of the current BAM models
#' rdat_AtMe <- add2rdat(rdat_AtlanticMenhaden)
#' rdat_BlSB <- add2rdat(rdat_BlackSeaBass)
#' rdat_BlTi <- add2rdat(rdat_BluelineTilefish)
#' rdat_Cobi <- add2rdat(rdat_Cobia)
#' rdat_GagG <- add2rdat(rdat_GagGrouper,fleet_sel_abb_key=data.frame("F"=c("cHL.D","rHB.D","rGN.D"),"sel"=c("D","D","D")))
#' rdat_GrTr <- add2rdat(rdat_GrayTriggerfish)
#' rdat_GrAm <- add2rdat(rdat_GreaterAmberjack)
#' rdat_ReGr <- add2rdat(rdat_RedGrouper)
#' rdat_RePo <- add2rdat(rdat_RedPorgy)
#' rdat_ReSn <- add2rdat(rdat_RedSnapper)
#' rdat_SnGr <- add2rdat(rdat_SnowyGrouper)
#' rdat_Tile <- add2rdat(rdat_Tilefish)
#' rdat_VeSn <- add2rdat(rdat_VermilionSnapper,fleet_sel_abb_key=data.frame("F"=c("rGN.D"),"sel"=c("rHB.D")))
#' }

add2rdat <- function(rdat,
                     fleet_sel_abb_key=data.frame("F"=c("cHL.D","rHB.D","rGN.D"),"sel"=c("D","D","D")),
                     F_names_exclude=c("F.Fmsy","F.full","F.sum","F.F30.ratio","F.F40.ratio")

)
  {

  rdat <- bamExtras::standardize_rdat(rdat)

  info <- rdat$info
  parms <- rdat$parms
  parm.cons <- rdat$parm.cons
  parm.tvec <- rdat$parm.tvec
  a.series <- rdat$a.series
  t.series <- rdat$t.series
  comp.mats <- rdat$comp.mats
  sel.age <- rdat$sel.age
  B.age <- rdat$B.age

  styr <- parms$styr
  endyr <- parms$endyr
  yr <- styr:endyr

  Name <- gsub(" ","",stringr::str_to_title(info$species))

  # MSEtool expects age-based data to begin with age 0
  if(min(a.series$age)>0){
    message(paste(Name,": Minimum age > 0. Age-based data (a.series) linearly extrapolated to age-0"))
    a.series <- bamExtras::data_polate(a.series,xout=0:max(a.series$age))
    a.series <- data_lim(a.series,xlim=c(0,Inf))
    a.series <- data_lim(a.series,xname=c("prop.female","prop.male","mat.female","mat.male"),xlim=c(0,1))
    a.series <- as.data.frame(a.series)
    rownames(a.series) <- a.series$age
  }
  age <- a.series$age

  t.series <- t.series[paste(styr:endyr),]

  for(i in names(sel.age)){
    sel.i <- sel.age[[i]]
    if(grepl("^sel.v",i)){
    if(min(as.numeric(names(sel.i)))>0){
      message(paste0(Name,": Minimum age of ",i," > 0. ", i, " linearly extrapolated to age-0"))
        sel.i <- cbind("age"=as.numeric(names(sel.i)),"sel"=sel.i)
        sel.i <- bamExtras::data_polate(sel.i,xout=0:max(sel.i[,"age"]))
        sel.i <- data_lim(sel.i,xname="sel",xlim=c(0,1))
        # sel.i <- as.data.frame(sel.i)
        sel.age[[i]] <- setNames(sel.i[,"sel"],sel.i[,"age"])
      }
    }
    if(grepl("^sel.m",i)){
      if(min(as.numeric(colnames(sel.i)))>0){
        message(paste0(Name,": Minimum age of ",i," > 0. ", i, " linearly extrapolated to age-0"))
        sel.i <- cbind("age"=as.numeric(colnames(sel.i)),t(sel.i))
        sel.i <- bamExtras::data_polate(sel.i,xout=0:max(sel.i[,"age"]))
        sel.i <- data_lim(sel.i,xname=colnames(sel.i)[colnames(sel.i)!="age"],xlim=c(0,1))
        sel.i <- t(sel.i)
        colnames(sel.i) <- sel.i["age",]
        sel.i <- sel.i[rownames(sel.i)!="age",]
        sel.age[[i]] <- sel.i
      }
    }
  }

  # Reformat sel.age to only include separate fleets and to format all elements as matrices (year, age)
  sel.age.fleet.m <- local({
    out <- sel.age
    a <- names(out)
    b <- a[!grepl("^sel.v.wgted",a)]
    out <- out[b]

    b.v <- b[grepl("^sel.v",b)]

    out_new_m <- list()
    for(i in b.v){
      out_i <- matrix(out[[i]],nrow=length(yr),ncol=length(age),dimnames=list(yr,age),byrow=TRUE)
      outName_i <- gsub("^sel.v","sel.m",i)
      out_new_m[[outName_i]] <- out_i
      out <- out[names(out)!=i]
    }
    out <- c(out_new_m,out)
  })

  ## Identify fleet-specific selectivity
  sel_names <- names(sel.age.fleet.m)
  sel_abb <- gsub("^sel.(m.)","",sel_names)


  ## Identify fleet-specific F
  F_names <- local({
    a <- names(t.series)[grepl("^(F\\.)",names(t.series))]
    a[!a%in%F_names_exclude]
  })
  F_abb <- gsub("^(F\\.)","",F_names)

  # Compute F-matrices (year,age) for each fleet
  # Match fleet-specific F to selectivities by abbreviation
  F.age.fleet <- list()
  for(i in F_abb){
    F_name_i <-   paste0("F.",i)
    sel_name_i <- paste0("sel.m.",i)
    F_i <- t.series[,F_name_i]
    # Check if sel_name_i is found in sel.age.fleet.m and then match it up
    if(!sel_name_i%in%names(sel.age.fleet.m)){
      ix_i <- match(i,fleet_sel_abb_key$F)
      if(is.na(ix_i)){
        message(paste(i, "does not match any abbreviation in fleet_sel_abb_key$F"))
      }
      sel_name_i <- paste0("sel.m.",fleet_sel_abb_key$sel[ix_i])
    }
    sel_i <- sel.age.fleet.m[[sel_name_i]]
    F.age.fleet[[i]] <- sel_i*F_i
  }

  rdat$F.age.fleet <- F.age.fleet

  return(rdat)
}
