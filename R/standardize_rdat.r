#' Change non-standard names of rdat objects to standard naming conventions. Compute standard objects.
#'
#' @param rdat rdat (list) object read in with dget(). Specifically, the BAM output file produced by the cxx file.
#' @param separator_key named vector indicating what separator(s) should be used in the rdat. By default, any "_" will be replaced with "." in names of rdat elements
#' @param separator_key_x names of elements in the rdat list to apply the fleet_key to.
#' @param fleet_replace Should fleet_key be applied to fleet_key_x to replace fleet_abbreviations in the rdat?
#' @param fleet_key list where values are patterns (character vector) to find in rdat list elements specified with fleet_key_x and the names are the replacements. Patterns are searched with regex to find these values at the beginning or end of character strings, followed or preceded by ".", or in the middle of a string preceded and followed by ".".
#' @param fleet_key_x names of elements in the rdat list to apply the fleet_key to.
#' @param parms_key named vector where values are patterns (character strings) to find in names(rdat$parms) and the names in parms_key are replacements
#' @param a_series_key named vector where values are patterns (character strings) to find in names(a.series$parms) and the names in a_series_key are replacements
#' @param t_series_key named vector where values are patterns (character strings) to find in names(t.series$parms) and the names in t_series_key are replacements
#' @param na_values Other values to change to NA
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' names(rdat_RedPorgy$a.series)
#' rdat_RedPorgy_std <- standardize_rdat(rdat_RedPorgy)
#' names(rdat_RedPorgy_std$a.series)
#' }

standardize_rdat <- function(rdat,
                             separator_key = c("."="_"),
                             separator_key_x = names(rdat),
                             fleet_replace = TRUE,
                             fleet_key=list("sCT"=c("CVT"),              # Chevron trap (could possibly include video)
                                            "sTV"=c("CVID","Mcvt"),      # Combined chevron trap/video data (sCT in Red Porgy; Mcvt used in Black Seabass bass for the combined index but also for the comps which are really only from the trap)
                                            "sVD"=c("VID"),              # Video data (from chevron trap survey)
                                            "sBT"=c("Mbft"),             # MARMAP blackfish trap (see Black Seabass)
                                            "sBL"=c("sM"),               # MARMAP bottom longline survey (see Tilefish)
                                            "sVL"=c("vll"),              # MARMAP vertical longline survey (see Snowy Grouper)
                                            "sFT"=c("FST"),              # MARMAP Florida Snapper Trap (see Vermilion Snapper)
                                            "sAN"=c("nad"),              # Northern Adult composite survey index (Menhaden)
                                            "sAM"=c("mad"),              # Middle Adult composite survey index (Menhaden)
                                            "sAS"=c("sad"),              # Southern Adult composite survey index (Menhaden)
                                            "sJA"=c("jai"),              # Juvenile Abundance composite survey index (Menhaden)
                                            "sME"=c("mareco"),           # MARMAP and ECOMON survey index (Menhaden)
                                            "cDV"=c("cD"),               # Commercial diving (see Gag)
                                            "cHL"=c("cH","cHl"),         # Commercial handlines (a.k.a. commercial lines)
                                            "cLL"=c("cL"),               # Commercial longlines (see Blueline Tilefish)
                                            "cOT"=c("cO"),               # Commercial other (see Red Grouper)
                                            "cPT"=c("cp","cP"),          # Commercial pots (see Black Seabass)
                                            "cTW"=c("cT","cTw", "cHTR"), # Commercial trawl (see Black Seabass, Red Porgy, Vermilion Snapper)
                                            "cGN"=c("comm"),             # Commercial all. general commercial (see Black Seabass "D.comm.ob")
                                            "cRN"=c("cRn"),              # Commercial Reduction fishery, North  (Menhaden)
                                            "cRS"=c("cRs"),              # Commercial Reduction fishery, South  (Menhaden)
                                            "cBN"=c("cBn"),              # Commercial Bait fishery, North  (Menhaden)
                                            "cBS"=c("cBs"),              # Commercial Reduction fishery, South  (Menhaden)
                                            "rHB"=c("HB","hb","rHb"),    # Recreational headboat
                                            "rHD"=c("hbd","HBD"),        # Recreational headboat discards (atypical abbreviation found in Black Sea Bass selectivity parameters)
                                            "rGN"=c("GR","mrip","rGe","rA")  # Recreational all (a.k.a. general recreational (i.e. not headboat)
                             ),
                           fleet_key_x=c("parms","parm.cons","t.series","comp.mats","sel.age","sel.parms"),
                           parms_key=c("R0"="BH.R0",
                                       "M.msst"="M.constant",
                                       "M.MSST"="M.constant"
                                       ),
                           a_series_key=c("mat.male.endyr"="mat.male",
                                          "mat.fem.endyr"="mat.female",
                                          "prop.female.endyr"="prop.female"
                                          ),
                           t_series_key=c("total.L.wgt.klb"="total.L.klb"),
                           na_values="-99999"
                           ){
  # rdat
  rdat_char <- deparse(rdat) # Convert rdat list into character vector preserving list structure
  for(i in na_values){
    replacement <- "NA"
    pattern <- i
    rdat_char <- gsub(pattern=pattern,replacement = replacement,x=rdat_char)
  }
  rdat <- eval(parse(text=rdat_char))

  # Replace "_" with "." in all rdat element names

  for(i in names(rdat)){
    oi <- rdat[[i]]
    nam <- names(oi)
    if(!is.null(nam)){
      nam <- gsub("_",".",nam)
      names(oi) <- nam
      rdat[[i]] <- oi
    }
  }

  if(fleet_replace){
    # apply general fleet key
    for(i in names(fleet_key)){
      pattern_beg_i <- paste0("^(",paste(fleet_key[[i]],collapse="|"),")\\.")
      pattern_mid_i <- paste0("\\.(",paste(fleet_key[[i]],collapse="|"),")\\.")
      pattern_end_i <- paste0("\\.(",paste(fleet_key[[i]],collapse="|"),")([0-9]*)$")

      for(j in fleet_key_x){
        if(j%in%names(rdat)){
          oj <- rdat[[j]]
          nam <- names(oj)
          nam <- gsub(pattern_beg_i,paste0(i,"."),nam)
          nam <- gsub(pattern_mid_i,paste0(".",i,"."),nam)
          nam <- gsub(pattern_end_i,paste0(".",i,"\\2"),nam)
          names(oj) <- nam
          rdat[[j]] <- oj
        }
      }
    }
  }

  # info
  info <- rdat$info
  info$species <- gsub("SA ","",info$species)
  rdat$info <- info

  # parms
  parms <- rdat$parms
  parms_names <- names(parms)
  names(parms) <- find_replace(names(parms_key),parms_key,parms_names)
  rdat$parms <- parms

  # a_series
  a_series <- rdat$a.series
  a_series_names <- names(a_series)
  names(a_series) <- find_replace(names(a_series_key),a_series_key,a_series_names)

  # If a_series includes prop.male instead of prop.female, compute prop.female
  if(!"prop.female"%in%names(a_series)&"prop.male"%in%names(a_series)){
    a_series$prop.female <- 1-a_series$prop.male
  }
  rdat$a.series <- a_series

  # t_series
  t_series <- rdat$t.series
  t_series_names <- names(t_series)
  names(t_series) <- find_replace(names(t_series_key),t_series_key,t_series_names)
  rdat$t.series <- t_series

  return(rdat)
}
