#' Plot results from Beaufort Assessment Model
#'
#' A generic plotting function that calls plotting functions from FishGraph and bamExtras.
#' @param spp bam output rdat object
#' @param graphics.type from FishGraph help: 'a vector of graphics file types to which graphics are saved. When NULL, no plots are saved.'
#' @param draft from FishGraph help: 'modifies plots for use in a report. When FALSE main titles are omitted.'
#' @param use.color from FishGraph help: 'plots are made in grayscale when FALSE'
#' @param years_plot vector of years to plot in time series. If NULL, all years are plotted
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Run MCBE, writing files to dir_bam_sim
#' plot_bam(rdat_BlackSeaBass)
#' }

plot_bam <- function(spp=NULL,
                     graphics.type="pdf", # NULL (no quotes) for no plots saved; other options: "pdf", "wmf", "eps"
                     draft=FALSE, # draft type
                     use.color=TRUE,   # color type
                     years_plot = NULL,
                     default_args=list(x=quote(spp), DataName="spp",draft=draft, use.color=use.color, graphics.type=graphics.type),
                     Bound.vec.plots.args = list(),
                     BSR.time.plots.args = list(),
                     CLD.total.plots.args = list(),
                     Cohort.plots.args = list(),
                     Comp.plots.args = list(),
                     Comp.yearly.plots.args = list(),
                     Eq.plots.args = list(),
                     F.time.plots.args = list(),
                     Growth.plots.args = list(),
                     Index.plots.args = list(),
                     Landings.plots.args = list(),
                     NFZ.age.plots.args = list(),
                     Parm.plots.args = list(),
                     PerRec.plots.args = list(),
                     Phase.plots.args = list(),
                     Selectivity.plots.args = list(),
                     StockRec.plots.args = list()
                     )
                     {

library(FishGraph)

  # set_ts_years()
  # Set years of time series in data frame, when rownames are years
  set_ts_years <- function(x,years){
    # Initialize data frame with desired dimensions
    x <- as.matrix(x)
    x.out <- matrix(NA,nrow=length(years),ncol=ncol(x),
                    byrow=TRUE,dimnames=list(years,colnames(x)))
    x.out[rownames(x)[rownames(x)%in%paste(years)],] <- x[rownames(x)[rownames(x)%in%paste(years)],]
    x.out
  }

styr <- spp$parms$styr
endyr <- spp$parms$endyr

if(is.null(years_plot)){
  years_plot <- as.numeric(styr:endyr)
}


  ### rdat file
  # Modify spp
  # This sorting is helpful because apparently Fishgraph gives you errors if the observed comps for a fleet aren't immediately followed by the predicted
  # comps for that same fleet.
  spp$comp.mats <- spp$comp.mats[sort(names(spp$comp.mats))]

  # Trim time series
  spp$N.age <- spp$N.age[rownames(spp$N.age)>=spp$parms$styr,]
  spp$N.age.mdyr <- spp$N.age.mdyr[rownames(spp$N.age.mdyr)>=spp$parms$styr,]
  spp$N.age.spawn <- spp$N.age.spawn[rownames(spp$N.age.spawn)>=spp$parms$styr,]

  # Modify time series to plot a certain set of years
  spp$t.series <- as.data.frame.matrix(set_ts_years(spp$t.series,years_plot))

  spp$CLD.est.mats <- lapply(spp$CLD.est.mats,FUN=set_ts_years,years=years_plot)

  spp$N.age <- set_ts_years(spp$N.age,years_plot)
  spp$N.age.mdyr <- set_ts_years(spp$N.age.mdyr,years_plot)
  spp$N.age.spawn <- set_ts_years(spp$N.age.spawn,years_plot)
  spp$B.age <- set_ts_years(spp$B.age,years_plot)

  # Truncate F-range in eq.series to make nicer plots
  spp$eq.series <- spp$eq.series[which(spp$eq.series$F.eq<=1),]

  # Calculate values to add to (CLD) plots
  spp$parms$Dmsy.num <- spp$parms$Dmsy.knum*1000 # Add value of Dmsy in numbers
  spp$parms$msy.num <- spp$parms$msy.knum*1000 # Add value of msy in numbers

  #Round Effective sample sizes for more enjoyable comp
  spp$t.series[,grepl(pattern="comp",names(spp$t.series))&grepl(pattern="neff",names(spp$t.series))] <-
    round(spp$t.series[,grepl(pattern="comp",names(spp$t.series))&grepl(pattern="neff",names(spp$t.series))],2)

  ###### PLOTAPALOOZA !!! ######
  # Generic plot function for calling various FishGraph plot functions, with a couple of tweaks
  plot_fn <- function(fn){
    cs <-   'if(!is.null(fn.args)){
    cat("\nRunning fn\n")
    default_args <- default_args[names(default_args)%in%names(formals(fn))]
    tryCatch(expr=do.call(FishGraph::fn,c(default_args,fn.args)),
             error = function(e) {message(paste("Error in fn:",e))}
    )

  }'
    eval(parse(text=gsub("fn",fn,cs)))
  }

  ### Bound.vec.plots ###
  plot_fn("Bound.vec.plots")

  ### BSR ###
  plot_fn("BSR.time.plots")

  ### CLD ###
  plot_fn("CLD.total.plots")

  ### comp ###
  plot_fn("Comp.plots")

  ### compyr ###
  plot_fn("Comp.yearly.plots")
  plot_fn("Cohort.plots")
  # if(!is.null(Cohort.plots.args)){
  #   cat("\nRunning Cohort.plots\n")
  #   # only include args that are in the function formals
  #   default_args_Cohort <- default_args[names(default_args)%in%names(formals(Cohort.plots))]
  #   tryCatch(expr=do.call(FishGraph::Cohort.plots,c(default_args_Cohort,Cohort.plots.args)),
  #            error = function(e) {message(paste("Error in Cohort.plots:",e))}
  #   )
  # }

  ### EQ ###
  plot_fn("Eq.plots")

  ### F ###
  plot_fn("F.time.plots")

  ### growth ###
  plot_fn("Growth.plots")

  ### index ###
  plot_fn("Index.plots")

  ### L ###
  plot_fn("Landings.plots")

  ### NFZ.age ###
  plot_fn("NFZ.age.plots")

  ### parms ###
  plot_fn("Parm.plots")

  ### phase ###
  plot_fn("Phase.plots")

  ### PR ###
  plot_fn("PerRec.plots")

  ### sel ###
  plot_fn("Selectivity.plots")

  ### SR ###
  plot_fn("StockRec.plots")

}
