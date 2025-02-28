#' Plot stochastic vectors
#'
#' Plot median values and confidence intervals (shaded region) for stochastic (e.g. bootstrapped) vectors.
#' @param data a numeric matrix or data frame that can be coerced to a numeric matrix where rows are runs and columns represent a  numeric sequence (e.g. age, year). Note that the function would probably work for any numeric matrix, but it's intended for bootstrapped time- or age-series data.
#' @param CIpct percentage to define the range of the confidence intervals to plot
#' @param runsToKeep optional vector of integer values corresponding to the rows of `data` to include in the plot
#' @param ref_x,ref_y x and y vectors for adding a reference line to the plot (e.g. base run values)
#' @param col_lines color for lines depicting central tendency and confidence limits. Expects hex color
#' @param add logical. Should a new plot be drawn (FALSE) or should it be added to an existing plot (TRUE)?
#' @param plot_runs_which Specify individual runs to add to the plot: "none" add none of the individual runs, "all" adds all of the runs. Or specify a numeric vector of the runs you want to plot, corresponding to columns of \code{data}
#' @param plot_runs_args list of arguments to pass to \code{matpoints} when adding runs to the plot
#' @param plot_runs_col_fn name of the color function to use when adding runs to plot. It should be a function like \code{rainbow} or \code{heat.colors} where the first argument is \code{n}
#' @param plot_runs_labels character vector of labels for runs to plot. Length should be equal to the number of individual runs being plotted
#' @param plot_runs_labels_args arguments to pass to \code{text} when labeling runs.
#' @param plot_runs_endpt_args arguments to pass to \code{points} for plotting the last non-NA value of a vector as a point (used in retrospective plots).
#' When pch=NA (default), points will be plotted invisibly.
#' @param polygon_args Optional arguments to pass to \code{polygon} for plotting shaded confidence bands.
#' @param plot_legend logical. Should a legend be added to the plot? If TRUE, legend_args can be specified to format the legend.
#' @param legend_args Optional arguments to pass to \code{legend}. By default,
#' the legend just shows a message about the error bands.
#' @param ... other arguments to pass to \code{matplot}
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' bootData <- t(replicate(n=200,random_walk(x=1:200,mu=0.1)))
#' plot_boot_vec(bootData)
#' plot_boot_vec(bootData,plot_runs_which = sample(1:nrow(bootData),10))
#' plot_boot_vec(bootData,plot_runs_which = sample(1:nrow(bootData),10),plot_shade=FALSE)
#' plot_boot_vec(bootData,plot_runs_which = sample(1:nrow(bootData),10),plot_shade=FALSE,plot_medCI = FALSE)
#' ids <- sample(1:nrow(bootData),10)
#' plot_boot_vec(bootData,plot_runs_which = ids,plot_shade=FALSE,plot_medCI = FALSE,
#' plot_runs_labels=ids,plot_runs_labels_args = list(pos=4,cex=0.75))
#' }

plot_boot_vec <- function(data,
                          CIpct = 95,
                          runsToKeep=NULL,
                          ref_x=NULL,
                          ref_y=NULL,
                          col_lines=rgb(0,0,0,1),
                          add=FALSE,
                          plot_shade=TRUE,
                          plot_medCI=TRUE,
                          plot_runs_which="none",
                          plot_runs_args=list(type="l",lty=1),
                          plot_runs_col_fn=rainbow,
                          plot_runs_labels = NULL,
                          plot_runs_labels_args = list(pos=4,cex=1),
                          plot_runs_endpt_args = list(),
                          polygon_args = list(),
                          plot_legend = FALSE,
                          legend_args = list(),
                          ...){
  # Set defaults which will be overridden by user supplied values
  plot_runs_endpt_args <- modifyList(list(type="p",pch=NA,cex=1.5), plot_runs_endpt_args)
  polygon_args <- modifyList(list(border=NA,lwd=2,col=bamExtras::color_trans(col_lines)), polygon_args)

  CItail <- ((100-CIpct)/2)/100

  if(is.null(runsToKeep)){
    runsToKeep <- 1:nrow(data)
  }
  #par(mar=c(2,2,2,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1.5,tck=-0.02)
  data2 <- t(data[runsToKeep,,drop=FALSE])
  if(is.null(rownames(data2))){rownames(data2) <- 1:nrow(data2)}
  y.median <- apply(X=data2,MARGIN=1,FUN=quantile,probs=0.5,na.rm=TRUE)
  y.lo <- apply(X=data2,MARGIN=1,FUN=quantile,probs=CItail,na.rm=TRUE)
  y.up <- apply(X=data2,MARGIN=1,FUN=quantile,probs=1-CItail,na.rm=TRUE)
  y <- data.frame(y.median,y.lo,y.up)
  x <- as.numeric(rownames(y))
  # draw empty plot
  if(!add){
    matplot(x,y,type="n",...)
    grid(lty=1,col="gray80")
  }

  if(plot_runs_which[1]!="none"){
    rtp_all <- 1:ncol(data2) # all possible runs to plot
    abort_plot_runs <- FALSE
    if(plot_runs_which[1]=="all"){
      rtp <- rtp_all}else{
        if(!is.numeric(plot_runs_which)){
          abort_plot_runs <- TRUE
          warning("runs not plotted. The value of plot_runs_which should be 'all', 'none', or numeric values of the runs you want to plot.")
        }
        rtp <- plot_runs_which[plot_runs_which%in%rtp_all] # only plot runs that exist

      }
    if(!abort_plot_runs){
      data3 <- data2[,rtp,drop=FALSE]
      if(!"col"%in%names(plot_runs_args)){
        plot_runs_args$col <- plot_runs_col_fn(ncol(data3))
      }
      do.call(matpoints,c(list(x=x,y=data3),plot_runs_args))
      # plot last value in vector as a point
        xy_end <- local({
          a <- apply(data3,2,function(yi){
            j <- tail(which(!is.na(yi)),1)
            xij <- x[j]
            yij <- yi[[j]]
            c("x"=xij,"y"=yij)
          })
          t(a)
        })
        if(!"col"%in%names(plot_runs_endpt_args)){
          plot_runs_endpt_args$col <- plot_runs_args$col
        }
        do.call(points,c(list(x=xy_end[,"x"], y=xy_end[,"y"]),plot_runs_endpt_args))

      if(!is.null(plot_runs_labels)){
        if(length(plot_runs_labels)!=length(rtp)){
          prl <- paste(rtp)
        }else{
          prl <- plot_runs_labels
        }
        y_txt <- tail(data3,1)
        x_txt <- rep(tail(x,1),length(rtp))
        do.call(text,c(list(x=x_txt,y=y_txt,labels=prl),plot_runs_labels_args))
        # text(x_txt,y_txt,labels=prl,pos=4)
      }
    }
  }

  ## draw confidence bands as shaded region
  # Identify groups of continuous years for plotting polygons correctly
  if(plot_shade){
    group <- rep(0,nrow(y))
    for(i in 1:nrow(y))
    {
      if(complete.cases(y[i,])) # If test i is true
      {
        if(i==1) # if i=1 then this is group 1
        {
          group[i] <- 1
        }
        if(i>1){ # if i>1
          if(complete.cases(y[i-1,])) # If test i-1 is true
          {
            group[i] <- group[i-1]
          }else{
            group[i] <- max(group)+1
          }
        }
      }
    }
    y$group <- group

    ycc <- y[complete.cases(y),]
    for(group.i in 1:max(ycc$group)){
      ycc.group.i <-  ycc[ycc$group==group.i,]
      do.call(polygon,c(list(x=c(rownames(ycc.group.i),rev(rownames(ycc.group.i))),y=c(ycc.group.i[,2],rev(ycc.group.i[,3]))),
                        polygon_args))
    }
  }

  if(plot_medCI){
    # draw median and confidence band lines
    matpoints(rownames(y),y[,c(1,2,3)],type="l",lwd=c(2,1,1),lty=c(2,1,1),col=col_lines)
    # draw reference lines (e.g. base run)
    if(!is.null(ref_x)&!is.null(ref_y)){
      points(ref_x,ref_y,type="o",lty=1,lwd=3,pch=16,col=col_lines)
    }
  }

  if(plot_medCI|plot_shade){
    legend_args <- modifyList(list(x="topright",legend=paste("error bands represent ",CIpct,"% CI",sep=""),bty="n",cex=0.75), legend_args)
    plot_legend <- TRUE
  }
  if(plot_legend){
    if(length(legend_args)==0){
      message("plot_boot_vec: legend_args has length zero. Legend will not be plotted")
    }else{
    do.call(legend,legend_args)
    }
  }
}
