#' Complete a length or age composition data frame so that all have desired dimensions. This is a modified version of comp_complete. The two functions should probably be combined with options to choose between sub-functions.
#'
#' @param D data frame of mean observations by year, where columns are different groups (e.g. fleets)
#' @param DErr error values for each value in D. Must be same dimensions as D
#' @param ErrStyle Style of displaying error around the central tendency time series D. Can be"bars" or "bands"
#' @param ErrUnits Units of DErr. Used only for displaying text on plot (e.g. "SE")
#' @param NErrUnits Multiplied by DErr to determine how far above and below D to to plot the errors
#' @param ylim range of y on plot
#' @param plotLabel character string for labeling plot
#' @param ... pass additional arguments to matplot
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' tseries_plot()
#' rdat <- rdat_BlackSeaBass
#' tser <- rdat$t.series
#' U <- tser[,grepl("^U.*ob$",names(tser))]
#' Ucv <- tser[,grepl("^cv.U.",names(tser))]
#' tseries_plot(U,Ucv,xlab="year",ylab="index of abundance",plotLabel=rdat$info$species)
#' }
#'

tseries_plot <- function(D,DErr,ErrStyle="bands",ErrUnits="",NErrUnits=2,
                         ylim=NULL,plotLabel=NULL,...){

  DLo <- D-NErrUnits*DErr
  DUp <- D+NErrUnits*DErr

  nSer <- ncol(D)

  if(is.null(ylim)){
    ylim <- c(range(c(0,D,DUp),na.rm=TRUE))
  }


  # plot
  par(mfrow=c(1,1),mar=c(3,3,1,3),mgp=c(1.5,.3,0),oma=c(0,0,0,0),tck=0.01,cex.lab=1.5,cex.axis=1.5)
  legend_text <- names(D)
  col <- rainbow(nSer)
  col_bands <- rainbow(nSer,alpha=0.3)
  type <- "o"
  lty <- rep(1,ncol(D))
  lwd <- 2
  pch <- 16

  matplot(x=as.numeric(rownames(D)),y=D,type=type,lty=lty,lwd=lwd,col=col,pch=pch,
          ylim=ylim,...)

  for(i in 1:nSer){
    x <- as.numeric(rownames(D))
    ylo <- DLo[,i]
    yup <- DUp[,i]

    cc <- complete.cases(x,ylo,yup)

    x <- x[cc]
    ylo <- ylo[cc]
    yup <- yup[cc]

    if(ErrStyle=="bars"){
      arrows(x0=x,
             y0=ylo,
             y1=yup,
             length=0,col=col[i]
      )
    }
    if(ErrStyle=="bands"){
      polygon(x=c(x,rev(x)),
              y=c(ylo,rev(yup)),
              col=col_bands[i],
              border=NA
      )
    }
  }

  legend("topright",legend=legend_text,
         bty="n", col=col, lty=lty, lwd=lwd)

  if(!is.null(plotLabel)){
    text_legend(plotLabel)
  }

  legend("right",bty="n",legend=paste("Error",ErrStyle,"\nplotted to",NErrUnits,ErrUnits))

  invisible(list("D"=D,"DErr"=DErr,"DLo"=DLo,"DUp"=DUp))
}
