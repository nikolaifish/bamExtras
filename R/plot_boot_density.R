#' Make density plots of bootstrapped parameters and add quantiles as vertical lines.
#'
#' Make density plots of bootstrapped parameters and add quantiles as vertical lines.
#' @param data a numeric matrix or data frame that can be coerced to a numeric matrix where rows are runs and columns represent parameters or variables.
#' @param x_name column name of a variable x to plot
#' @param density_args arguments passed to \code{stats::density} \code{\link[stats]{density}} (excluding x)
#' @param plot_args arguments passed to \code{base::plot} \code{\link[base]{plot}} (excluding x and y)
#' @param quantile_draw logical. Should quantiles of x be added to the plot as vertical lines? (e.g. confidence limits and median)
#' @param quantile_args arguments passed to \code{stats::quantile} \code{\link[stats]{quantile}} for computing quantiles of x to add to plot as vertical lines
#' @param quantile_abline_args arguments to pass to \code{graphics::abline} \code{\link[graphics]{abline}} when plotting quantiles as vertical lines (excluding v)
#' @param col_shade color for shading density region
#' @param add_text_nsim logical. Should the number of simulation runs (nsim) in `data` be added as text to the plots?
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' nsim <- 100
#' data <- data.frame("M"=runif(nsim,0.19, 2.1),"lwb"=rnorm(nsim,mean=3))
#' par(mfrow=c(2,2))
#' plot_boot_density(data,"M")
#' plot_boot_density(data,"M",add_text_nsim = TRUE)
#' plot_boot_density(data,"lwb")
#' plot_boot_density(data,"lwb",add_text_nsim = TRUE,quantile_draw=TRUE,col_shade=rgb(0,0,1,0.5))
#' }

plot_boot_density <- function(data,
                              x_name,
                              density_args = list(from=NULL,to=NULL),
                              plot_args = list(xlab=NULL,ylab="density",main="",type="l"),
                              quantile_draw=FALSE,
                              quantile_args = list(probs=c(0.05,0.50,0.95)),
                              quantile_abline_args = list(lty=c(3,2,3),lwd=c(2,2,2),lend="butt"),
                              col_shade="gray50",
                              add_text_nsim=FALSE){
  par_i <- data[,x_name]

  if(is.null(density_args$from)){density_args$from <- min(par_i)}
  if(is.null(density_args$to)){density_args$to <- max(par_i)}
  if(is.null(plot_args$xlab)){plot_args$xlab=x_name}

  quantile <- do.call(quantile,c(list(x=par_i),quantile_args))
  den <- do.call(stats::density,args=c(list(x=par_i),density_args))

  x <- den$x
  y <- den$y
  do.call(plot,args=c(list(x=den$x,y=den$y),plot_args))
  polygon(x=c(x,rev(x)),y=c(y,rep(0,length(y))),col=col_shade)
  if(quantile_draw){
    do.call(abline,c(list(v=quantile),quantile_abline_args))
  }
  if(add_text_nsim){
    text_legend(legend=paste("n sim =",length(par_i)),"topright")
  }
}
