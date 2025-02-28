#' Plot results from MCBE parameter profile analysis
#'
#' @param sim_summary Output object from summarize_MCBE
#' @param nm_par Name of parameter to plot profile
#' @param nm_obj Name of object where par should be found in sim_summary. By default,
#' plot_profile will look in parms and parm.cons and use the first instance of nm_par it finds.
#' It works with data frames in sim_summary like parms and is programmed to work
#' with parm.cons as a special case, but may not work with other objects within sim_summary.
#' @param like_pattern pattern to use to identify likelihood components from \code{like} object in sim_summary.
#' Internally, the function uses \code{grepl} to identify names in \code{like}. The default "^lk."
#' uses regular expression to identify all likelihood components.
#' @param like_groups names of groups to use when plotting grouped likelihoods. The function adds "^lk."
#' to the group names to find parameters to group in \code{like}.
#' @param like_shift method to use to scale likelihood components for plotting:
#' "min" shifts each likelihood component by the overall minimum,
#' "min_par" shifts each likelihood component by its own minimum,
#' "none" returns raw likelihoods
#' @param like_exclude Names of columns in \code{like} to exclude from plotting profiles
#' @param like_exclude_stacked Names of columns in \code{like} to exclude from plotting stacked barplot
#' @param like_plot_last Name of the likelihood component to
#' plot last so that it is on top in the plots.
#' @param par_args arguments to pass to \code{par}
#' @param matplot_args_all arguments to pass to \code{matplot} for plot of all likelihoods
#' @param matplot_args_group arguments to pass to \code{matplot} for plot of grouped likelihoods
#' @param plot_profiles logical. Plot likelihood profiles as lines, possibly shifted with like_shift.
#' @param plot_stacked logical. Plot likelihood profiles as stacked areas (\code{ggplot})
#' columns not in like_prop_exclude.
#' @param filter_info list of ranges used to filter out sim runs. set to NULL to retain all results.
#' @param par_val_ref reference value to add to profile plots (e.g. parameter value from base model)
#' @param col_sub Specify colors to substitute for colors in the ungrouped profile plot.
#'  named vector where the names match the names match the names in the plot legend.
#'
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' }

plot_profile <- function(sim_summary,
                         nm_par="steep",
                         nm_obj=c("parms","parm.cons"),
                         like_pattern = "^lk.",
                         like_groups = c("L","D","U","lenc","agec","SR"),
                         like_shift = "min_par",
                         like_exclude = c("lk.unwgt.data","lk.Nage.init","lk.fullF","lk.Ftune","lk.SRearly","lk.SRend"),
                         like_exclude_stacked = c("lk.total","lk.unwgt.data"),
                         like_plot_last="lk.total",
                         par_args = list(mfrow=c(2,1),mar=c(3,3,0,0),mgp=c(1,0.3,0),tck=-0.01),
                         matplot_args_all = list(),
                         matplot_args_group = list(),
                         plot_profiles = TRUE,
                         plot_stacked = TRUE,
                         filter_info = list("gradient.max"  = c(0,0.01)),
                         par_val_ref = NULL,
                         col_sub=c("total"=rgb(0,0,0))
){
ss <- sim_summary
parms <- ss$parms
like <- local({
  a <- ss$like
  a <- a[,which(!names(a)%in%like_exclude)]
  if(!is.null(like_plot_last)){
  b <- a[,c(which(names(a)!=like_plot_last),which(names(a)==like_plot_last))]
  }else{b <- a}
  b
})
gradient.max <- ss$like$gradient.max
nsim <- length(parms$modRunName)
simID2 <- simID <- 1:nsim

if(!is.null(filter_info)){
filter_tests <- data.frame(gradient.max = (gradient.max>filter_info$gradient.max[1] & gradient.max<filter_info$gradient.max[2]))
sim_pass <- which(apply(filter_tests,1,all))
sim_fail <- which(!apply(filter_tests,1,all))
n_fail <- length(sim_fail)
if(n_fail>0){message(paste(n_fail,"sim runs were outside of the filter range and were excluded from plots."))}
simID2 <- simID[sim_pass]
                         }

# Get names of all par in obj
#nm_obj_par <- setNames(lapply(nm_obj,function(x){names(ss[[x]])}),nm_obj)

# look for nm_par in each obj
for(nm_obj_i in nm_obj){
  nm_par_obj_i <- names(ss[[nm_obj_i]])
  if(nm_par%in%nm_par_obj_i){
    parx <- if(nm_obj_i=="parm.cons"){
      ss[[nm_obj_i]][[nm_par]][,8]
    }else{
      ss[[nm_obj_i]][[nm_par]]
    }
  parx <- parx[simID2] # filter runs
  }
}

# Compute delta likelihood for each column in like
y <- like[simID2,grepl(like_pattern,names(like))]
ylong <- local({
  a <- y[,!dimnames(y)[[2]]%in%like_exclude_stacked]
  b <- cbind(data.frame("parx"=parx),a)
  dimnames(b)[[2]] <- gsub("^lk.","",dimnames(b)[[2]])
  pivot_longer(data=b,cols=!parx,names_to = "component",values_to = "lk")
})

yd <- apply(y,2,function(x){x-min(x)})
yd2 <- y-min(y)
y_plot <- switch(like_shift,
                 none=y,
                 min_par=yd,
                 min=yd2
                 )
colnames(y_plot) <- gsub("^lk\\.","",colnames(y_plot))
nm_lk <- colnames(y_plot)

# Plot grouped delta likelihood profiles
yg <- ydg <- yd2g <- matrix(NA,nrow=nrow(y),ncol=length(like_groups),dimnames=list(rownames(y),like_groups))

for(i in like_groups){
  # lki <- paste0(like_pattern,i)
  lki <- paste0("^lk.",i)
  yg[,i] <- rowSums(y[,grepl(lki,dimnames(y)[[2]]),drop=FALSE])
}
yg <- yg[,colSums(yg)!=0,drop=FALSE]
ydg <- apply(yg,2,function(x){x-min(x)})
yd2g <- yg-min(yg)
yg_plot <- switch(like_shift,
                 none=yg,
                 min_par=ydg,
                 min=yd2g
)

yglong <- local({
  a <- yg[,!dimnames(yg)[[2]]%in%like_exclude_stacked]
  b <- cbind(data.frame("parx"=parx),a)
  dimnames(b)[[2]] <- gsub("^lk.","",dimnames(b)[[2]])
  pivot_longer(data=b,cols=!parx,names_to = "component",values_to = "lk")
})


# Plot likelihood profiles
ylab <- switch(like_shift,
               none="like",
               min_par="delta like (shift by par min)",
               min="delta like (shift by all min)"
)


col <- setNames(rainbow(ncol(y),start=0,end=0.8),nm_lk)
col[names(col_sub)] <- col_sub

if(plot_profiles){
do.call(par,par_args)

# Plot all components
maa_user <- matplot_args_all
maa_default <- list(type="l",xlab=nm_par,ylab=ylab,lty=1,lwd=2,col=col)
maa <- modifyList(maa_default,maa_user)
do.call(matplot,c(list(x=parx,y=y_plot),maa))
if(!is.null(par_val_ref)){abline(v=par_val_ref,lwd=2,lty=2)}
legend("top",legend=nm_lk,col=maa$col,lty=maa$lty,lwd=maa$lwd,ncol=ceiling(ncol(y_plot)/10),bty="n")

# Plot grouped components
colg <- rainbow(ncol(yg_plot),start=0,end=0.8)

mag_user <- matplot_args_group
mag_default <- list(type="l",xlab=nm_par,ylab=ylab,lty=1,lwd=2,col=colg)
mag <- modifyList(mag_default,mag_user)
do.call(matplot,c(list(x=parx,y=yg_plot),mag))
if(!is.null(par_val_ref)){abline(v=par_val_ref,lwd=2,lty=2)}
legend("top",legend=colnames(yg_plot),col=mag$col,lty=mag$lty,lwd=mag$lwd,ncol=ceiling(ncol(yg_plot)/10),bty="n")
}

if(plot_stacked){
ggy <- ggplot(ylong, aes(x = parx, y = lk, fill=component)) +
  geom_area(position = 'stack') +
  theme_classic() +
  ggtitle("stacked likelihood profiles") +
  xlab(nm_par)

ggyg <- ggplot(yglong, aes(x = parx, y = lk, fill=component)) +
  geom_area(position = 'stack') +
  theme_classic() +
  ggtitle("stacked likelihood profiles") +
  xlab(nm_par)

gridExtra::grid.arrange(ggy,ggyg)
}

invisible(list("parx"=parx,"y"=y_plot,"yg"=yg_plot))
}
