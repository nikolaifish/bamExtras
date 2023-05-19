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
#' @param like_prop_exclude Names of columns in \code{like} to exclude when computing
#' @param plot_profiles logical. Plot likelihood profiles as lines, possibles shifted with like_shift?
#' @param plot_stacked logical. Plot likelihood profiles as stacked areas (\code{ggplot})?
#' columns not in like_prop_exclude.
#' @param par_args arguments to pass to \code{par}
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
                         like_prop_exclude = c("lk.total","lk.unwgt.data"),
                         par_args = list(mfrow=c(2,1),mar=c(3,3,0,0),mgp=c(1,0.3,0),tck=-0.01),
                         plot_profiles = TRUE,
                         plot_stacked = TRUE
){
ss <- sim_summary
parms <- ss$parms
like <- ss$like

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
  }
}

# Compute delta likelihood for each column in like
y <- like[,grepl(like_pattern,names(like))]
ylong <- local({
  a <- y[,!dimnames(y)[[2]]%in%like_prop_exclude]
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
  a <- yg[,!dimnames(yg)[[2]]%in%like_prop_exclude]
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

col <- rainbow(ncol(y),start=0,end=0.8)

if(plot_profiles){
do.call(par,par_args)

# Plot all components
matplot(parx,y_plot,type="l",xlab=nm_par,ylab=ylab,lty=1,lwd=2,col=col)
legend("top",legend=colnames(y_plot),col=col,lty=1,lwd=2,ncol=ceiling(ncol(y_plot)/10),bty="n")

# Plot grouped components
colg <- rainbow(ncol(yg_plot),start=0,end=0.8)
matplot(parx,yg_plot,type="l",xlab=nm_par,ylab=ylab,lty=1,lwd=2,col=colg)
legend("top",legend=colnames(yg_plot),col=colg,lty=1,lwd=2,ncol=ceiling(ncol(yg_plot)/10),bty="n")
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

invisible(list("y"=y_plot,"yg"=yg_plot))
}
