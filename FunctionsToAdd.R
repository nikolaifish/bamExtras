# Note
# Run this in RStudio Cloud to build online documentation page/site
# pkgdown seems pretty cool
# pkgdown::build_site(override = list(destination = "~/docs"))
# apparently you can/should use the NMFS pkgdown template: https://github.com/nmfs-fish-tools/pkgdownTemplate


## FUNCTIONS FROM 00_Setup.Rmd from SEDAR 66 MCB

# color function like rainbow() but without yellow (or indigo)
ROGBP <- colorRampPalette(c("red","orange","green","blue","purple"))
ROGBP.alpha <- function(alpha){
  x <- colorRampPalette(c(
    rgb(1,0,0,alpha), # red
    rgb(1,0.5,0.25,alpha), # orange
    rgb(0,1,0,alpha), # green 
    rgb(0,0,1,alpha), # blue
    rgb(0.5,0,0.5,alpha)), # purple
    alpha=TRUE)
  return(x)
}

#qtrunc and rtrunc are from:
#Nadarajah&Kotz. 2006. R Programs for Computing Truncated Distributions. J of Stat Software
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
}


# Formula for a circle
# x: x values
# r: radius
# center: coordinates of the center of the circle
# diff.digits: number of digits to round differences to to avoid weird and obnoxious differences in precision
#              I was sometimes getting, resulting in NaN when trying to take the sqrt of a tiny negative number
circle <- function(x,r,center=c(0,0),diff.digits=10){
  sqrt(round(r^2-(x-center[1])^2,diff.digits))+center[2]}

# Random generation of a semi-circular distribution
# n: number of random values to draw
# r: radius of circle
# x.center: x.coordinate of center of a circle (y set to zero)
# method: which method should be used: "filter.matrix"=randomly draws from rectangular area and then filters out values that fall
#         within semi-circular area. "discrete.sample"=uses sample() function based on discrete x values with prob values based on
#         circle formula
# n.discrete.x: number of discrete values of x to try if using discrete.sample method.
rcircle <- function(n,r,x.center=0,method="discrete.sample",n.discrete.x=10000){
  if(method=="discrete.sample"){
    xlim <- x.center+c(-1,1)*r
    x.discrete <- seq(xlim[1],xlim[2],length=n.discrete.x) # Vector of possible x values
    x.prob <- circle(x=x.discrete,r=r,center=c(x.center,0))
    x.out <- sample(x=x.discrete,size=n,prob = x.prob,replace=TRUE)
    x.out      
  }
  
  if(method=="filter.matrix"){
    n.try <- n*3
    xlim <- x.center+r*c(-1,1)
    x.rand <- runif(n.try,min=xlim[1],max=xlim[2])
    y.rand <- runif(n.try,min=0,max=r)
    y.pred <- circle(x=x.rand,r=r,center=c(x.center,0))
    
    x.try <- x.rand[which(y.rand<=y.pred)]
    x.out <- x.try[1:n]
    x.out  
  }
  return(x.out)  
}

source("ErrorBand.r")

# Label plot with a legend that is just text
text.legend <- function(legend,x="topleft",bty="n",adj=c(0.1,0),cex=1.3,text.col=rgb(0,0,0.3,0.7),...){
  legend(x=x,legend=legend,
         adj=adj,bty=bty,cex=cex,text.col=text.col,...)
}

# plotBootDensity()    
plotBootDensity <- function(data,par.name,xlab,ylab="density",main="",probs=c(0.05,0.50,0.95),from,to,
                            dataScale=1,drawCL=FALSE,col.shade=rgb(0,0,0,0.25),add.text.Nboot=FALSE,...){
  if(missing(xlab)){xlab=par.name}
  par.i <- data[,par.name]/dataScale
  CL <- quantile(par.i,probs=probs)
  #CL99 <- quantile(par.i,probs=c(0.005,0.995)) # 99% confidence limits
  if(missing(from)){from <- min(par.i)}
  if(missing(to)){to <- max(par.i)}
  den <- density(par.i,from=from,to=to)
  x <- den$x
  y <- den$y
  plot(x,y,type="l",
       xlab=xlab,
       ylab=ylab,
       main=main,
       ...)
  polygon(x=c(x,rev(x)),y=c(y,rep(0,length(y))),col=col.shade)
  if(drawCL){
    abline(v=CL,lty=c(3,2,3),lwd=2,lend="butt")
  }
  if(add.text.Nboot){
    text.legend(legend=paste("n boot =",length(par.i)),"topright")
  }
}

# plotBootSeries
# dataBoot: data set where rows are bootstrap runs, columns might be age or years
# runsToKeep  
plotBootSeries <- function(dataBoot,CIpct = 95,runsToKeep,x.ref,y.ref,col.shade=rgb(0,0,0,0.25),...){
  if(missing(runsToKeep)){
    runsToKeep <- 1:nrow(dataBoot)
  }
  CItail <- ((100-CIpct)/2)/100
  #par(mar=c(2,2,2,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1.5,tck=-0.02)
  x <- t(dataBoot[runsToKeep,])
  x.median <- apply(X=x,MARGIN=1,FUN=quantile,probs=0.5,na.rm=TRUE)
  x.lo <- apply(X=x,MARGIN=1,FUN=quantile,probs=CItail,na.rm=TRUE)
  x.up <- apply(X=x,MARGIN=1,FUN=quantile,probs=1-CItail,na.rm=TRUE)
  x2 <- data.frame(x.median,x.lo,x.up)
  
  # draw empty plot
  matplot(rownames(x2),x2,type="n",...)
  grid(lty=1,col="gray80")
  
  ## draw confidence bands as shaded region
  # Identify groups of continuous years for plotting polygons correctly
  group <- rep(0,nrow(x2))
  for(i in 1:nrow(x2))
  {
    if(complete.cases(x2[i,])) # If test i is true
    {
      if(i==1) # if i=1 then this is group 1
      {
        group[i] <- 1 
      }
      if(i>1){ # if i>1   
        if(complete.cases(x2[i-1,])) # If test i-1 is true
        {
          group[i] <- group[i-1]  
        }else{
          group[i] <- max(group)+1
        }
      }
    }
  }
  x2$group <- group
  
  x2cc <- x2[complete.cases(x2),]
  for(group.i in 1:max(x2cc$group)){
    x2cc.group.i <-  x2cc[x2cc$group==group.i,] 
    polygon(x=c(rownames(x2cc.group.i),rev(rownames(x2cc.group.i))),y=c(x2cc.group.i[,2],rev(x2cc.group.i[,3])),col=col.shade,border=col.shade,lwd=2)
  }
  
  # draw median and confidence band lines
  matpoints(rownames(x2),x2[,c(1,2,3)],type="l",lwd=c(2,1,1),lty=c(2,1,1),col="black")
  # draw reference lines (e.g. base run)
  if(!missing(x.ref)&!missing(y.ref)){
    points(x.ref,y.ref,type="o",lty=1,lwd=3,pch=16,col="black")
  }
  legend("topright",legend=paste("error bands represent ",CIpct,"% CI",sep=""),bty="n",cex=0.75)
}

# plotBootPhase()
plotBootPhase <- function(x.Boot,y.Boot,x.ref,y.ref,
                          col=rgb(0,0,0,0.2),probs=c(0.05,0.95),
                          text.x.adj=c(0,0,0,0),
                          text.y.adj=c(0,0,0,0),
                          ...){
  quad1=round(length(x.Boot[x.Boot>1 & y.Boot>1])/length(x.Boot),3)*100
  quad2=round(length(x.Boot[x.Boot<1 & y.Boot>1])/length(x.Boot),3)*100
  quad3=round(length(x.Boot[x.Boot<1 & y.Boot<1])/length(x.Boot),3)*100
  quad4=round(length(x.Boot[x.Boot>1 & y.Boot<1])/length(x.Boot),3)*100
  
  plot(x.Boot, y.Boot,
       xlim=range(c(0,1,x.Boot)), ylim=range(c(0,1,y.Boot)),
       pch=16,col=col,...)
  abline(h=1,lty=2)
  abline(v=1, lty=2)
  points(rep(x.ref,2), quantile(y.Boot,probs=probs), type="l", lty=1, col="green4", lwd=2)      
  points(quantile(x.Boot,probs=probs),rep(y.ref,2), type="l", lty=1, col="green4", lwd=2)
  
  usr <- par("usr")
  
  text(x=c(rep((usr[1]+1)/2,2),rep((usr[2]+1)/2,2))+text.x.adj,
       y=rep(c((usr[3]+1)/2,(usr[4]+1)/2),2)+text.y.adj,
       labels=paste(c(quad3,quad2,quad4,quad1),"%",sep=""),col="blue")
}

## Compute mean squared errors (MSE) for indices
# x: should equal spp
compute.mse <- function(x){
  
  t.ser <- x$t.series
  t.ser[t.ser==-99999.00] <- NA
  t.ser.nm <- names(t.ser)
  
  comp.mats <- x$comp.mats
  comp.mats.nm <- names(comp.mats)
  
  # Indices
  U.nm.ob <- t.ser.nm[substr(t.ser.nm,1,2)=="U." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]
  U.nm.pr <- t.ser.nm[substr(t.ser.nm,1,2)=="U." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="pr"]
  U.nm <- gsub(".ob","",U.nm.ob)
  
  U.ob <- t.ser[,sort(U.nm.ob),drop=FALSE]
  U.pr <- t.ser[,sort(U.nm.pr),drop=FALSE]
  
  sq.errors <- (U.ob-U.pr)^2
  
  U.mse <- apply(sq.errors,2,function(x){mean(x,na.rm=TRUE)})
  names(U.mse) <- gsub(".ob","",names(U.mse))
  
  # Landings
  L.nm.ob <- t.ser.nm[substr(t.ser.nm,1,2)=="L." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]
  L.nm.pr <- t.ser.nm[substr(t.ser.nm,1,2)=="L." & substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="pr"]
  
  L.ob <- t.ser[,sort(L.nm.ob),drop=FALSE]
  L.pr <- t.ser[,sort(L.nm.pr),drop=FALSE]
  
  sq.errors <- (L.ob-L.pr)^2
  
  L.mse <- apply(sq.errors,2,function(x){mean(x,na.rm=TRUE)})
  names(L.mse) <- gsub(".ob","",names(L.mse))
  
  # comps
  comp.nm.ob <- comp.mats.nm[substr(comp.mats.nm,1,6)%in%c("acomp.","lcomp.") & substr(comp.mats.nm,nchar(comp.mats.nm)-1,nchar(comp.mats.nm))=="ob"]
  comp.nm.pr <- comp.mats.nm[substr(comp.mats.nm,1,6)%in%c("acomp.","lcomp.") & substr(comp.mats.nm,nchar(comp.mats.nm)-1,nchar(comp.mats.nm))=="pr"]
  
  comp.ob <- comp.mats[sort(comp.nm.ob)]
  comp.pr <- comp.mats[sort(comp.nm.pr)]
  
  comp.mse <- rep(NA,length(comp.nm.ob))
  names(comp.mse) <- gsub(".ob","",comp.nm.ob)
  
  for(i in 1:length(comp.nm.ob)){
    comp.mse[i] <- mean((comp.ob[[i]]-comp.pr[[i]])^2)
  }
  
  out <- c(U.mse,L.mse,comp.mse)
  names(out) <- paste("mse",names(out),sep=".")
  
  return(out)
}

# Check to see if any estimated parameters are near bounds
# x: BAM rdat object (i.e. spp)
# p.cutoff: set minimum acceptable value for how close to a bound the estimated parameter can be
#           before considering it too close to a bound, as a proportion of the range of bounds
#            supplied to BAM
checkBounds <- function(x,p.cutoff=0.01){
  pc <- spp$parm.cons # parameter matrix
  pc <- pc[,pc[4,]>0] # only consider parameters that are being estimated (not fixed)
  pc.range <- pc[3,]-pc[2,] # parameter range (width of bounds)
  pc.mindist <- apply(abs(pc[2:3,,drop=FALSE]-pc[c(8,8),,drop=FALSE]),2,min) # distance from estimate to nearest bound
  pc.pdist <- pc.mindist/pc.range # distance to nearest bound as a proportion of bound range
  pc.nearbound <- pc.pdist<=p.cutoff
  return(list("p.cutoff"=p.cutoff,"pdist"=pc.pdist,"nearbound"=pc.nearbound))
}

