##################################################################################
###   Sample R program for testing and demonstrating the "FishGraph" collection
###     of R graphics functions for analysis of stock-assessment results.
##################################################################################

##### Start fresh ###########
#rm(list=ls(all=TRUE))
graphics.off()
.SavedPlots <- NULL
# If
figPaths <- file.path("spp-figs",list.files("spp-figs",recursive = TRUE)) # paths to FishGraph figures
if(length(figPaths)>0){file.remove(figPaths)} # If FishGraph figures exist from a previous run, remove them

library(FishGraph)
library(plot.matrix)
library(stringr)
library(corrplot)
library(Hmisc)

# Fish Graph beta functions
source("Comp.plots.beta.r")
source("Comp.yearly.plots.beta.R")
source("Index.plots.beta.r")
source("StockRec.plots.beta.R")

load("_05_BAM_ConvertToSpeciesSpecificData.RData")

##### FUNCTIONS #####

# setTsYears()
# Set years of time series in data frame, when rownames are years
setTsYears <- function(x,years){
  # Initialize data frame with desired dimensions
  x <- as.matrix(x)
  x.out <- matrix(NA,nrow=length(years),ncol=ncol(x),
                  byrow=TRUE,dimnames=list(years,colnames(x)))
  x.out[rownames(x)[rownames(x)%in%paste(years)],] <- x[rownames(x)[rownames(x)%in%paste(years)],]
  x.out
}

# Data input timeline
data.timeline.plot <- function(spp,
                               dataTypeAbbkey=c("L"="landings","D"="discards","U"="index","acomp"="ageComp","lcomp"="lengthComp"),
                               fleetAbbKey=c("cHl"="com handline", "cTw" = "com trawl",
                                             "rHb" = "rec headboat", "rGe"="rec MRIP",
                                             "sCT"="chevron trap"),
                               dataTypes = c("L","D","U","acomp","lcomp")
                               ){
  t.ser <- spp$t.series
  t.ser.nm <- names(t.ser)
  
  # Landings
  L.nm <- t.ser.nm[substr(t.ser.nm,1,2)=="L." & 
                        substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]

  # Discards
  D.nm <- t.ser.nm[substr(t.ser.nm,1,2)=="D." & 
                        substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]
  
  # Comps
  acomp.nm <- t.ser.nm[substr(t.ser.nm,1,5)=="acomp" &
                               substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))==".n"]
  lcomp.nm <- t.ser.nm[substr(t.ser.nm,1,5)=="lcomp" &
                                    substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))==".n"]
  
  # Indices
  U.nm <- t.ser.nm[substr(t.ser.nm,1,2)=="U." & 
                        substr(t.ser.nm,nchar(t.ser.nm)-1,nchar(t.ser.nm))=="ob"]
  
  # Combine data
  nm <- as.character(unlist(mget(paste(dataTypes,"nm",sep="."))))
  
  t.ser2 <- t.ser[,nm]
  t.ser2[t.ser2==-99999] <- NA
  t.ser2.logical <- !is.na(t.ser2)
  
  # data set where values are years
    t.ser2.years <- apply(t.ser2.logical,2,function(x){
    out <- rep(NA,length(x))
    out[which(x)] <- as.numeric(rownames(t.ser2.logical)[x])
    out
    })
  rownames(t.ser2.years) <- rownames(t.ser2)
  
  fleetAbb <- gsub(pattern="L.|D.|acomp.|lcomp.|U.|.ob|.n",replacement="",x=colnames(t.ser2.years))
  dataTypeAbb <- unlist(lapply(strsplit(colnames(t.ser2.years),split=".",fixed=TRUE),function(x){x[1]}))
  
  fleetName <- fleetAbbKey[fleetAbb]
  dataTypeName <- dataTypeAbbkey[dataTypeAbb]
  
  colnames(t.ser2.years) <- paste(dataTypeName,fleetName)
  
  # data set where values indicate cell color
  t.ser2.fleet <- t.ser2.years
  for(i in 1:ncol(t.ser2.fleet)){
    x <- t.ser2.fleet[,i]
    x[which(!is.na(x))] <- fleetName[i]
    t.ser2.fleet[,i] <- x
  }
  
# Plot using plot.matrix package
  # years <- as.numeric(rownames(t.ser2.years))
  # dataNames <- colnames(t.ser2.years)
  plot(t(t.ser2.fleet),xlab="",ylab="",col=rainbow,main="")


 # for(i in 1:ncol(t.ser2.years)){
 #   xi <- t.ser2.years[,i]
 #   yi <- rep(i,length(xi))
 #   
 #   xlim <- range(years)
 #   ylim <- c(0,ncol(t.ser2.years))
 #   if(i==1){
 #     plot(xi,yi,xlim=xlim,ylim=ylim,type="p",pch=15,cex=2,
 #          lwd=lwd,
 #          xaxt="n",yaxt="n",
 #          xlab="",ylab="")
 #   }else{
 #     points(xi,yi,type="p",lwd=lwd,pch=15)
 #   }
 #   axis(side=1,at=years,labels=years,las=2)
 #   axis(side=2,at=1:ncol(t.ser2.years),labels=dataNames,las=2)
 # }
  
}
F.stacked.NK <- function(x,start.drop = 0, legend.pos = "topleft", PlotTitle="",
                         Fnames.exclude=c("F.Fmsy","F.F30.ratio","F.full"),colvec){
  yrndx <- (start.drop + 1):nrow(x$t.series)
  Fcols <- grep("^F\\.", names(x$t.series))
  Fdata <- x$t.series[yrndx, Fcols, drop = FALSE]
  Fdata <- Fdata[,!names(Fdata)%in%Fnames.exclude]
  
  Fmat <- t(as.matrix(Fdata))
  bar.names <- sub("F\\.", "", rownames(Fmat))
  if(missing(colvec)){
    colvec <- FGGetPal(length(bar.names), T)
  }
  
  barplot(Fmat, beside = FALSE, axis.lty = 1, col = colvec, 
          main = PlotTitle, ylab = "Fishing mortality rate", 
          xlab = "Year", las = 1)
  legend(legend.pos, inset = c(0.05, 0.1), rev(bar.names), 
         fill = rev(colvec), bg = "white")
}

############

# Years to plot
plotYears <- 1972:2017

# tseries.plot()
text.legend <- env.Spe$text.legend
tseries.plot <- env.Spe$tseries.plot
plotLabel <- ""

##### Read in the data from the ASCII .rdat file: #####
dADMB <- paste(wd,"ADMBFiles",sep="/")
spp <- dget(file.path(dADMB,SperdatName))

# Run OY.figs
source("Analyses-OY.r")

SpecorName <- gsub(pattern=".rdat",replacement=".cor",x=SperdatName)
dir.cor.figs <- "cor-figs"
if(!dir.exists(dir.cor.figs)){dir.create(dir.cor.figs)}
if(!dir.exists(file.path(dir.cor.figs,"cor"))){dir.create(file.path(dir.cor.figs,"cor"))}
if(!dir.exists(file.path(dir.cor.figs,"sd"))){dir.create(file.path(dir.cor.figs,"sd"))}
if(!dir.exists(file.path(dir.cor.figs,"series"))){dir.create(file.path(dir.cor.figs,"series"))}

### cor file
# Read in ADMB cor file and create usable matrix
D.BAM.cor <- local({
  a <- readLines(file.path(dADMB,SpecorName))
  b <- str_squish(a) # remove leading and trailing whitespace, and reduce internal whitespace to single space
  c <- b[-1]
  cor.columnnames <- unlist(strsplit(c[1],split=" "))
  cor.rownames <- cor.columnnames[!cor.columnnames%in%c("index","name","value","std.dev")]
  M.cor <- matrix(NA,nrow=length(cor.rownames),ncol=length(cor.columnnames),dimnames=list(cor.rownames,cor.columnnames))
  for(i in 1:(length(c)-1)){
    vals.i <- unlist(strsplit(c[i+1],split=" "))
    M.cor[i,1:length(vals.i)] <- vals.i
  }
  out <- as.data.frame(M.cor)

  # Convert most columns to numeric values
  cor.columnnames.numeric <- cor.columnnames[cor.columnnames!="name"]
  for(col.i in cor.columnnames.numeric){
    out[,col.i] <- as.numeric(as.character(out[,col.i]))
  }
  out$CV <- out$std.dev/out$value
  out <- out[,c("index","name","value","std.dev","CV",paste(out$index))]
  out$name <- as.character(out$name)
  return(out)
  })

# Match dev vectors with useful labels when available (years or ages)
D.BAM.cor$dev.x <- local({
  parname <- as.character(D.BAM.cor$name)
  years <- spp$t.series$year
  
  dev.x <- rep(NA,nrow(D.BAM.cor))
  
  # Assign age to Nage_dev
  dev.x[grepl(pattern="Nage_dev",x=parname)] <- paste("age",sprintf("%02.0f",spp$a.series$age[-1]),sep="")

  # Assign year to rec_dev (assumes continuous vector of rec_dev)
  years.rec_dev <- years[min(which(spp$t.series$logR.dev!=0)):max(which(spp$t.series$logR.dev!=0))]
  dev.x[grepl(pattern="rec_dev",x=parname)] <- paste(years.rec_dev)
  
  # Assign year to F_dev vectors
  F_dev.names <- unique(parname[grepl(pattern="log_F_dev_",x=parname)])
  t.series <- spp$t.series
  ts.names <- names(t.series)
  for(i in F_dev.names){
    root.i <- gsub(pattern="log_F_dev_",replacement = "",x=i)
    if(grepl(pattern="L_",x=root.i)){
      match.spp.i <- paste("L",gsub(pattern="L_",replacement="",x=root.i),"ob",sep=".")
      #not.match.spp.i <- paste("F",gsub(pattern="L_",replacement="",x=root.i),"D",sep=".")
      name.ts.i <- ts.names[grepl(pattern=match.spp.i,x=ts.names)]#&(!grepl(pattern=not.match.spp.i,x=ts.names))]
    }
    if(grepl(pattern="D_",x=root.i)){
      match.spp.i <- paste("D",gsub(pattern="D_",replacement="",x=root.i),"ob",sep=".")
      name.ts.i <- ts.names[grepl(pattern=match.spp.i,x=ts.names)]
    }
    years.i <- years[min(which(t.series[,name.ts.i]!=0)):max(which(t.series[,name.ts.i]!=0))]
    dev.x[grepl(pattern=i,x=parname)] <- paste(years.i)
  }
  dev.x
})

D.BAM.cor$name2 <- local({
a <- D.BAM.cor
par.name <- a$name
dev.x <- a$dev.x
par.id <- par.name
par.id[!is.na(dev.x)] <- paste(par.name[!is.na(dev.x)],dev.x[!is.na(dev.x)],sep=".")
par.id
})


D.BAM.cor.rec_dev <- local({
  parname <- as.character(D.BAM.cor$name)
  a <- D.BAM.cor[grepl(pattern="rec_dev",x=parname),c("index","name","value","std.dev","CV")]
  rownames(a) <- spp$t.series$year[which(spp$t.series$logR.dev!=0)]
  a
})

M.BAM.cor <- local({
  a <- D.BAM.cor
  par.name2 <- a$name2
  par.index <- a$index
  b <- as.matrix(a[,paste(par.index)]) # Select only columns with correlation values

  dimnames(b) <- list(par.name2,par.name2)
  b
})

M.BAM.cor.dev <- M.BAM.cor[grepl(pattern="dev",x=rownames(M.BAM.cor)),grepl(pattern="dev",x=colnames(M.BAM.cor))]
M.BAM.cor.nodev <- M.BAM.cor[!grepl(pattern="dev",x=rownames(M.BAM.cor)),!grepl(pattern="dev",x=colnames(M.BAM.cor))]

dev.types <- unique(D.BAM.cor$name[grepl(pattern="dev",x=D.BAM.cor$name)])
for(dev.type.i in dev.types){
M.out <-   M.BAM.cor[grepl(pattern=dev.type.i,x=rownames(M.BAM.cor)),grepl(pattern=dev.type.i,x=colnames(M.BAM.cor))]
rownames(M.out) <- gsub(pattern=paste(dev.type.i,".",sep=""),replacement="",x=rownames(M.out))
colnames(M.out) <- gsub(pattern=paste(dev.type.i,".",sep=""),replacement="",x=colnames(M.out))
assign(x=paste("M.BAM.cor",dev.type.i,sep="."),value=M.out)
}

D.BAM.cor.sub <- D.BAM.cor[,c("name2","value")]
D.BAM.cor.sub$name2 <- gsub(pattern="_",replacement=".",x=D.BAM.cor.sub$name2)

write.csv(D.BAM.cor,file="BAM.cor.csv",row.names = FALSE)
write.csv(D.BAM.cor.sub,file="BAM.cor.sub.csv",row.names = TRUE)

dum <-  latex(D.BAM.cor.sub,
            file="tab-admb-par.tex",where="!htbp",
            caption=paste("Names and estimated values of parameters estimated in the base run of the Beaufort Assessment Model." ,sep=""),      
            caption.lot="Parameter estimates from the BAM base run",      
            label="tab-admb-par",booktabs=TRUE,size="small",landscape=FALSE,
            colheads =c("Parameter","Value"),
            longtable=TRUE,
            rowlabel="ID",
            col.just=rep("r",ncol(D.BAM.cor.sub)),blank.dot=TRUE)

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
spp$t.series <- as.data.frame.matrix(setTsYears(spp$t.series,plotYears))

spp$CLD.est.mats <- lapply(spp$CLD.est.mats,FUN=setTsYears,years=plotYears)

spp$N.age <- setTsYears(spp$N.age,plotYears)
spp$N.age.mdyr <- setTsYears(spp$N.age.mdyr,plotYears)
spp$N.age.spawn <- setTsYears(spp$N.age.spawn,plotYears)
spp$B.age <- setTsYears(spp$B.age,plotYears)

# Truncate F-range in eq.series to make nicer plots
spp$eq.series <- spp$eq.series[which(spp$eq.series$F.eq<=1),]

# Calculate values to add to (CLD) plots
spp$parms$Dmsy.num <- spp$parms$Dmsy.knum*1000 # Add value of Dmsy in numbers
spp$parms$msy.num <- spp$parms$msy.knum*1000 # Add value of msy in numbers

##### Common arguments in FishGraph functions. Convenient to define them once.
ptype="pdf" #NULL (no quotes) for no plots saved; other options: "pdf", "wmf", "eps"
dtype=FALSE  #draft type
ctype=TRUE  #color type
########## Open a graphics device ###########
windows(width = 8, height = 8, record = TRUE)
#pdf(file="RG11fits-web1.pdf", width = 8, height = 8)

#Round Effective sample sizes for more enjoyable comp
spp$t.series[,grepl(pattern="comp",names(spp$t.series))&grepl(pattern="neff",names(spp$t.series))] <-
  round(spp$t.series[,grepl(pattern="comp",names(spp$t.series))&grepl(pattern="neff",names(spp$t.series))],2)

# spp$comp.mats$lcomp.cTw.ob <- matrix(rep(spp$comp.mats$lcomp.cTw.ob,2),
#                                      nrow=2,
#                                      ncol=ncol(spp$comp.mats$lcomp.cTw.ob),
#                                      dimnames=list(rep(rownames(spp$comp.mats$lcomp.cTw.ob),2),colnames(spp$comp.mats$lcomp.cTw.ob)),
#                                      byrow=TRUE)
# spp$comp.mats$lcomp.cTw.pr <- matrix(rep(spp$comp.mats$lcomp.cTw.pr,2),
#                                      nrow=2,
#                                      ncol=ncol(spp$comp.mats$lcomp.cTw.pr),
#                                      dimnames=list(rep(rownames(spp$comp.mats$lcomp.cTw.pr),2),colnames(spp$comp.mats$lcomp.cTw.pr)),
#                                      byrow=TRUE)

########## Call the functions #########################
### boundvecs ###
Bound.vec.plots(spp, draft=dtype, graphics.type=ptype)

### BSR ###
BSR.time.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, legend.pos="top",
               BSR.references=list("Bmsy", "SSBmsy", "Rmsy"))

windows(width = 7, height = 5, record = TRUE)
### CLD ###
CLD.total.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
                units.CLD.w="1000 lb whole",
                CLD.w.references=list("msy.klb", "msy.klb", "Dmsy.klb"),
                CLD.n.references=list("msy.num","msy.num","Dmsy.num"),
                plot.proportion = TRUE)

### comp ###
windows(width = 10, height = 8, record = TRUE)
Comp.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, c.min=0.2)

### compyr ###
windows(width = 10, height = 8, record = TRUE)
Comp.yearly.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, plot.neff=FALSE, print.neff=TRUE,
                  compact = TRUE, print.n=TRUE)
Cohort.plots(spp,draft=dtype, graphics.type=ptype)

### EQ ###
windows(width = 7, height = 5, record = TRUE)
Eq.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
         F.references=list("Fmsy"), user.Eq=list("L.eq.wholeklb", "L.eq.knum", "D.eq.knum"))

### F ###
windows(width = 7, height = 5, record = TRUE)
F.time.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
             F.references=list("Fmsy"="Fmsy"), F.additional=c("F.F30.ratio"="F.F30.ratio"))

pdf("spp-figs/F/spp.F.stacked.Fmsy.pdf")
F.stacked.NK(spp)
abline(h=spp$parms$Fmsy,lty=2,lwd=2)
text(x=par("usr")[1]+diff(par("usr")[1:2])*.1,y=spp$parms$Fmsy,labels="Fmsy",pos=3)
dev.off()

### growth ###
windows(width = 7, height = 5, record = TRUE)
Growth.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, plot.all = TRUE)

### index ###
windows(width = 7, height = 5, record = TRUE)
Index.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, resid.analysis=F)

### L ###
windows(width = 7, height = 5, record = TRUE)
Landings.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
               L.units=c("1000 lb","1000 lb","1000 fish","1000 fish"), D.units="1000 dead fish")

### NFZ ###
windows(width = 7, height = 5, record = TRUE)
NFZ.age.plots(spp,draft=dtype, use.color=ctype, graphics.type=ptype,
              user.plots="N.age.mdyr", start.drop=0, plot.CLD=F)

### parms ###
windows(width = 8, height = 6, record = TRUE)
Parm.plots(spp, graphics.type=ptype)

### phase ###
windows(width = 7, height = 5, record = TRUE)
Phase.plots(spp, start.drop=0, draft=dtype, use.color=ctype, graphics.type=ptype, Xaxis.F=F, year.pos=3,F.B.references=list("Fmsy","msst"))

### PR ###
windows(width = 7, height = 5, record = TRUE)
PerRec.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
             user.PR = list("SPR", "ypr.lb.whole"),
             F.references=list("Fmsy"))
### sel ###
windows(width = 7, height = 5, record = TRUE, xpos = 10, ypos = 10)
Selectivity.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, plot.points=T, compact=TRUE)

### SR ###
StockRec.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
               draw.lowess = FALSE, start.drop = 0, units.rec="number age-1 fish")

graphics.off()

### cor file plots ###
pdf(file.path(dir.cor.figs,"cor","par.pdf"),width=60,height=66)
par(oma=c(1,2,2,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.01)
corrplot(M.BAM.cor,method="circle",type="lower",
         mar=c(0,0,2,0),main="par")
dev.off()

pdf(file.path(dir.cor.figs,"cor","par.nodev.pdf"),width=10,height=11)
par(oma=c(1,2,2,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.01)
corrplot(M.BAM.cor.nodev,method="circle",type="lower",
         mar=c(0,0,2,0),main="par.nodev")
dev.off()

M.BAM.cor.dev.names <- paste("M.BAM.cor",dev.types,sep=".")
for(M.name.i in M.BAM.cor.dev.names){
  M.i <- get(M.name.i)
  pdf.name.i <- gsub(pattern="M.BAM.cor.",replacement="",M.name.i)
  pdf(file.path(dir.cor.figs,"cor",paste(pdf.name.i,"pdf",sep=".")),width=10,height=11)
  par(oma=c(0,2,2,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
  corrplot(M.i,method="circle",type="lower",
           mar=c(0,0,2,0),main=pdf.name.i)
  dev.off()
}

### std.dev plots
pdf(file.path(dir.cor.figs,"sd","par.nodev.pdf"),width=6,height=6)
M.i <- M.BAM.cor.nodev  
par(mfrow=c(1,1),mar=c(10,3,2,2),oma=c(0,0,0,0),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.005,lend="butt")
var.names <- rownames(M.i)
std.dev <- D.BAM.cor[D.BAM.cor$name2%in%var.names,"std.dev"]
value <- D.BAM.cor[D.BAM.cor$name2%in%var.names,"value"]
#barplot(height=std.dev,names.arg = var.names,las=2,ylab="std.dev",col="black")
x <- 1:length(std.dev)
plot(x=x,y=std.dev,ylab="std.dev",main="par.nodev",col="black",type="h",lwd=10,
     xaxt="n",xlab="")
axis(side=1,at=x,labels=var.names,las=2)
grid()
par(new=TRUE)
plot(x=x,value,type="p",pch=16,col="red",
     xaxt="n",yaxt="n",xlab="",ylab="")
axis(side=4,at=pretty(value),col="red",col.axis="red")
mtext(text="value",side = 4,line=par("mgp")[1],col="red")
dev.off()

M.BAM.cor.dev.names <- paste("M.BAM.cor",dev.types,sep=".")
for(M.name.i in M.BAM.cor.dev.names){
  M.name.root.i <- gsub(pattern="M.BAM.cor.",replacement = "",x=M.name.i)  
  pdf(file.path(dir.cor.figs,"sd",paste(M.name.root.i,"pdf",sep=".")),width=8,height=8)  
  M.i <- get(M.name.i)
  par(mfrow=c(1,1),mar=c(5,3,2,2),oma=c(0,0,0,0),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.005,lend="butt")
  var.names <- rownames(M.i)
  name2 <- paste(M.name.root.i,var.names,sep=".")
  std.dev <- D.BAM.cor[D.BAM.cor$name2%in%name2,"std.dev"]
  value <- D.BAM.cor[D.BAM.cor$name2%in%name2,"value"]
  bar.label <- D.BAM.cor[D.BAM.cor$name2%in%name2,"dev.x"]
  x <- 1:length(std.dev)
  plot(x=x,y=std.dev,ylab="std.dev",main=M.name.root.i,col="black",type="h",lwd=10,
       xaxt="n",xlab="")
  axis(side=1,at=x,labels=bar.label,las=2)
  grid()
  par(new=TRUE)
  plot(x=x,value,type="o",pch=16,col="red",
  xaxt="n",yaxt="n",xlab="",ylab="")
  axis(side=4,at=pretty(value),col="red",col.axis="red")
  mtext(text="value",side = 4,line=par("mgp")[1],col="red")
  
  dev.off()
}

# 
# par(mfrow=c(length(M.BAM.cor.names)/2,2),mar=c(2,2,0,0),mgp=c(1,.3,0),cex.lab=0.5,cex.axis=1,cex=0.5,tck=-0.01)
# for(M.name.i in M.BAM.cor.names){
#   M.i <- get(M.name.i)  
#   var.names <- rownames(M.i)
#   std.dev <- D.BAM.cor[D.BAM.cor$name2%in%var.names,"std.dev"]
#   barplot(height=std.dev,names.arg = var.names,las=2)
# }
# 
# M.BAM.cor.names <- c("M.BAM.cor.nodev",M.BAM.cor.dev.names)
# par(mfrow=c(length(M.BAM.cor.names)/2,2),mar=c(2,2,0,0),mgp=c(1,.3,0),cex.lab=0.5,cex.axis=1,cex=0.5,tck=-0.01)
# for(M.name.i in M.BAM.cor.names){
# M.i <- get(M.name.i)  
# var.names <- rownames(M.i)
# std.dev <- D.BAM.cor[D.BAM.cor$name2%in%var.names,"std.dev"]
# barplot(height=std.dev,names.arg = var.names,las=2)
# }

### t.series plots ###
pdf(file.path(dir.cor.figs,"series","rec_dev.pdf"),width=10,height=11)
tseries.plot(D=D.BAM.cor.rec_dev[,"value",drop=FALSE],
             DErr=D.BAM.cor.rec_dev[,"std.dev",drop=FALSE],
             xlab="year",ylab="rec_dev")
abline(h=0)
dev.off()

### Data timeline plot
pdf(file.path("data-timeline.pdf"),width=12,height=5)
par(mar=c(3,10,1,8),mgp=c(1,0.3,0),lend="butt",tck=0.01,las=2)
data.timeline.plot(spp)
dev.off()

pdf(file.path("data-timeline-comp.pdf"),width=12,height=2)
par(mar=c(3,10,1,8),mgp=c(1,0.3,0),lend="butt",tck=0.01,las=2)
data.timeline.plot(spp,dataTypes = c("acomp","lcomp"))
dev.off()
