#' Convert BAM .cor file to csv
#'
#' Convert BAM .cor file to .csv file. Also invisibly returns the cor matrix.
#' @param cor_file 	cor file path
#' @param rdat 	BAM output rdat (list) object read in with dget(). Optional, but the function can do more if you supply the rdat
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' }
#'
cor2csv <- function(cor_file,rdat=NULL){
library(stringr)

# dir.cor.figs <- "cor-figs"
# if(!dir.exists(dir.cor.figs)){dir.create(dir.cor.figs)}
# if(!dir.exists(file.path(dir.cor.figs,"cor"))){dir.create(file.path(dir.cor.figs,"cor"))}
# if(!dir.exists(file.path(dir.cor.figs,"sd"))){dir.create(file.path(dir.cor.figs,"sd"))}
# if(!dir.exists(file.path(dir.cor.figs,"series"))){dir.create(file.path(dir.cor.figs,"series"))}

### cor file
# Read in ADMB cor file and create usable matrix
D.BAM.cor <- local({
  a <- readLines(cor_file)
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

if(!is.null(rdat)){
# Match dev vectors with useful labels when available (years or ages)
D.BAM.cor$dev.x <- local({
  parname <- as.character(D.BAM.cor$name)
  years <- rdat$t.series$year

  dev.x <- rep(NA,nrow(D.BAM.cor))

  # Assign age to Nage_dev
  dev.x[grepl(pattern="Nage_dev",x=parname)] <- paste("age",sprintf("%02.0f",rdat$a.series$age[-1]),sep="")

  # Assign year to rec_dev (assumes continuous vector of rec_dev)
  years.rec_dev <- with(rdat$parm.tvec,{year[min(which(!is.na(log.rec.dev))):max(which(!is.na(log.rec.dev)))]})
  dev.x[grepl(pattern="rec_dev",x=parname)] <- paste(years.rec_dev)

  # Assign year to F_dev vectors
  F_dev.names <- unique(parname[grepl(pattern="log_F_dev_",x=parname)])
  t.series <- rdat$t.series
  ts.names <- names(t.series)
  for(i in F_dev.names){
    root.i <- gsub(pattern="log_F_dev_",replacement = "",x=i)
    if(!grepl(pattern="D_",x=root.i)){
      match.rdat.i <- paste("L",gsub(pattern="L_",replacement="",x=root.i),"ob",sep=".")
      #not.match.rdat.i <- paste("F",gsub(pattern="L_",replacement="",x=root.i),"D",sep=".")
      name.ts.i <- ts.names[grepl(pattern=match.rdat.i,x=ts.names)]#&(!grepl(pattern=not.match.rdat.i,x=ts.names))]
    }
    if(grepl(pattern="D_",x=root.i)){
      match.rdat.i <- paste("D",gsub(pattern="D_",replacement="",x=root.i),"ob",sep=".")
      name.ts.i <- ts.names[grepl(pattern=match.rdat.i,x=ts.names)]
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
  a <- D.BAM.cor[grepl(pattern="rec_dev|dev_rec",x=parname),c("index","name","value","std.dev","CV")]
  rownames(a) <- with(rdat$parm.tvec,{year[min(which(!is.na(log.rec.dev))):max(which(!is.na(log.rec.dev)))]})#rdat$t.series$year[which(rdat$t.series$logR.dev!=0)]
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
write.csv(D.BAM.cor.sub,file="BAM.cor.sub.csv",row.names = TRUE)
} # end if(!is.null(rdat))


write.csv(D.BAM.cor,file="BAM.cor.csv",row.names = FALSE)

# dum <-  latex(D.BAM.cor.sub,
#               file="tab-admb-par.tex",where="!htbp",
#               caption=paste("Names and estimated values of parameters estimated in the base run of the Beaufort Assessment Model." ,sep=""),
#               caption.lot="Parameter estimates from the BAM base run",
#               label="tab-admb-par",booktabs=TRUE,size="small",landscape=FALSE,
#               colheads =c("Parameter","Value"),
#               longtable=TRUE,
#               rowlabel="ID",
#               col.just=rep("r",ncol(D.BAM.cor.sub)),blank.dot=TRUE)

if(is.null(rdat)){
invisible(list("D.BAM.cor"=D.BAM.cor))
}else{
invisible(list("D.BAM.cor"=D.BAM.cor,"D.BAM.cor.rec_dev"=D.BAM.cor.rec_dev,"M.BAM.cor"=M.BAM.cor,"M.BAM.cor.dev"=M.BAM.cor.dev,"M.BAM.cor.nodev"=M.BAM.cor.nodev))
}
}
