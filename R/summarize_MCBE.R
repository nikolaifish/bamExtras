#' Summarize results from MCBE uncertainty analysis
#'
#' @param dir_bam_sim Name of directory where run_MCBE results files have been saved, relative to the working directory.
#' @param obj_collect Names of objects in MCBE rdat files to collect and summarize
#' @param coresUse Number of cores to use when running parallel processes
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Run MCBE, writing files to dir_bam_sim
#' run_MCBE("GrayTriggerfish", dir_bam_sim="GrTr_sim", dir_bam_base="GrTr_base")
#' # Summarize results of MCBE and assign to object
#' MCBE_GrTr <- summarize_MCBE(dir_bam_sim="sim_GrTr")
#' }



summarize_MCBE <- function(dir_bam_sim="sim",
                           dir_bam_base="base",
                           obj_collect = c("info","parms","like","spr.brps","mse",
                                           "a.series","t.series","N.age","Z.age",
                                           "sel.age"),
                           coresUse=NULL
                     ){

# Test args
# fileName <-  "bam"
# dir_bam <- "base"
# dir_bam_sim="sim_GrTr"
# obj_collect = c("info","parms","like","spr.brps","a.series","t.series","sel.age","N.age","Z.age","mse")

#######################
  library(doParallel)
  library(foreach)

#### parallel setup
if(is.null(coresUse)){
coresAvail <- detectCores()
coresUse <- coresAvail-1
}
coresUse  <- min(c(coresUse,coresAvail))
cl <- makeCluster(coresUse)
registerDoParallel(cl)

### Do stuff with base run
# spp <- dget(file.path(dir_bam,paste(filename,'rdat', sep=".")))


# nm_sim <- sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim)
#
#
# ## Get base model stuff
#
# # Identify working directory
#   wd <- getwd()
#   message(paste("working directory:",wd))
#
simFiles <- list.files(dir_bam_sim)
simFilesRdat <- simFiles[grepl(x=simFiles,pattern=".rdat$")]
modRunName <- local({
  a <- gsub(x=simFilesRdat,pattern=".rdat",replacement = "",fixed=TRUE)
  num.parts <- unique(unlist(lapply(strsplit(a,split="-"),length))) # number of parts to file names separated by hyphens
  b <- matrix(unlist(strsplit(a,split="-")),ncol=num.parts,byrow=TRUE)
  c <- b[,1]
  if(num.parts==3){ # Give longer names, when sets of runs exist (e.g. profiling)
    c <- paste(b[,1],b[,3])
  }
  return(c)})

spp1 <- dget(paste(dir_bam_sim,simFilesRdat[1],sep="/"))

a.series.names <- names(spp1$a.series)
t.series.names <- names(spp1$t.series)
ages <- spp1$a.series$age
years <- spp1$t.series$year

# Retrieve values from simstrap rdat files
simSummary <- foreach(iter.j=seq_along(simFilesRdat),
                      .packages="bamExtras",
                       .multicombine=TRUE
) %dopar% {
  # dget each rdat file and pull stuff out of it
  rdat.j.name <- simFilesRdat[iter.j]
  modRunName.j <- modRunName[iter.j]#unlist(strsplit(rdat.j.name, "\\-|\\."))[1]
  modRunGroup.j <- unlist(strsplit(modRunName.j,split=" "))[2]
  rdat.j <- dget(paste(dir_bam_sim,rdat.j.name,sep="/"))

  ## Calculate some values needed in the benchmarks table
  eq <- rdat.j$eq.series
  Fmsy <- rdat.j$parms$Fmsy
  F.eq <- eq$F.eq
  #L.eq <- eq$L.eq.wholeklb
  L.eq <- eq$L.eq.gutklb
  # D.eq <- eq$D.eq.knum
  Fmsy.65 <- 0.65*Fmsy
  Fmsy.75 <- 0.75*Fmsy
  Fmsy.85 <- 0.85*Fmsy

  Lmsy.65.klb  <- L.eq[which.min(abs(F.eq-Fmsy.65))]
  Lmsy.75.klb  <- L.eq[which.min(abs(F.eq-Fmsy.75))]
  Lmsy.85.klb  <- L.eq[which.min(abs(F.eq-Fmsy.85))]
  # Dmsy.65.knum <- D.eq[which.min(abs(F.eq-Fmsy.65))]
  # Dmsy.75.knum <- D.eq[which.min(abs(F.eq-Fmsy.75))]
  # Dmsy.85.knum <- D.eq[which.min(abs(F.eq-Fmsy.85))]

  rdat.j$parms$Lmsy.65.klb <- Lmsy.65.klb
  rdat.j$parms$Lmsy.75.klb <- Lmsy.75.klb
  rdat.j$parms$Lmsy.85.klb <- Lmsy.85.klb
  # rdat.j$parms$Dmsy.65.knum <- Dmsy.65.knum
  # rdat.j$parms$Dmsy.75.knum <- Dmsy.75.knum
  # rdat.j$parms$Dmsy.85.knum <- Dmsy.85.knum

  # Compute mean squared errors associated with data fitted
  rdat.j$mse <- mse_calc(rdat.j)

  # Check for bounding issues
  rdat.j$parms$nearbound <-  any(check_bounds(rdat.j)$nearbound)

  # Create output list from each iteration of foreach
  # NOTE: This output will be combined with the .combine function
  out <- list("modRunName"=modRunName.j, "modRunGroup"=modRunGroup.j)
  out <- c(out,rdat.j[obj_collect])

  return(out)
}##### END FOREACH LOOP ######

# Redefine objects for collecting values from each simstrap run
names(simSummary) <- unlist(lapply(simSummary,FUN=function(x){x[["modRunName"]]}))

## Vectors

if("info"%in%obj_collect){
L.sim.info <- lapply(simSummary,FUN=function(x){x[["info"]]})    # Generate lists with potentially uneven levels
L.sim.info.itemnames <- unique(unlist(lapply(L.sim.info,names))) # Identify the set of item names among all levels in L.sim.info
# Fill in all item names in L.sim.info so that all levels have the same items in the same order
L.sim.info <- lapply(L.sim.info,FUN=function(x){
  missing.itemnames <- L.sim.info.itemnames[!L.sim.info.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.info.itemnames]
  a
})
L.sim.info <- lapply(1:length(L.sim.info),FUN=function(x){c("modRunName"=names(L.sim.info)[x],L.sim.info[[x]])})
D.sim.info <- do.call(rbind.data.frame, L.sim.info)
rownames(D.sim.info) <- D.sim.info$modRunName
out.info <- D.sim.info
}

if("parms"%in%obj_collect){
L.sim.parms <- lapply(simSummary,FUN=function(x){x[["parms"]]})
L.sim.parms.itemnames <- unique(unlist(lapply(L.sim.parms,names)))
L.sim.parms <- lapply(L.sim.parms,FUN=function(x){
  missing.itemnames <- L.sim.parms.itemnames[!L.sim.parms.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.parms.itemnames]
  a
})
L.sim.parms <- lapply(1:length(L.sim.parms),FUN=function(x){c("modRunName"=names(L.sim.parms)[x],L.sim.parms[[x]])})
D.sim.parms <- do.call(rbind.data.frame, L.sim.parms)
rownames(D.sim.parms) <- D.sim.parms$modRunName
out.parms <- D.sim.parms
}

if("like"%in%obj_collect){
L.sim.like <- lapply(simSummary,FUN=function(x){as.list(x[["like"]])})
L.sim.like.itemnames <- unique(unlist(lapply(L.sim.like,names)))
L.sim.like <- lapply(L.sim.like,FUN=function(x){
  missing.itemnames <- L.sim.like.itemnames[!L.sim.like.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.like.itemnames]
  a
})
L.sim.like <- lapply(1:length(L.sim.like),FUN=function(x){c("modRunName"=names(L.sim.like)[x],L.sim.like[[x]])})
D.sim.like <- do.call(rbind.data.frame, L.sim.like)
rownames(D.sim.like) <- D.sim.like$modRunName
out.like <- D.sim.like
}

if("spr.brps"%in%obj_collect){
L.sim.spr.brps <- lapply(simSummary,FUN=function(x){x[["spr.brps"]]})
L.sim.spr.brps.itemnames <- unique(unlist(lapply(L.sim.spr.brps,names)))
L.sim.spr.brps <- lapply(L.sim.spr.brps,FUN=function(x){
  missing.itemnames <- L.sim.spr.brps.itemnames[!L.sim.spr.brps.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.spr.brps.itemnames]
  a
})
L.sim.spr.brps <- lapply(1:length(L.sim.spr.brps),FUN=function(x){c("modRunName"=names(L.sim.spr.brps)[x],L.sim.spr.brps[[x]])})
D.sim.spr.brps <- do.call(rbind.data.frame, L.sim.spr.brps)
rownames(D.sim.spr.brps) <- D.sim.spr.brps$modRunName
out.spr.brps <- D.sim.spr.brps
}

if("mse"%in%obj_collect){
L.sim.mse <- lapply(simSummary,FUN=function(x){as.list(x[["mse"]])})
L.sim.mse.itemnames <- unique(unlist(lapply(L.sim.mse,names)))
L.sim.mse <- lapply(L.sim.mse,FUN=function(x){
  missing.itemnames <- L.sim.mse.itemnames[!L.sim.mse.itemnames%in%names(x)]
  missing.items <- rep(NA,length(missing.itemnames))
  names(missing.items) <- missing.itemnames
  a <- c(x,missing.items)
  a <- a[L.sim.mse.itemnames]
  a
})
L.sim.mse <- lapply(1:length(L.sim.mse),FUN=function(x){c("modRunName"=names(L.sim.mse)[x],L.sim.mse[[x]])})
D.sim.mse <- do.call(rbind.data.frame, L.sim.mse)
rownames(D.sim.mse) <- D.sim.mse$modRunName
out.mse <- D.sim.mse
}

## Matrices
if("a.series"%in%obj_collect){
# a.series
L.sim.a.series <- foreach(i=seq_along(a.series.names),
                      .multicombine=TRUE
) %dopar% {
  a <- lapply(simSummary,FUN=function(x){x[["a.series"]][[a.series.names[i]]]})
  b <- do.call(rbind.data.frame,a)
  names(b) <- ages # Rename columns
  return(b)
}
names(L.sim.a.series) <- a.series.names
out.a.series <- L.sim.a.series
}

if("t.series"%in%obj_collect){
# t.series
L.sim.t.series <- foreach(i=seq_along(t.series.names),
                          .multicombine=TRUE
) %dopar% {
  a <- lapply(simSummary,FUN=function(x){x[["t.series"]][[t.series.names[i]]]})
  b <- do.call(rbind.data.frame,a)
  names(b) <- years # Rename columns
  return(b)
}
names(L.sim.t.series) <- t.series.names
out.t.series <- L.sim.t.series
}

if("N.age"%in%obj_collect){
# N.age
N.age.dn1 <- dimnames(spp1$N.age)[[1]]
N.age.dn2 <- dimnames(spp1$N.age)[[2]]
L.sim.N.age <- foreach(i=seq_along(N.age.dn2),
                          .multicombine=TRUE
) %dopar% {
  a <- lapply(simSummary,FUN=function(x){x[["N.age"]][,N.age.dn2[i]]})
  b <- do.call(rbind.data.frame,a)
  names(b) <- N.age.dn1 # Rename columns
  return(b)
}
names(L.sim.N.age) <- N.age.dn2
out.N.age <- L.sim.N.age
}

if("Z.age"%in%obj_collect){
  # Z.age
  Z.age.dn1 <- dimnames(spp1$Z.age)[[1]]
  Z.age.dn2 <- dimnames(spp1$Z.age)[[2]]
  L.sim.Z.age <- foreach(i=seq_along(Z.age.dn2),
                         .multicombine=TRUE
  ) %dopar% {
    a <- lapply(simSummary,FUN=function(x){x[["Z.age"]][,Z.age.dn2[i]]})
    b <- do.call(rbind.data.frame,a)
    names(b) <- Z.age.dn1 # Rename columns
    return(b)
  }
  names(L.sim.Z.age) <- Z.age.dn2
  out.Z.age <- L.sim.Z.age
}

## Lists
if("sel.age"%in%obj_collect){
  # sel.age
  sel.age.names <- names(spp1$sel.age)
  L.sim.sel.age <- list()

  for(i in seq_along(sel.age.names)){
    if(is.null(dim(spp1$sel.age[[i]]))){
      # vectors
      a <- lapply(simSummary,FUN=function(x){x[["sel.age"]][[sel.age.names[i]]]})
      b <- do.call(rbind,a)
    }else{
      # matrices
      sel.age.dn1 <- dimnames(spp1$sel.age[[i]])[[1]]
      sel.age.dn2 <- dimnames(spp1$sel.age[[i]])[[2]]
      b <- foreach(j=seq_along(sel.age.dn2),
                             .multicombine=TRUE
      ) %dopar% {
        a <- lapply(simSummary,FUN=function(x){x[["sel.age"]][[i]][,sel.age.dn2[j]]})
        b <- do.call(rbind,a)
        return(b)
      }
      names(b) <- sel.age.dn2
    }
    L.sim.sel.age[[i]] <- b

  }
  names(L.sim.sel.age) <- sel.age.names
  out.sel.age <- L.sim.sel.age
}


# Assign out to object
out <- mget(paste0("out.",obj_collect))
names(out) <- gsub("^out.","",names(out))


# parallel shut down
  stopCluster(cl)

invisible(out)
}
