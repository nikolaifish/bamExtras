# Data input timeline
data.timeline.plot <- function(spp,
                               dataTypeAbbkey=c("L"="landings",#"D"="discards",
                                                "U"="index","acomp"="ageComp","lcomp"="lengthComp"),
                               fleetAbbKey=c("cH"="com handline", "cL" = "com longline",
                                             "rA" = "rec all",
                                             "sM"="MARMAP longline"),
                               dataTypes = c("L","U","acomp","lcomp")
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

  #fleetAbb <- gsub(pattern="L.|D.|acomp.|lcomp.|U.|.ob|.n",replacement="",x=colnames(t.ser2.years))
  fleetAbb <- unlist(strsplit(colnames(t.ser2.years),split=".",fixed=TRUE))[c(F,T,F)]
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
  plot(t(t.ser2.fleet),xlab="year",ylab="",col=rainbow,main="")
}
