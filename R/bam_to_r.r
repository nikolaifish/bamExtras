#' BAM inputs to R object conversion function
#'
#' This script reads in BAM dat, tpl, and cxx files, and converts them to R objects (character vectors of each line of the code). It also builds an R list on the tpl init objects, from the tpl and dat files.
#' @param dat.file dat file path
#' @param tpl.file tpl file path
#' @param cxx.file cxx file path
#' @param L.init.user user supplied L.init object. If missing, the function builds one from values in the dat and tpl files.
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' bam_to_r()

bam_to_r <- function(dat.file,tpl.file,cxx.file,L.init.user){
  # Read in dat, tpl, and cxx files
  dat <- readLines(con=dat.file)
  dat <- trimws(dat) # Remove leading and/or trailing whitespace from character strings
  tpl <- readLines(con=tpl.file)
  cxx <- readLines(con=cxx.file)

  ### tpl ###

  # Isolate tpl DATA_SECTION
  tplDATA <- tpl[which(tpl=="DATA_SECTION"):(which(tpl=="PARAMETER_SECTION")-1)]

  # Identify which items in tplDATA correspond to initialization of variables corresponding to values in the dat file
  tplDATA.items <- tplDATA[substr(tplDATA,start=1,stop=5)=="init_"]

  # Create a data frame that can be used to parse the dat file and create an R list
  D.tplVar <- local({
    # Split items into meaningful parts and assemble a data frame
    sp1 <- strsplit(tplDATA.items,split=";")

    sp2 <- lapply(sp1,FUN=function(x){x[[1]]})

    sp2ab <- strsplit(unlist(sp2),split=" ",fixed=TRUE)
    type <- unlist(lapply(sp2ab,FUN=function(x){x[[1]]}))
    nameArg <- unlist(lapply(sp2ab,FUN=function(x){x[[2]]}))
    sp.nameArg <- strsplit(nameArg,split="(",fixed=TRUE)
    name <- unlist(lapply(sp.nameArg,FUN=function(x){x[[1]]}))
    arg <- unlist(lapply(sp.nameArg,FUN=function(x){ifelse(length(x)==2,x[[2]],"")}))
    arg <- gsub(x=arg,pattern=")",replacement="",fixed=TRUE)

    sp3 <- unlist(lapply(sp1,FUN=function(x){ifelse(length(x)==2,x[[2]],"")}))
    comment <- unlist(lapply(strsplit(sp3,split="//",fixed=TRUE),FUN=function(x){ifelse(length(x)==2,x[[2]],"")}))

    # Create a data frame that can be used to parse the dat file and create an R list
    return(
      data.frame("type"=type,"name"=name,"arg"=arg,"comment"=comment,stringsAsFactors=FALSE)
    )
  })

  #tpl <- gsub(x=tpl,pattern=PastcxxName,replacement=ModcxxName)

  #### LnonInitDdatnum ####

  # Identify what type of information is in each row of the dat file (numeric or non-numeric)
  dat.num.rows.logical <- substr(dat,start="1",stop="1")%in%c(paste(0:9),"-",".")
  dat.num.rows <- which(dat.num.rows.logical)
  dat.nonNum.rows <- which(!dat.num.rows.logical) # Rows of dat file with non-numeric row values (i.e. comments and whitespace)

  # List of non-numeric row values from dat file
  L.nonInit <- local({
    a <- as.list(dat[dat.nonNum.rows])
    names(a) <- paste("datRow",sprintf("%05.0f", dat.nonNum.rows))
    a
  })

  # Create a data frame that includes only numeric values from the dat file that the tpl reads, and adjacent comments
  D.dat.num <- local({
    dat.num.row <- dat[dat.num.rows] # dat rows that begin with a numeric value

    sp1 <- strsplit(dat.num.row,split="#")
    dat.value <- unlist(lapply(sp1,FUN=function(x){x[[1]]}))
    dat.comment <- trimws(unlist(lapply(sp1,FUN=function(x){ifelse(length(x)==2,x[[2]],"")})))
    data.frame("value"=dat.value,"datRow"=dat.num.rows,"comment"=dat.comment,stringsAsFactors=FALSE)
  })

  #### Build L.init ####
  # Build R list with names equal to BAM tpl init object names and values equal to init objects from dat file
  initNames <- D.tplVar[,"name"]
  L.init <- vector("list", length(initNames))
  names(L.init) <- initNames
  datRowInit <- local({a <- rep(NA,length(L.init)) # Row of dat file associated with the first row of each init object
  names(a) <- names(L.init)
  return(a)
  })
  D.tplVar<<-D.tplVar
  for(i in 1:nrow(D.tplVar)){
    if(i==1){
      datRow.i <- i-1 # Initialize dat row counter
    }

    type.i <- D.tplVar[i,"type"]
    name.i <- D.tplVar[i,"name"]
    arg.i <- D.tplVar[i,"arg"]

    datRow.first <- datRow.i + 1 # dat row counter associated with the first row of init object i
    datRowInit[name.i] <- D.dat.num[datRow.first,"datRow"]

    if(type.i%in%c("init_int","init_number")){
      datRow.i <- datRow.first
      varValue.i <- D.dat.num$value[datRow.i]
      if(type.i=="init_int"){
        #varValue.i <- as.integer(varValue.i)
        varValue.i <- as.character(varValue.i)
      }
      if(type.i=="init_number"){
        #varValue.i <- as.numeric(varValue.i)
        varValue.i <- as.character(varValue.i)
      }
    }

    if(type.i%in%c("init_vector","init_ivector")){
      datRow.i <- datRow.first
      # Initialize vector
      vector.length <- local({
        out <- strsplit(arg.i,split=",",fixed=TRUE)[[1]][2]
        out <- type.convert(out,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(out)){out <- get(out)}
        return(as.numeric(out))
      })


      if(type.i=="init_vector"){
        varValue.i <- vector(mode = "numeric", length = vector.length)
      }
      if(type.i=="init_ivector"){
        varValue.i <- vector(mode = "integer", length = vector.length)
      }

      # Assign values to vector
      #varValue.i <- as.numeric(scan(text=D.dat.num$value[i],what=""))
      varValue.i <- as.character(strsplit(D.dat.num$value[datRow.i], "\\s+")[[1]])
    }

    if(type.i=="init_matrix"){
      datRow.first <- datRow.first # First row of matrix
      matrix.nrow <- local({
        rowStart <- strsplit(arg.i,split=",",fixed=TRUE)[[1]][1]
        rowStart <- type.convert(rowStart,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(rowStart)){rowStart <- get(rowStart)}

        rowEnd <- strsplit(arg.i,split=",",fixed=TRUE)[[1]][2]
        rowEnd <- type.convert(rowEnd,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(rowEnd)){rowEnd <- get(rowEnd)}

        nrow <- length(rowStart:rowEnd)

        return(as.numeric(nrow))
      })

      matrix.ncol <- local({
        colStart <- strsplit(arg.i,split=",",fixed=TRUE)[[1]][3]
        colStart <- type.convert(colStart,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(colStart)){colStart <- get(colStart)}

        colEnd <- strsplit(arg.i,split=",",fixed=TRUE)[[1]][4]
        colEnd <- type.convert(colEnd,as.is=TRUE) # Converts character vector to numeric if it can. Otherwise leaves it as character
        if(is.character(colEnd)){colEnd <- get(colEnd)}

        ncol <- length(colStart:colEnd)

        return(as.numeric(ncol))
      })

      datRow.last <- datRow.first+matrix.nrow-1 # last row of matrix

      varValue.i <- matrix(
        as.character(unlist(strsplit(D.dat.num$value[datRow.first:datRow.last], "\\s+"))),
        byrow=TRUE,
        nrow=matrix.nrow, ncol=matrix.ncol)
      datRow.i <- datRow.last
    }
    varValue.i <- trimws(varValue.i) # Remove extra whitespace from values

    assign(name.i,varValue.i) # Assign object to global environment (important for defining some later objects)
    L.init[[name.i]] <- varValue.i # Assign value to L.init
  }

  # If L.init.user is supplied, override L.init.
  if(!missing(L.init.user)){
    L.init <- L.init.user
  } # End "if(missing(L.init))"

  #### Create L.dat ####
  # Create new list of dat file elements incorporating L.init

  L.dat <- local({

    L.init2 <- L.init
    names(L.init2) <- paste("datRow",sprintf("%05.0f", datRowInit[names(L.init)]))

    L.dat <- c(L.init2,L.nonInit)

    L.dat <- L.dat[order(names(L.dat))]

    # convert all items to make one long character vector
    L.dat <- lapply(L.dat,FUN=function(x){
      if(class(x)=="matrix"){
        out <- apply(x,1,FUN=function(x){paste(x,collapse="\t")})
      }else{
        out <- paste(x,collapse="\t")
      }
    })
    return(L.dat)
  })
  # Add available comments back to rows with numeric values
  L.dat <- local({
    datRowsWithComments <- D.dat.num[which(D.dat.num$comment!=""),"datRow"]
    datRowsWithCommentsNames <- paste("datRow",sprintf("%05.0f", datRowsWithComments)) # numeric dat rows with comments
    datRowsWithCommentsComments <- D.dat.num$comment[D.dat.num$datRow%in%datRowsWithComments]
    L.dat[datRowsWithCommentsNames] <- paste(L.dat[datRowsWithCommentsNames],datRowsWithCommentsComments,sep=" # ")

    L.dat
  })

  return(list("dat"=unlist(L.dat),"tpl"=tpl,"cxx"=cxx,"L.init"=L.init,"L.dat"=L.dat))
}
