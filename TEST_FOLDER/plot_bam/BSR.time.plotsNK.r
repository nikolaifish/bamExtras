BSR.time.plotsNK <- function (x, DataName = deparse(substitute(x)), draft = TRUE,
          start.drop = 0, graphics.type = NULL, use.color = TRUE,
          units.b = x$info$units.biomass, units.ssb = x$info$units.ssb,
          units.r = x$info$units.rec, legend.pos = "topright", from.zero = TRUE,
          BSR.references = list("Bmsy", "SSBmsy", "Rmsy"))
{
  if (!("t.series" %in% names(x))) {
    Errstring = (paste("Component ", deparse(substitute(x)),
                       "$t.series not found.", sep = ""))
    warning(Errstring, immediate. = TRUE)
    return(invisible(-1))
  }
  if (!("parms" %in% names(x))) {
    Errstring = (paste("Component ", deparse(substitute(x)),
                       "$parms not found.", sep = ""))
    warning(Errstring, immediate. = TRUE)
    return(invisible(-1))
  }
  h.ref = ifelse(is.list(BSR.references), TRUE, FALSE)
  ref.parms = rep(FALSE, length = 3)
  if (h.ref) {
    ref.parms = BSR.references %in% names(x$parms)
    for (i in 1:3) {
      if ((!is.null(BSR.references[[i]])) && (!ref.parms[i])) {
        Errstring = (paste("BSR.reference element ",
                           BSR.references[[i]], " not found.", sep = ""))
        warning(Errstring, immediate. = TRUE)
        return(invisible(-1))
      }
    }
  }
  plot.options = FGGetOptions()
  ts <- x$t.series[(start.drop + 1):nrow(x$t.series), ]
  parms <- x$parms
  lab.x <- "Year"
  savepar <- FGSetPar(draft)
  PlotTitle = ""
  if (length(graphics.type > 0)) {
    write.graphs <- TRUE
    GraphicsDirName <- paste(DataName, "-figs/BSR", sep = "")
    GraphicsDirName<<-GraphicsDirName
  }
  else {
    write.graphs <- FALSE
  }

  if ("B" %in% names(ts)) {
    lab.y <- FGMakeLabel("Total biomass", units.b)
    ylim <- range(ts$B, na.rm = TRUE)
    if (from.zero)
      ylim <- range(ylim, 0)
    if ((!is.null(BSR.references[[1]])) && h.ref) {
      hrefstring <- BSR.references[[1]]
      hrefindex <- which(names(x$parms) == BSR.references[[1]])
      href <- unlist(x$parms[hrefindex])
      hrefnames <- BSR.references[[1]]
      ylim <- range(ylim, href)
    }
    else {
      href <- NULL
      hrefnames <- NULL
    }
    if (draft)
      PlotTitle <- FGMakeTitle("Biomass", DataName)
    FGTimePlot(ts$year, ts$B, lab.x = lab.x, lab.y = lab.y,
               href = href, hrefnames = hrefnames, use.color = use.color,
               ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "B.total",
                 graphics.type)
  }
  if ("B" %in% names(ts) && ref.parms[1]) {
    hrefstring <- BSR.references[[1]]
    hrefindex <- which(names(x$parms) == BSR.references[[1]])
    Bref <- unlist(x$parms[hrefindex])
    ylim <- range(ts$B/Bref, 1.1, na.rm = TRUE)
    lab.y <- paste("B / ", hrefstring, sep = "")
    if (from.zero)
      ylim <- range(ylim, 0)
    if (draft)
      PlotTitle <- FGMakeTitle(lab.y, DataName)
    href <- if (h.ref) {
      1
    }
    else {
      NULL
    }
    FGTimePlot(ts$year, ts$B/Bref, lab.x = lab.x, lab.y = lab.y,
               href = href, hrefnames = NULL, use.color = use.color,
               ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "B.Bref",
                 graphics.type)
  }
  if ("B" %in% names(ts) && "B0" %in% names(parms)) {
    lab.y <- "B / B0"
    if (ref.parms[1] && h.ref) {
      href <- Bref/parms$B0
      hrefnames <- paste(BSR.references[[1]], " / B0",
                         sep = "")
      ylim <- range(ts$B/parms$B0, Bref/parms$B0, 1, na.rm = TRUE)
    }
    else {
      href <- NULL
      hrefnames <- NULL
      ylim <- range(ts$B/parms$B0, 1, na.rm = TRUE)
    }
    if (from.zero)
      ylim <- range(ylim, 0)
    if (draft)
      PlotTitle <- FGMakeTitle("B/B0", DataName)
    FGTimePlot(ts$year, ts$B/parms$B0, lab.x = lab.x, lab.y = lab.y,
               href = href, hrefnames = hrefnames, use.color = use.color,
               ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "B.B0",
                 graphics.type)
  }
  if ("SSB" %in% names(ts)) {
    lab.y <- FGMakeLabel("Spawning stock", units.ssb)
    ylim <- range(ts$SSB, na.rm = TRUE)
    if (from.zero)
      ylim <- range(ylim, 0)
    href <- as.numeric(c(NA, NA))
    hrefnames <- as.character(c(NA, NA))
    if ((!is.null(BSR.references[[2]])) && h.ref) {
      hrefstring <- BSR.references[[2]]
      hrefindex <- which(names(x$parms) == BSR.references[[2]])
      href[1] <- unlist(x$parms[hrefindex])
      hrefnames[1] <- BSR.references[[2]]
      ylim <- range(ylim, href[1])
      if ("msst" %in% names(parms)) {
        href[2] <- parms$msst
        hrefnames[2] <- "M*S*S*T"
        ylim <- range(ylim, href[2])
      }
      href <- href[!is.na(href)]
      hrefnames <- parse(text = hrefnames[!is.na(hrefnames)])
    }
    else {
      href <- NULL
      hrefnames <- NULL
    }
    if (draft)
      PlotTitle <- FGMakeTitle("Spawning biomass", DataName)
    FGTimePlot(ts$year, ts$SSB, lab.x = lab.x, lab.y = lab.y,
               href = href, hrefnames = hrefnames, legend.pos = legend.pos,
               ylim = ylim, use.color = use.color, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "SSB",
                 graphics.type)
  }
  if ("SSB" %in% names(ts) && ref.parms[2]) {
    hrefstring <- BSR.references[[2]]
    hrefindex <- which(names(x$parms) == BSR.references[[2]])
    Sref <- unlist(x$parms[hrefindex])
    ylim <- range(ts$SSB/Sref, 1.2, na.rm = TRUE)
    lab.y <- paste("SSB / ", hrefstring, sep = "")
    if (from.zero)
      ylim <- range(ylim, 0)
    if (draft)
      PlotTitle <- FGMakeTitle(lab.y, DataName)
    href <- if (h.ref) {
      1
    }
    else {
      NULL
    }
    FGTimePlot(ts$year, ts$SSB/Sref, lab.x = lab.x, lab.y = lab.y,
               href = href, hrefnames = NULL, use.color = use.color,
               ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "SSB.SSBref",
                 graphics.type)
  }
  if ("SSB" %in% names(ts) && "SSB0" %in% names(parms)) {
    lab.y <- "SSB / SSB0"
    if (ref.parms[2] && h.ref) {
      href <- Sref/parms$SSB0
      hrefnames <- paste(BSR.references[[2]], " / SSB0",
                         sep = "")
      ylim <- range(ts$SSB/parms$SSB0, Sref/parms$SSB0,
                    1, na.rm = TRUE)
    }
    else {
      href <- NULL
      hrefnames <- NULL
      ylim <- range(ts$SSB/parms$SSB0, 1, na.rm = TRUE)
    }
    if (from.zero)
      ylim <- range(ylim, 0)
    if (draft)
      PlotTitle <- FGMakeTitle("SSB/SSB0", DataName)
    FGTimePlot(ts$year, ts$SSB/parms$SSB0, lab.x = lab.x,
               lab.y = lab.y, href = href, hrefnames = hrefnames,
               use.color = use.color, ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "SSB.SSB0",
                 graphics.type)
  }
  if ("recruits" %in% names(ts)) {
    lab.y <- FGMakeLabel("Recruitment", units.r)
    ylim <- range(ts$recruits, na.rm = TRUE)
    if (from.zero)
      ylim <- range(ylim, 0)
    if ((!is.null(BSR.references[[3]])) && h.ref) {
      hrefstring <- BSR.references[[3]]
      hrefindex <- which(names(x$parms) == BSR.references[[3]])
      href <- unlist(x$parms[hrefindex])
      hrefnames <- BSR.references[[3]]
      ylim <- range(ylim, href)
    }
    else {
      href <- NULL
      hrefnames <- NULL
    }
    if (draft)
      PlotTitle <- FGMakeTitle("Recruitment", DataName)
    FGTimePlot(ts$year, ts$recruits, lab.x = lab.x, lab.y = lab.y,
               href = href, hrefnames = hrefnames, use.color = use.color,
               ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "R",
                 graphics.type)
  }
  if ("recruits" %in% names(ts) && ref.parms[3]) {
    hrefstring <- BSR.references[[3]]
    hrefindex <- which(names(x$parms) == BSR.references[[3]])
    Rref <- unlist(x$parms[hrefindex])
    ylim <- range(ts$recruits/Rref, 1.2, na.rm = TRUE)
    lab.y <- paste("Recruits / ", hrefstring, sep = "")
    if (from.zero)
      ylim <- range(ylim, 0)
    if (draft)
      PlotTitle <- FGMakeTitle(lab.y, DataName)
    href <- if (h.ref) {
      1
    }
    else {
      NULL
    }
    FGTimePlot(ts$year, ts$recruits/Rref, lab.x = lab.x,
               lab.y = lab.y, href = href, hrefnames = NULL, use.color = use.color,
               ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "R.Rref",
                 graphics.type)
  }
  if ("recruits" %in% names(ts) && "R0" %in% names(parms)) {
    lab.y <- "Recruits / R0"
    if (ref.parms[3] && h.ref) {
      href <- Rref/parms$R0
      hrefnames <- paste(BSR.references[[3]], " / R0",
                         sep = "")
      ylim <- range(ts$recruits/parms$R0, Rref/parms$R0,
                    1, na.rm = TRUE)
    }
    else {
      href <- NULL
      hrefnames <- NULL
      ylim <- range(ts$recruits/parms$R0, 1, na.rm = TRUE)
    }
    if (from.zero)
      ylim <- range(ylim, 0)
    if (draft)
      PlotTitle <- FGMakeTitle("Recruits/R0", DataName)
    FGTimePlot(ts$year, ts$recruits/parms$R0, lab.x = lab.x,
               lab.y = lab.y, href = href, hrefnames = hrefnames,
               use.color = use.color, ylim = ylim, main = PlotTitle)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "R.R0",
                 graphics.type)
  }
  if ("logR.dev" %in% names(ts)) {
    a <- loess(logR.dev ~ year, data = ts, span = 0.6)
    ts$lRd.smooth <- rep(NA, nrow(ts))
    ts$lRd.smooth[!is.na(ts$logR.dev)] <- a$fitted
    ylim <- c(-1, 1) * max(abs(c(ts$logR.dev, ts$lRd.smooth)),
                           na.rm = TRUE)
    if (use.color)
      FGO <- plot.options$color
    else FGO <- plot.options$bw
    href <- if (h.ref) {
      0
    }
    else {
      NULL
    }
    if (draft)
      PlotTitle <- FGMakeTitle("Recruitment deviations",
                               DataName)
    FGTimePlot(x = ts$year, y2 = ts$logR.dev, y = ts$lRd.smooth,
               lab.x = "Year", Y1Col = FGO$clr.lightline, Y2Col = FGO$clr.line,
               href = href, use.color = use.color, lab.y = "log Recruitment deviations + loess",
               hrefnames = NULL, main = PlotTitle, FGtype = "linepointnodots",
               ylim = ylim)
    if (write.graphs)
      FGSavePlot(GraphicsDirName, DataName, GraphName = "R.logRdev",
                 graphics.type)
  }
  return(invisible(NULL))
}
