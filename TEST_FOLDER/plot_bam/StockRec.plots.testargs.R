x=rdat_BlackSeaBass
DataName = "spp"
draft = FALSE
graphics.type = "pdf"
use.color = TRUE
start.drop = 0
units.ssb = x$info$units.ssb
units.rec = x$info$units.rec
rec.model = x$info$rec.model
draw.model = TRUE
draw.lowess = FALSE
draw.time = TRUE
year.pos = 1


if (!("t.series" %in% names(x))) {
  Errmsg <- paste("Data frame 't.series' not found in data object:",
                  deparse(substitute(x)))
  warning(Errmsg, immediate. = TRUE)
  return(invisible(-1))
}
if (!("parms" %in% names(x))) {
  Errmsg <- paste("List 'parms' not found in data object:",
                  deparse(substitute(x)))
  warning(Errmsg, immediate. = TRUE)
  return(invisible(-1))
}
if (!("rec.lag" %in% names(x$parms))) {
  Errmsg <- paste("Value 'parms$rec.lag' not found in data object:",
                  deparse(substitute(x)))
  warning(Errmsg, immediate. = TRUE)
  return(invisible(-1))
}
ts <- x$t.series
rec.lag <- x$parms$rec.lag
draw.bias <- FALSE
sndx = as.integer(rep(0, 2))
sndx[1] <- start.drop + 1
sndx[2] <- nrow(ts) - rec.lag
if (sndx[1] >= sndx[2]) {
  Errmsg <- "Length of S-R data series is <= 1"
  warning(Errmsg, immediate. = TRUE)
  return(invisible(-1))
}
sndx <- sndx[1]:sndx[2]
rndx <- sndx + x$parms$rec.lag
plot.options = FGGetOptions()
PlotTitle <- ""
savepar <- FGSetPar(draft)
if (!is.null(graphics.type)) {
  write.graphs <- TRUE
  GraphicsDirName <- paste(DataName, "-figs/SR", sep = "")
}else {write.graphs <- FALSE}
if (use.color){
  parlist <- plot.options$color
}else {parlist <- plot.options$bw}
clr.points <- parlist$clr.line
clr.lowess <- parlist$clr.lightline
ssb.max <- max(ts$SSB[sndx], na.rm = TRUE)
r.max <- max(ts$recruits[rndx], na.rm = TRUE)
lab.x <- FGMakeLabel("Spawning stock", units.ssb)
lab.y <- FGMakeLabel("Recruitment", units.rec)
rec.sim <- as.numeric(NA)
lim.x = range(0, 1.1 * ssb.max)
lim.y = range(0, 1.1 * r.max)
if (is.null(rec.model)) {
  draw.model <- FALSE
  rec.model <- "none"
  curve.OK <- TRUE
}
if (draw.model) {
  stock.sim <- seq(from = 0.001 * ssb.max, to = 1.1 *
                     ssb.max, length = 200)
  curve.OK <- FALSE
  if (rec.model == "BH") {
    alpha <- x$parms$BH.alpha
    beta <- x$parms$BH.beta
    bias.corr <- x$parms$BH.biascorr
    rec.sim <- alpha * stock.sim/(1 + beta * stock.sim)
    curve.OK <- TRUE
  }
  if (rec.model == "BH-steep") {
    R0 <- x$parms$BH.R0
    h <- x$parms$BH.steep
    phi0 <- x$parms$BH.Phi0
    bias.corr <- x$parms$BH.biascorr
    rec.sim <- (0.8 * R0 * h * stock.sim)/(0.2 * phi0 *
                                             R0 * (1 - h) + (h - 0.2) * stock.sim)
    curve.OK <- TRUE
  }
  if (rec.model == "Ricker") {
    alpha <- x$parms$Ricker.alpha
    beta <- x$parms$Ricker.beta
    bias.corr <- x$parms$Ricker.biascorr
    rec.sim <- alpha * stock.sim * exp(-beta * stock.sim)
    curve.OK <- TRUE
  }
  if (rec.model == "Ricker-steep") {
    R0 <- x$parms$Ricker.R0
    h <- x$parms$Ricker.steep
    phi0 <- x$parms$Ricker.Phi0
    bias.corr <- x$parms$Ricker.biascorr
    rec.sim <- stock.sim/phi0 * exp(h * (1 - stock.sim/(phi0 *
                                                          R0)))
    curve.OK <- TRUE
  }
  if (!curve.OK) {
    warning("Bad value of argument rec.model\n")
    draw.model <- FALSE
    rec.model <- "none"
    bias.corr <- NULL
  }
}
if (draw.model) {
  if (is.numeric(bias.corr) && bias.corr != 1)
    draw.bias <- TRUE
  {
    if (draw.bias)
      lim.y <- range(0, max(bias.corr * rec.sim, ts$recruits[rndx],
                            na.rm = TRUE))
    else lim.y <- range(0, max(rec.sim, ts$recruits[rndx],
                               na.rm = TRUE))
  }
}
if (draft)
  PlotTitle <- FGMakeTitle("Stock-recruitment", DataName)
FGTimePlot(x = ts$SSB[sndx], y = ts$recruits[rndx], y2 = NULL,
           lab.x = lab.x, lab.y = lab.y, FGtype = "circles", main = PlotTitle,
           use.color = use.color, xlim = lim.x, ylim = lim.y)
if (draw.time) {
  lines(ts$SSB[sndx], ts$recruits[rndx], col = clr.points,
        lty = 3)
  if (year.pos > 0) {
    text(x = ts$SSB[sndx[1]], y = ts$recruits[rndx[1]],
         labels = ts$year[rndx[1]], pos = year.pos, cex = 0.85,
         offset = 0.35)
    text(x = ts$SSB[sndx[length(sndx)]], y = ts$recruits[rndx[length(rndx)]],
         labels = ts$year[rndx[length(rndx)]], pos = year.pos,
         cex = 0.85, offset = 0.35)
  }
}
redraw <- FALSE
if (draw.lowess) {
  lines(lowess(ts$SSB[sndx], ts$recruits[rndx], f = 0.55),
        col = clr.lowess, lwd = 2)
  redraw <- TRUE
}
if (draw.model && curve.OK) {
  lines(stock.sim, rec.sim, lwd = 2)
  if (draw.bias) {
    lines(stock.sim, bias.corr * rec.sim, lwd = 2, lty = "dashed")
    leg.text <- c(rec.model, "Expected")
    legend("topright", legend = leg.text, lty = c("solid",
                                                  "dashed"), lwd = 2, inset = 0.01, bg = "white")
    redraw <- TRUE
  }
}
if (redraw) {
  points(ts$SSB[sndx], ts$recruits[rndx], col = clr.points)
}
if (write.graphs)
  FGSavePlot(GraphicsDirName, DataName, GraphName = paste("SR.",
                                                          rec.model, sep = ""), graphics.type)
if (draft)
  PlotTitle <- FGMakeTitle("Stock v log(R/S)", DataName)
logRS = log(ts$recruits[rndx]/ts$SSB[sndx])
lab.y2 <- "log(recruits/spawner)"
if (draw.model && curve.OK) {
  logRS.sim = log(rec.sim/stock.sim)
  lim.y <- range(logRS, logRS.sim, na.rm = TRUE)
}else {
  lim.y <- range(logRS, na.rm = TRUE)
}
if (any(logRS <= 0)) {
  Errmsg <- "Warning: attempted to take log of a non-positive R/S"
  warning(Errmsg, immediate. = TRUE)
  return(invisible(-1))
}
lim.y[1] = 0.9 * lim.y[1]
lim.y[2] = 1 * lim.y[2]
FGTimePlot(x = ts$SSB[sndx], y = logRS, y2 = NULL, lab.x = lab.x,
           lab.y = lab.y2, FGtype = "circles", main = PlotTitle,
           use.color = use.color, xlim = lim.x, ylim = lim.y)
if (draw.time) {
  lines(ts$SSB[sndx], logRS, col = clr.points, lty = 3)
  if (year.pos > 0) {
    text(x = ts$SSB[sndx[1]], y = logRS[1], labels = ts$year[rndx[1]],
         pos = year.pos, cex = 0.85, offset = 0.35)
    text(x = ts$SSB[sndx[length(sndx)]], y = logRS[length(logRS)],
         labels = ts$year[rndx[length(rndx)]], pos = year.pos,
         cex = 0.85, offset = 0.35)
  }
}
redraw <- FALSE
if (draw.lowess) {
  lines(lowess(ts$SSB[sndx], logRS, f = 0.55), col = clr.lowess,
        lwd = 2)
  redraw = TRUE
}
if (draw.model && curve.OK) {
  lines(stock.sim, logRS.sim, lwd = 2)
}
if (redraw) {
  points(ts$SSB[sndx], logRS, col = clr.points)
}
if (write.graphs)
  FGSavePlot(GraphicsDirName, DataName, GraphName = paste("SlogRS.",
                                                          rec.model, sep = ""), graphics.type)
par(savepar)
