# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmEscape",
  description = "This Escapes fire(s) from an initial set of loci returned by an ignition module, and readies the results for use by scfmSpread",
  keywords = c("fire Escape BEACONs"),
  authors = c(person(c("Steven", "G"), "Cumming",
                     email = "stevec@sbf.ulaval.ca", role = c("aut"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "scfmEscape.Rmd"),
  reqdPkgs = list("raster","data.table","magrittr"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("p0", "numeric", 0.1, 0, 1, "probability of an ignition spreading to an unburned immediate neighbour"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "time at which the first plot event should occur"),
    #defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "time at which the first save event should occur"),
    #defineParameter(".saveInterval", "numeric", NA, NA, NA, "time at which the first save event should occur"),
    defineParameter("returnInterval", "numeric", NA, NA, NA, "This specifies the time interval between Escape events"),
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of cell immediate neighbours")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmDriverPars", objectClass = "list", desc = "fire modules' parameters"),
    expectsInput(objectName = "ignitionLoci", objectClass = "numeric", desc = ""),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = ""),
    expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer", desc = "template raster for raster GIS operations. Must be supplied by user")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "spreadState", objectClass = "data.table", desc = "")
  )
))

## event types
#   - type `init` is required for initiliazation

doEvent.scfmEscape = function(sim, eventTime, eventType, debug = FALSE){
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmEscape", "escape")
      sim <- scheduleEvent(sim, params(sim)$scfmEscape$.plotInitialTime, "scfmEscape", "plot")
      #sim <- scheduleEvent(sim, params(sim)$scfmEscape$.saveInitialTime, "scfmEscape", "save")
    },
    plot = {
      values(sim$tmpRaster)[sim$spreadState[, indices]] <- 2 #this reference method is believed to be faster
      values(sim$tmpRaster)[sim$ignitionLoci] <- 1          #mark the initials specially
      Plot(sim$tmpRaster)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmEscape", "plot")
    },
    escape = {
      sim <- Escape(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "scfmEscape", "escape")
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initilization
Init <- function(sim) {
  sim$spreadState <- NULL

  if ("scfmDriverPars" %in% ls(sim)) {
    if (length(sim$scfmDriverPars) > 1) {
      p0 <- raster(sim$flammableMap)
      for (x in names(sim$scfmDriverPars)) {
        p0[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmDriverPars[[x]]$p0
      }
      p0[] <- p0[] * (1 - sim$flammableMap[])
    } else {
      p0 <- sim$scfmDriverPars[[1]]$p0
    }
  } else {
    p0 <- P(sim)$p0
  }
  sim$p0 <- p0

  return(invisible(sim))
}

Escape <- function(sim) {
  if (length(sim$ignitionLoci) > 0) {
    # print(paste("Year",time(sim), "loci = ", length(sim$ignitionLoci)))
    maxSizes <- unlist(lapply(sim$scfmDriverPars, function(x) x$maxBurnCells))
    maxSizes <- maxSizes[sim$cellsByZone[sim$ignitionLoci, "zone"]]

    sim$spreadState <- SpaDES.tools::spread2(landscape = sim$flammableMap,
                                            start = sim$ignitionLoci,
                                            iterations = 1,
                                            spreadProb = sim$p0,
                                            directions = P(sim)$neighbours,
                                            asRaster = FALSE,
                                            maxSize = maxSizes)


  }

  return(invisible(sim))
}

#same model as scfmIgnition to enable standalone execution
.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  #This module has many dependencies that aren't sourced in .inputObjects.
  if (!suppliedElsewhere("flammableMap", sim)) {
    sim$flammableMap <- sim$rasterToMatch
    sim$flammableMap[] <- sim$ageMap[] * 0
  }
  return(invisible(sim))
}
