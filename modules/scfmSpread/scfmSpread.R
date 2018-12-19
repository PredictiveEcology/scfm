# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmSpread",
  description = "insert module description here",
  keywords = c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email = "email@example.com", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.1.0.9002"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "scfmSpread.Rmd"),
  reqdPkgs = list("raster","data.table","magrittr"),
  parameters = rbind(
    defineParameter("pSpread", "numeric", 0.23, 0, 1, desc = "spread module for BEACONs"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA, desc = "Time interval between burn events"),
    defineParameter("startTime", "numeric", 0, NA, NA, desc = "Simulation time at which to initiate burning"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    #defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    #defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of immediate cell neighbours")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmDriverPars", objectClass = "list", desc = "fire modules' parameters"),
    expectsInput(objectName = "spreadState", objectClass = "data.table", desc = ""),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "burnMap", objectClass = "RasterLayer", desc = "cumulative burn map"),
    createsOutput(objectName = "burnDT", objectClass = "data.table", desc = "data table with pixel IDs of most recent burn"),
    createsOutput(objectName = "rstCurrentBurn", object = "RasterLayer", desc = "annual burn map"),
    createsOutput(objectName = "pSpread", object = "RasterLayer", desc = "spread probability applied to flammabiliy Map")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmSpread", "burn")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmSpread", "plot")
    },
    plot = {
      Plot(sim$burnMap, legend = FALSE)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmSpread", "plot")
    },
    burn = {
      if (!is.null(sim$spreadState)) {
        ## we really want to test if the data table has any rows
        if (NROW(sim$spreadState[state == "activeSource"]) > 0)
          sim <- Burnemup(sim)
      }
      sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmSpread$returnInterval, "scfmSpread", "burn")
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {

  sim$burnMap <- sim$flammableMap
  sim$burnMap[] <- sim$flammableMap[] * 0  # 0 * NA = NA

  if ("scfmDriverPars" %in% ls(sim)) {
    if (length(sim$scfmDriverPars) > 1) {
      pSpread <- raster(sim$flammableMap)
      for (x in names(sim$scfmDriverPars)) {
        pSpread[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmDriverPars[[x]]$pSpread
      }
      pSpread[] <- pSpread[] * (sim$flammableMap[])
    } else {
      pSpread <- sim$flammableMap * sim$scfmDriverPars[[1]]$pSpread
    }
  } else {
    pSpread <- sim$pSpread * sim$flammableMap
  }
  sim$pSpread <- pSpread
  setColors(sim$burnMap, n = 2) <- colorRampPalette(c("transparent", "red"))(2)
  return(invisible(sim))
}

Burnemup <- function(sim){ #name is a homage to Walters and Hillborne

  # maxSizes <- unlist(lapply(sim$scfmDriverPars, function(x) x$maxBurnCells))
  # activeLoci <- unique(sim$spreadState$initialLocus) # indices[sim$spreadState$active]
  #we prevent multiple ignitions, which shouldn't happen anyway.
  # maxSizes <- maxSizes[sim$cellsByZone[activeLoci, "zone"]]

  sim$burnDT <- SpaDES.tools::spread2(sim$flammableMap,
                                      start = sim$spreadState,
                                     spreadProb = sim$pSpread,
                                     #spreadState = sim$spreadState,
                                     directions = P(sim)$neighbours,
                                     # maxSize = maxSizes,  #not sure this works
                                     asRaster = FALSE)

  sim$rstCurrentBurn <- sim$vegMap #This preserves NAs
  sim$rstCurrentBurn[!is.na(sim$rstCurrentBurn)] <- 0 #reset annual burn
  sim$rstCurrentBurn[sim$burnDT$pixels] <- 1 #update annual burn
  sim$burnMap[sim$burnDT$pixels] <- 1 #update cumulative burn
  sim$ageMap[sim$burnDT$pixels] <- 0 #update age
  return(invisible(sim))
}
