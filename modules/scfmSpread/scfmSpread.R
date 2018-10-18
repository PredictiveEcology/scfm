
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmSpread",
  description = "insert module description here",
  keywords = c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email="email@example.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.1.0.9002"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "scfmSpread.Rmd"),
  reqdPkgs = list("raster","data.table","magrittr"),
  parameters = rbind(
    defineParameter("pSpread","numeric", 0.23, 0, 1, desc="spread module for BEACONs"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA, desc="Time interval between burn events"),
    defineParameter("startTime", "numeric", 0, NA, NA, desc="Simulation time at which to initiate burning"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    #defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    #defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of immediate cell neighbours")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmPars", objectClass = "list", desc = ""),
    expectsInput(objectName = "spreadState", objectClass = "data.table", desc = ""),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "burnMap", objectClass = "RasterLayer", desc= "")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <-Init(sim)
      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmSpread", "burn")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmSpread", "plot")
    },
    plot = {
      browser()
      Plot(sim$burnMap, title="Fire map", legendRange=c(0,1), new=FALSE)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmSpread", "plot")
    },
    burn = { 
      if (!is.null(sim$spreadState)){  #we really want to test if the data table has any rows
        if(any(sim$spreadState$active))
          sim<-Burnemup(sim)
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
  
  if("scfmPars" %in% ls(sim)) {
    if(length(sim$landscapeAttr) > 1) {
      pSpread <- raster(sim$flammableMap)
      for(x in names(sim$landscapeAttr)) {
        pSpread[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmPars[[x]]$pSpread
      }     
      pSpread[] <- pSpread[] * (1-sim$flammableMap[])
    } else {
      pSpread <- sim$scfmPars[[1]]$pSpread
    }
    
  } else {
    pSpread <- P(sim)$pSpread
  }
  sim$pSpread <- pSpread
  setColors(sim$burnMap,n=2) <- colorRampPalette(c("transparent", "red"))(2)
  return(invisible(sim))
}

Burnemup <- function(sim){ #name is a homage to Walters and Hillborne

  maxSizes <- unlist(lapply(sim$scfmPars, function(x) x$maxBurnCells))
  activeLoci <- unique(sim$spreadState$initialLocus) # indices[sim$spreadState$active]
  #we prevent multiple ignitions, which shouldn't happen anyway.
  maxSizes <- maxSizes[sim$cellsByZone[activeLoci,"zone"]]
  
  sim$burnDT <- SpaDES.tools::spread(sim$flammableMap,
                                     spreadProb = sim$pSpread,
                                     spreadState = sim$spreadState,
                                     directions = P(sim)$neighbours,
                                     maxSize = maxSizes,  #not sure this works 
                                     returnIndices = TRUE, 
                                     id = TRUE)
  sim$burnMap[sim$burnDT$indices] <- 1
  sim$ageMap[sim$burnDT$indices] <- 0
  return(invisible(sim))
}

