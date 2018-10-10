
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
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur")
  ),
  inputObjects = data.frame(
    objectName = c("scfmPars","spreadState","flammableMap"),
    objectClass = c("list","data.table","RasterLayer"),
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = "burnMap",
    objectClass = "RasterLayer",
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmSpread = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    sim <- sim$scfmSpreadInit(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, params(sim)$scfmSpread$startTime, "scfmSpread", "burn")
    sim <- scheduleEvent(sim, params(sim)$scfmSpread$.plotInitialTime, "scfmSpread", "plot")
    sim <- scheduleEvent(sim, params(sim)$scfmSpread$.saveInitialTime, "scfmSpread", "save")
  } else if (eventType == "plot") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event
    # sim$burnMap[sim$burnMap[]>0]<-1

    Plot(sim$burnMap, title="Fire map", legendRange=c(0,1), 
         #cols=c("white", "red"), 
         new=FALSE)#length(sim$ignitionLoci)))
    sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmSpread$.plotInterval,
                         "scfmSpread", "plot")
   
    # ! ----- STOP EDITING ----- ! #
  } else if (eventType == "save") {
    
  } else if (eventType == "burn") { 
    if (!is.null(sim$spreadState)){ #we really want to test if the data table has any rows
      if(any(sim$spreadState$active))
        sim<-scfmSpreadBurnemup(sim)
    }
    
    sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmSpread$returnInterval, "scfmSpread", "burn")

  }  else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
scfmSpreadInit <- function(sim) {
  
  sim$burnMap <- sim$flammableMap * 0  # 0 * NA = NA
  
  if("scfmPars" %in% ls(sim)) {
    if(length(sim$landscapeAttr) > 1) {
      pSpread <- raster(sim$flammableMap)
      for(x in names(sim$landscapeAttr)) {
        pSpread[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmPars[[x]]$pSpread
        pSpread[] <- pSpread[] * (1-sim$flammableMap[])
      }
    } else {
      pSpread <- sim$scfmPars[[1]]$pSpread
    }
    
  } else {
    pSpread <- P(sim)$pSpread
  }
  sim$pSpread <- pSpread
  
  
  # pSpread <- ifelse("scfmPars" %in% ls(sim),
  #            sim$scfmPars$pSpread,
  #            params(sim)$scfmSpread$pSpread
  #           )
  
  setColors(sim$burnMap,n=2) <- colorRampPalette(c("transparent", "red"))(2)
  
  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

scfmSpreadBurnemup <- function(sim){ #name is a homage to Walters and Hillborne

  maxSizes <- unlist(lapply(sim$scfmPars, function(x) x$maxBurnCells))
  activeLoci <- unique(sim$spreadState$initialLocus) # indices[sim$spreadState$active]
  maxSizes <- maxSizes[sim$cellsByZone[activeLoci,"zone"]]
  
  sim$burnDT <- SpaDES.tools::spread(sim$flammableMap, 
                 spreadProb=sim$pSpread,
                 spreadState=sim$spreadState,
                 #mask=sim$flammableMap, #this should work but it don't
                 directions=globals(sim)$neighbours,
                 maxSize=maxSizes, returnIndices = TRUE, 
                 id=TRUE)
  sim$burnMap[sim$burnDT$indices] <- 1
  idx<-which(sim$burnMap[] != 0)
  if (length(idx)>0){
      sim$ageMap[idx] <- 0
  }
  return(invisible(sim))
}
### template for save events
scfmSpreadSave <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

### template for plot events
scfmSpreadPlot <- function(sim) {
  return(invisible(sim))
}


### add additional events as needed by copy/pasting from above

