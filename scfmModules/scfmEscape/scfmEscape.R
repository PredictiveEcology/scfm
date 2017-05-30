# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmEscape",
  description = "This Escapes fire(s) from an initial set of loci returned by an ignition module,\ 
                 and readies the results for use by scfmSpread",
  keywords = c("fire Escape BEACONs"),
  authors = c(person(c("Steve", "X"), "Cumming", 
                     email="stevec@sbf.ulaval.ca", role=c("aut"))),
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
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter("returnInterval", "numeric", NA, NA, NA, "This specifies the time interval between Escape events")
  ),
  inputObjects = data.frame(
    objectName = c("scfmPars","ignitionLoci","flammableMap"),
    objectClass = c("list","numeric", "RasterLayer"),
    sourceURL = c("","",""),
    other = c(NA_character_,NA_character_,NA_character_), 
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = c("spreadState"),
    objectClass =c("data.table"),
    other = c(NA_character_),
    stringsAsFactors = FALSE
  )
))


## event types
#   - type `init` is required for initiliazation

doEvent.scfmEscape = function(sim, eventTime, eventType, debug = FALSE){
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)
    
    # do stuff for this event
    sim <- sim$scfmEscapeInit(sim)
    
    # schedule future event(s)
    sim <- scheduleEvent(sim, params(sim)$scfmEscape$startTime, "scfmEscape", "escape")
    #sim <- scheduleEvent(sim, params(sim)$scfmEscape$.plotInitialTime, "scfmEscape", "plot")
    #sim <- scheduleEvent(sim, params(sim)$scfmEscape$.saveInitialTime, "scfmEscape", "save")
  } else if (eventType == "plot") {
    
    values(sim$tmpRaster)[sim$spreadState[,indices]]<-2 #this reference method is believed to be faster
    values(sim$tmpRaster)[sim$loci]<-1    #mark the initials specialy
    Plot(sim$tmpRaster) 
    sim <- scheduleEvent(sim, time(sim)+params(sim)$scfmEscape$.plotInterval, "scfmEscape", "plot")
    
  } else if (eventType == "escape") {
    sim <- sim$scfmEscapeEscape(sim)
    sim <- scheduleEvent(sim, time(sim)+params(sim)$scfmEscape$returnInterval, "scfmEscape", "escape")
    
    
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initilization
scfmEscapeInit <- function(sim) {
  
  sim$spreadState <- NULL
  
  if("scfmPars" %in% ls(sim)) {
    if(length(sim$landscapeAttr) > 1) {
      p0 <- raster(sim$flammableMap)
      for(x in names(sim$landscapeAttr)) {
        p0[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmPars[[x]]$p0
        p0[] <- p0[] * (1-sim$flammableMap[])
      }
    } else {
      p0 <- sim$scfmPars[[1]]$p0
    }
    
  } else {
    p0 <- P(sim)$p0
  }
  sim$p0 <- p0
  
  
  return(invisible(sim))
}


scfmEscapeEscape <- function(sim) {
  if (length(sim$ignitionLoci) > 0){
    print(paste("Year",time(sim), "loci = ", length(sim$ignitionLoci)))
    # pull 
    maxSizes <- unlist(lapply(sim$scfmPars, function(x) x$maxBurnCells))
    maxSizes <- maxSizes[sim$cellsByZone[sim$ignitionLoci,"zone"]]
    
    sim$spreadState <- SpaDES::spread(landscape=sim$flammableMap,
                                      loci=sim$ignitionLoci,
                                      iterations=1,
                                      spreadProb=sim$p0,
                                      #mask=sim$flammableMap,
                                      directions=globals(sim)$neighbours,
                                      maxSize=maxSizes,
                                      returnIndices=TRUE, id=TRUE)
  } 
  return(invisible(sim))
}
