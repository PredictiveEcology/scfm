
# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "ageModule",
  description = "Creates and maintains a raster called ageMap",
  keywords = c("forest age", "modelling course", "Lab 5"),
  authors = c(person(c("Steve", "G"), "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.9.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "ageModule.Rmd"),
  reqdPkgs = list("raster","RColorBrewer"),
  parameters = rbind(
    defineParameter("initialAge", "numeric", 99.0, 0, 1e4, desc =  "initial age"),
    defineParameter("maxAge","numeric", 200, 0, 2**16-1, desc = "maximum age for plotting"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA, desc = "Time interval between aging aevents"),
    defineParameter("startTime", "numeric", 0, NA, NA, desc = "Simulation time at which to initiate aging"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "map of flammability vegetation")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "ageMap", objectClass = "RasterLayer", desc = "map of vegetation age")
  )
))

## event types
#   - type `init` is required for initiliazation

doEvent.ageModule = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    sim <- sim$ageModuleInit(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, params(sim)$ageModule$startTime, "ageModule", "age")
    sim <- scheduleEvent(sim, params(sim)$ageModule$.plotInitialTime, "ageModule", "plot")
    sim <- scheduleEvent(sim, params(sim)$ageModule$.saveInitialTime, "ageModule", "save")
  } else if (eventType=="age") {
    # do stuff for this event
    sim <- ageModuleAge(sim)
    
    # schedule the next event
    sim <- scheduleEvent(sim, time(sim) + params(sim)$ageModule$returnInterval,
                         "ageModule", "age")
  } else if (eventType == "plot") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event
    Plot(sim$ageMap, legendRange=c(0,params(sim)$ageModule$maxAge))
    sim <- scheduleEvent(sim, 
                         time(sim) + params(sim)$ageModule$.plotInterval,
                         "ageModule", "plot")
    # ! ----- STOP EDITING ----- ! #
  }  else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}


### template initilization
ageModuleInit <- function(sim) {

  # # ! ----- EDIT BELOW ----- ! #
  
  if (!("ageMap" %in% ls(sim))){
    N <- sim$mapDim
    x <- raster::extent(c(0,N-1,0,N-1))
    sim$ageMap <- raster(x, nrows = N, ncols = N, vals = P(sim)$ageModule$initialAge) %>%
      setColors(n = 10, colorRampPalette(c("LightGreen", "DarkGreen"))(10))
  }
  else {
    # we will use our colour choices, not whatever may have come with the loaded map.
    setColors(sim$ageMap, n = 10, colorRampPalette(c("LightGreen", "DarkGreen"))(10))
  }
  
   #temporary until we buid the rest of the modules

  return(invisible(sim))
}

### template for save events
ageModuleSave <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

ageModuleAge <- function(sim) {
  
  sim$ageMap <- setValues(sim$ageMap, 
                          pmin(P(sim)$ageModule$maxAge, getValues(sim$ageMap)+
                                             P(sim)$ageModule$returnInterval))
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- inputPath(sim)
  
  if (!suppliedElsewhere(ageMap, sim)) {
    message("age map not supplied. Using default")
    
    ageMapFilename <- file.path(dPath, "age.tif")
    ageMap <- Cache(prepInputs, targetFile = ageMapFilename,
                      studyArea = sim$studyArea,
                      rasterToMatch = sim$vegMap,
                      destinationPath = file.path(dPath, "age"))
    
    sim$ageMap <- ageMap
    
  }
  return(invisible(sim))
}
