# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmIgnition",
  description = "start a random number of fires",
  keywords = c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email="email@example.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.1.0.9002"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "scfmIgnition.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    #need a Flash parameter controlling fixed number of fires, a la Ratz (1995)
    defineParameter("pIgnition", "numeric", 0.001, 0, 1, desc = "per cell and time ignition probability"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA, desc = "interval between main events"),
    defineParameter("startTime", "numeric", 0, NA, NA, desc = "Simulation time at which to initiate ignitions"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, desc = "time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, desc = "time at which the first plot event should occur")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmPars", objectClass = "list", desc = ""),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "map of flammability"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc ="")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "ignitionLoci", objectClass = "numeric", desc = "")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmIgnition = function(sim, eventTime, eventType, debug = FALSE) {
  switch (
    eventType, 
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmIgnition", "ignite")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmIgnition", "plot")
    },
    plot = {
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmIgnition", "plot")
    },
    ignite = {
      sim <- Ignite(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "scfmIgnition", "ignite")
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}



Init <- function(sim) {
  #browser()
  
  #if either of these is a map, it needs to have NAs in the right place
  #and be conformant with flammableMap
  if("scfmPars" %in% ls(sim)) {
    if(length(sim$landscapeAttr) > 1) {
      pIg <- raster(sim$flammableMap)
      for(x in names(sim$landscapeAttr)) {
        pIg[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmPars[[x]]$pIgnition
      }
      pIg[] <- pIg[] * (1-sim$flammableMap[]) #apply non-flammmable 1s and NAs
    } else {
      pIg <- sim$scfmPars[[1]]$pIgnition #and pIg is a constant from scfmDriver
    }
    
  } else {
    pIg <- P(sim)$pIgnition #and pIg is a constant from the parameter list
  }
  sim$pIg <- pIg
  
  sim$ignitionLoci <- NULL
  
  return(invisible(sim))
}


### template for your event1
Ignite <- function(sim) {
  sim$ignitionLoci <- NULL #initialise FFS
  ignitions <- lapply(names(sim$landscapeAttr), function(polygonType) {
    cells <- sim$landscapeAttr[[polygonType]]$cellsByZone
    if(length(sim$pIg)>1) {
      cells[which(runif(length(cells)) < sim$pIg[cells])]
    } else {
      cells[which(runif(length(cells)) < sim$pIg)]
    }
  })
  #resample generates a random permutation of the elements of ignitions
  #so that we don't always sequence in map index order. EJM pointed this out.
  sim$ignitionLoci <- SpaDES.tools:::resample(unlist(ignitions)) 
  
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
  if (!suppliedElsewhere(flammableMap, sim)){
    sim$flammableMap <- sim$ageMap
    sim$flammableMap[] <- sim$ageMap[] * 0
  }
  
  return(invisible(sim))
}

