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
    defineParameter("pIgnition", "numeric", 0.001, 0, 1, desc="per cell and time ignition probability"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA, desc="interval between main events"),
    defineParameter("startTime", "numeric", 0, NA, NA, desc="Simulation time at which to initiate ignitions"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, desc="This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, desc="This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, desc="This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, desc="This describes the simulation time at which the first save event should occur")
  ),
  inputObjects = data.frame(
    objectName = c("scfmPars","flammableMap","landscapeAttr"),
    objectClass = c("list","RasterLayer","list"),
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = "ignitionLoci",
    objectClass = "numeric",
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmIgnition = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    sim <- sim$scfmIgnitionInit(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, params(sim)$scfmIgnition$startTime, "scfmIgnition", "ignite")
    sim <- scheduleEvent(sim, params(sim)$scfmIgnition$.plotInitialTime, "scfmIgnition", "plot")
    sim <- scheduleEvent(sim, params(sim)$scfmIgnition$.saveInitialTime, "scfmIgnition", "save")
  } else if (eventType == "plot") {
    sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmIgnition$.plotInterval, "scfmIgnition", "plot")
  } else if (eventType == "save") {
  } else if (eventType == "ignite") {
    sim <- scfmIgnitionIgnite(sim)
    sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmIgnition$returnInterval, "scfmIgnition", "ignite")
    # ! ----- STOP EDITING ----- ! #
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


scfmIgnitionInit <- function(sim) {
  #browser()
  if (!("flammableMap" %in% ls(sim))){
    if ("ageMap" %in% ls(sim)){
      sim$flammableMap<-sim$ageMap
      sim$flammableMap[]<-0
    }
    else {
      stop("need to give me something!")
    }

  }
  
  #if either of these is a map, it needs to have NAs in the right place
  #and be conformant with flammableMap
  if("scfmPars" %in% ls(sim)) {
    if(length(sim$landscapeAttr) > 1) {
      pIg <- raster(sim$flammableMap)
      for(x in names(sim$landscapeAttr)) {
        pIg[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmPars[[x]]$pIgnition
        pIg[] <- pIg[] * (1-sim$flammableMap[])
      }
    } else {
      pIg <- sim$scfmPars[[1]]$pIgnition
    }
    
  } else {
    pIg <- params(sim)$scfmIgnition$pIgnition
  }
  sim$pIg <- pIg
  
  sim$ignitionLoci <- NULL
  
  return(invisible(sim))
}

### template for save events
scfmIgnitionSave <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

### template for your event1
scfmIgnitionIgnite <- function(sim) {
  
  ignitions <- lapply(names(sim$landscapeAttr), function(polygonType) {
    cells <- sim$landscapeAttr[[polygonType]]$cellsByZone
    if(length(sim$pIg)>1) {
      cells[which(runif(length(cells)) < sim$pIg[cells])]
    } else {
      cells[which(runif(length(cells)) < sim$pIg)]
    }
  })
  #randomise sequence to remove any effects of correlation between location in map and sequence of spread. 
  sim$ignitionLoci <- SpaDES:::resample(unlist(ignitions)) # beware case where only 1 fire... see ?sample example
  
  return(invisible(sim))
}


