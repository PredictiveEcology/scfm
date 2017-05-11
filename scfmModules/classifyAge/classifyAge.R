# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "classifyAge",
  description = "produce a binary map of old forest cells",
  keywords = c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email="email@example.com", role=c("aut", "cre"))),
  childModules = character(0),
  version = numeric_version("1.3.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "classifyAge.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("returnInterval", "numeric", 10, 1, NA, "Probably should not be an annual event"),
    defineParameter("startTime", "numeric", 1, 0, NA, "Probably should not be an annual event"),
    defineParameter("howOldIsOld", "numeric", 100, 1, NA, "Age of old forest is context dependent"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "ageMap", objectClass = "RasterLayer", desc = "the age Map to classify", sourceURL = NA)
   #expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "ids cells that can't burn", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "ageClassMap", objectClass = "RasterLayer", desc = "the classified binary map")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.classifyAge = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    sim <- sim$classifyAgeInit(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, P(sim)$startTime, "classifyAge", "classify")
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "classifyAge", "plot")
  } else if (eventType == "plot") {
    Plot(sim$ageClassMap) # uncomment this, replace with object to plot
    sim <- scheduleEvent(sim, time(sim)+P(sim)$.plotInterval, "classifyAge", "plot")
  } else if (eventType == "classify") {
    sim$ageClassMap[] <-ifelse(sim$ageMap[] > P(sim)$howOldIsOld, 1, 0)
    
    sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "classifyAge", "classify")
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}


### template initialization
classifyAgeInit <- function(sim) {
  if (!("ageClassMap" %in% names(objs(sim)))){
      if ("ageMap" %in% names(objs(sim))){
        sim$ageClassMap<-sim$ageMap
        sim$ageClassMap[]<-NA
      }
      else {
        stop("need to give me something!")
      }
    }
 return(invisible(sim))
}

