# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmIgnition",
  description = "start a random number of fires",
  keywords = c("fire scfm ignition"),
  authors = c(
    person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = numeric_version("1.1.0.9002"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmIgnition.Rmd"), ## same file
  reqdPkgs = list("raster", "SpaDES.tools", "PredictiveEcology/LandR",
                  "PredictiveEcology/scfmutils (>= 1.0.0)"),
  loadOrder = list(after = c("scfmLandcoverInit", "scfmRegime", "scfmDriver"),
                   before = c("scfmEscape", "scfmSpread")),
  parameters = rbind(
    ## TODO: need a Flash parameter controlling fixed number of fires, a la Ratz (1995)
    defineParameter("pIgnition", "numeric", 0.001, 0, 1,
                    "per cell and time ignition probability"),
    defineParameter("startTime", "numeric", start(sim), NA, NA,
                    "simulation time of first ignition"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA,
                    "interval between main events"),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bindrows(
    expectsInput("fireRegimePolys", "sf", "fire regime polys with ignition rate"),
    expectsInput("fireRegimeRas", "SpatRaster", "rasterized version of fire regime polys"),
    expectsInput("flammableMap", "SpatRaster", desc = "map of flammability")
    ),
  outputObjects = bindrows(
    createsOutput("ignitionLoci", "numeric", desc = "vector of ignition locations"),
    createsOutput("pIg", "numeric", desc = "ignition probability raster")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmIgnition = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmIgnition", "ignite", eventPriority = 7.5)
    },
    ignite = {
      sim <- Ignite(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "scfmIgnition", "ignite", eventPriority = 7.5)
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  ## if either of these is a map, it needs to have NAs in the right place
  ##   and be conformant with flammableMap
  if (!is.null(sim$fireRegimePolys$pIgnition)) {

    igValues <- data.table(PolyID = sim$fireRegimePolys$PolyID,
                           pIg = sim$fireRegimePolys$pIgnition)
    igRas <- data.table(PolyID = as.vector(sim$fireRegimeRas),
                        flam = as.vector(sim$flammableMap))
    igValues <- igValues[igRas, on = c("PolyID")]
    igValues[flam != 1, pIg := NA]
    sim$pIg <- rast(sim$fireRegimeRas)
    sim$pIg <- setValues(sim$pIg, igValues$pIg)
  } else {
    warning("using default pIgnition as no 'ignitionRate' column found in fireRegimePolys")
    sim$pIg <- P(sim)$pIgnition #and pIg is a constant from the parameter list
  }

  sim$ignitionLoci <- NULL

  return(invisible(sim))
}

### template for your event1
Ignite <- function(sim) {

  ## TODO: this calcIgnitions could be simpler
  sim$ignitionLoci <- NULL #initialise FFS

  ignitions <- calcIgnitions(fireRegimePolys = sim$fireRegimePolys,
                             pIg = sim$pIg,
                             fireRegimeRas = sim$fireRegimeRas)
  ## resample generates a random permutation of the elements of ignitions
  ## so that we don't always sequence in map index order. EJM pointed this out.
  sim$ignitionLoci <- SpaDES.tools:::resample(unlist(ignitions))

  return(invisible(sim))
}

calcIgnitions <- function(fireRegimePolys, pIg, fireRegimeRas) {

  fireRegimePixels <- 1:ncell(fireRegimeRas)
  fireRegimePixels <- fireRegimePixels[!is.na(as.vector(fireRegimeRas))]

  if (inherits(pIg, "SpatRaster")) {
    stopifnot(ncell(fireRegimeRas) == ncell(pIg))
    igs <- which(runif(ncell(fireRegimeRas)) < as.vector(pIg)) #NA are allowed and omitted

  } else {
    igs <- fireRegimePixels[which(runif(length(fireRegimePixels)) <= pIg)]
  }

  return(igs)
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  mod$dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", mod$dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  if (!suppliedElsewhere("flammableMap", sim)) {
    message("you should run scfmIgnition with scfmLandcoverInit")
    vegMap <- prepInputsLCC(
      year = 2010,
      destinationPath = dPath,
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      userTags = c("prepInputsLCC", "studyArea")
    )
    vegMap[] <- asInteger(vegMap[])
    sim$flammableMap <- defineFlammable(vegMap,
                                        mask = sim$rasterToMatch,
                                        nonFlammClasses = c(13L, 16L:19L)
    )
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
