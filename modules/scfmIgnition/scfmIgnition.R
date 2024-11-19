defineModule(sim, list(
  name = "scfmIgnition",
  description = "start a random number of fires",
  keywords = c("fire scfm ignition"),
  authors = c(
    person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = list(scfmIgnition = "2.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmIgnition.Rmd"), ## same file
  reqdPkgs = list(
    "sf", "SpaDES.tools", "terra",
    "PredictiveEcology/LandR (>= 1.1.1)",
    "PredictiveEcology/scfmutils@development (>= 2.0.1)"
  ),
  loadOrder = list(after = c("scfmLandcoverInit", "scfmRegime", "scfmDriver"),
                   before = c("scfmEscape", "scfmSpread")),
  parameters = rbind(
    ## TODO: need a Flash parameter controlling fixed number of fires, a la Ratz (1995)
    defineParameter("dataYear", "numeric", 2011, 1985, 2020,
                    paste("used to select the year of landcover data used to create",
                          "`flammableMap` if the object is unsupplied")),
    defineParameter("pIgnition", "numeric", 0.001, 0, 1,
                    "default per-cell and time ignition probability if unsupplied."),
    defineParameter("startTime", "numeric", start(sim), NA, NA,
                    "simulation time of first ignition"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA,
                    "interval between main events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bindrows(
    expectsInput("fireRegimePolys", "sf", "`fireRegimePolys` with ignition rate attribute"),
    expectsInput("fireRegimeRas", "SpatRaster", "rasterized version of `fireRegimePolys`"),
    expectsInput("flammableMap", "SpatRaster", desc = "map of flammability")
    ),
  outputObjects = bindrows(
    createsOutput("ignitionLoci", "numeric", desc = "vector of ignition locations"),
    createsOutput("pIg", "SpatRaster", desc = "ignition probability raster")
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
    vegMap <- prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      maskTo = sim$studyArea,
      cropTo = sim$rasterToMatch,
      projectTo = sim$rasterToMatch,
      userTags = c("prepInputs_NTEMS_LCC_FAO", "studyArea")
    )
    vegMap[] <- asInteger(vegMap[])
    sim$flammableMap <- defineFlammable(vegMap,
                                        mask = sim$rasterToMatch,
                                        nonFlammClasses = c(20, 31, 32, 33)
    )
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
