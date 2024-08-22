# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmEscape",
  description = paste("'Escapes' fire(s) from an initial set of loci returned by `scfmIgnition`,",
                      "and prepares the results for use by `scfmSpread`."),
  keywords = c("fire escape", "BEACONs"),
  authors = c(
    person(c("Steven", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
    person("Ian MS", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = list(scfmEscape = "2.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmEscape.Rmd"),
  reqdPkgs = list("data.table",
                  "PredictiveEcology/LandR (>= 1.1.1)",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/scfmutils@development (>= 2.0.1)",
                  "sf",
                  "PredictiveEcology/SpaDES.tools@development",
                  "terra"),
  parameters = rbind(
    defineParameter("dataYear", "numeric", 2011, 1985, 2020,
                    paste("used to select the year of landcover data used to create",
                          "`flammableMap` if the obejct is unsupplied.")),
    defineParameter("neighbours", "integer", 8L, 4L, 8L,
                    "Number of cell immediate neighbours (one of `4L` or `8L`)."),
    defineParameter("p0", "numeric", 0.1, 0, 1,
                    "probability of an ignition spreading to an unburned immediate neighbour"),
    defineParameter("returnInterval", "numeric", 1, NA, NA,
                    "This specifies the time interval between Escape events"),
    defineParameter("startTime", "numeric", start(sim, "year"), NA, NA,
                    "simulation time of first escape"),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    "Internal. Can be names of events or the whole module name; these will be cached by SpaDES.")
  ),
  inputObjects = bindrows(
    expectsInput("fireRegimePolys", "sf",
                 desc = "fire regime polys with ignition rate"),
    expectsInput("fireRegimeRas", "SpatRaster",
                 desc = "rasterized version of `fireRegimePolys`"),
    expectsInput("flammableMap", "SpatRaster",
                 desc = "map of flammability"),
    expectsInput("ignitionLoci", "numeric",
                 desc = "pixel IDs where ignition occurs")
  ),
  outputObjects = bindrows(
    createsOutput("spreadState", "data.table", desc = "stores the current fire spread state"),
    createsOutput("p0", "SpatRaster", desc = "escape probability raster")
  )
))

## event types
#   - type `init` is required for initiliazation

doEvent.scfmEscape = function(sim, eventTime, eventType, debug = FALSE){
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmEscape", "escape", eventPriority = 7.5)
    },
    escape = {
      sim <- Escape(sim)

      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "scfmEscape", "escape", eventPriority = 7.5)
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### template initilization
Init <- function(sim) {
  sim$spreadState <- NULL

  if (!is.null(sim$fireRegimePolys$p0)) {
    escValues <- data.table(PolyID = sim$fireRegimePolys$PolyID,
                            p0 = sim$fireRegimePolys$p0)
    escRas <- data.table(PolyID = as.vector(sim$fireRegimeRas),
                         flam = as.vector(sim$flammableMap))
    escValues <- escValues[escRas, on = c("PolyID")]
    escValues[flam != 1, p0 := NA]
    sim$p0 <- rast(sim$fireRegimeRas)
    sim$p0 <- setValues(sim$p0, escValues$p0)
  } else {
    warning("using default escape prob as no `xxx` column found in fireRegimePolys")
    sim$p0 <- P(sim)$p0
  }

  return(invisible(sim))
}

Escape <- function(sim) {

  if (length(sim$ignitionLoci) > 0) {
    sim$spreadState <- spread2(landscape = sim$flammableMap,
                               start = sim$ignitionLoci,
                               iterations = 1,
                               spreadProb = sim$p0,
                               directions = P(sim)$neighbours,
                               asRaster = FALSE)
  }

  return(invisible(sim))
}

## same model as scfmIgnition to enable standalone execution
.inputObjects <- function(sim) {
  ## TODO: This module has other dependencies that aren't created: scfmDriverPars and ignitionLoci
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)

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
  return(invisible(sim))
}
