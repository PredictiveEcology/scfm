# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfmEscape",
  description = paste("'Escapes' fire(s) from an initial set of loci returned by an ignition module,",
                      "and prepares the results for use by `scfmSpread`."),
  keywords = c("fire escape", "BEACONs"),
  authors = c(
    person(c("Steven", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
    person("Ian MS", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = numeric_version("0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmEscape.Rmd"),
  reqdPkgs = list("data.table", "raster",
                  "PredictiveEcology/LandR (>= 1.1.1)",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/SpaDES.tools@development"),
  parameters = rbind(
    defineParameter("neighbours", "integer", 8L, 4L, 8L, "Number of cell immediate neighbours (one of `4L` or `8L`)."),
    defineParameter("p0", "numeric", 0.1, 0, 1, "probability of an ignition spreading to an unburned immediate neighbour"),
    defineParameter("returnInterval", "numeric", 1, NA, NA, "This specifies the time interval between Escape events"),
    defineParameter("startTime", "numeric", start(sim, "year"), NA, NA, "simulation time of first escape"),
    defineParameter(".plotInitialTime", "numeric", start(sim, "year") + 1, NA, NA, "time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "time at which the first plot event should occur"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name; these will be cached by SpaDES.")
  ),
  inputObjects = bindrows(
    expectsInput("flammableMap", "RasterLayer", desc = "binary map of landscape flammability"),
    expectsInput("ignitionLoci", "numeric", desc = "Pixel IDs where ignition occurs"),
    expectsInput("rasterToMatch", "RasterLayer", desc = "template raster for raster GIS operations. Must be supplied by user"),
    expectsInput("scfmDriverPars", "list", desc = "fire modules' parameters")
  ),
  outputObjects = bindrows(
    createsOutput("spreadState", "data.table", desc = ""),
    createsOutput("p0", "raster", desc = "")
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

      if ("screen" %in% P(sim)$.plots) {
        sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmEscape", "plot", eventPriority = 7.5)
      }
    },
    plot = {
      if (time(sim) > start(sim)) {
        escapeRaster <- rast(sim$rasterToMatch)
        values(escapeRaster)[sim$spreadState[, pixels]] <- 2 # this reference method is believed to be faster
        values(escapeRaster)[sim$ignitionLoci] <- 1          # mark the initials specially
        Plot(escapeRaster, title = paste0("Annual fire escapes", time(sim)))
      }
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmEscape", "plot", eventPriority = 7.5)
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

  if ("scfmDriverPars" %in% ls(sim)) {
    if (length(sim$scfmDriverPars) > 1) {
      p0 <- rast(sim$flammableMap)
      for (x in names(sim$scfmDriverPars)) {
        values(p0)[sim$landscapeAttr[[x]]$cellsByZone] <- sim$scfmDriverPars[[x]]$p0
      }
      p0 <- p0 * sim$flammableMap
    } else {
      p0 <- sim$scfmDriverPars[[1]]$p0
    }
  } else {
    p0 <- P(sim)$p0
  }
  sim$p0 <- p0

  return(invisible(sim))
}

Escape <- function(sim) {
  if (length(sim$ignitionLoci) > 0) {
    # print(paste("Year",time(sim), "loci = ", length(sim$ignitionLoci)))

    maxSizes <- unlist(lapply(sim$scfmDriverPars, function(x) x$maxBurnCells))
    cellsByZone <- as.character(sim$fireRegimeRas[sim$ignitionLoci])
    maxSizes <- maxSizes[cellsByZone]
    sim$spreadState <- SpaDES.tools::spread2(landscape = sim$flammableMap,
                                             start = sim$ignitionLoci,
                                             iterations = 1,
                                             spreadProb = sim$p0,
                                             directions = P(sim)$neighbours,
                                             asRaster = FALSE,
                                             maxSize = maxSizes)
  }

  return(invisible(sim))
}

## same model as scfmIgnition to enable standalone execution
.inputObjects <- function(sim) {
  ## TODO: This module has other dependencies that aren't created: scfmDriverPars and ignitionLoci

  if (!suppliedElsewhere("flammableMap", sim)) {
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
  return(invisible(sim))
}
