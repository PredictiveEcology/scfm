defineModule(sim, list(
  name = "ageModule",
  description = "Creates and maintains a raster called ageMap",
  keywords = c("forest age", "modelling course", "Lab 5"),
  authors = c(
    person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = list(LandR = "0.0.11.9000", ageModule = "0.9.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "ageModule.Rmd"),
  reqdPkgs = list("RColorBrewer", "sf", "terra",
                  "PredictiveEcology/LandR@development",
                  "PredictiveEcology/scfmutils@development (>= 0.0.13.9003)"),
  parameters = rbind(
    defineParameter("initialAge", "numeric", 99.0, 0, 1e4, desc =  "initial age"),
    defineParameter("maxAge", "numeric", 200, 0, 2**16 - 1, desc = "maximum age for plotting"),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA, desc = "Time interval between aging events"),
    defineParameter("startTime", "numeric", start(sim), NA, NA, desc = "Simulation time at which to initiate aging"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by Plots function, which can be optionally used here")
  ),
  inputObjects = bindrows(
    expectsInput("ageMap", "SpatRaster",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureStandVolume.tar"),
    expectsInput("studyArea", "sf",
                 desc = "study area template",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = "template raster for raster GIS operations. Must be supplied by user."),
    expectsInput("rstCurrentBurn", "SpatRaster",
                 desc = "annual burn map created by `scfmSpread`.")
  ),
  outputObjects = bindrows(
    createsOutput("ageMap", "SpatRaster", desc = "map of vegetation age")
  )
))

doEvent.ageModule = function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$startTime, "ageModule", "age", eventPriority = 7.5)

      if (anyPlotting(P(sim)$.plots)) {
        sim <- scheduleEvent(sim, P(sim)$startTime, "ageModule", "plot", eventPriority = 7.5)
      }
    },
    age = {
      # do stuff for this event
      sim <- Age(sim)

      # schedule the next event
      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "ageModule", "age")
    },
    plot = {
      Plots(sim$ageMap, fn = scfmutils::plot_ageMap, type = P(sim)$.plots,
            filename = paste0("ageMap_year_", time(sim)),
            title = paste0("Age map: year ", time(sim)),
            maxAge = P(sim)$maxAge)

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "ageModule", "plot")
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
               "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  ## TODO: remove this workaround -- why isn't this 'being 'sticking' when done in .inputObjects??
  if (!compareRaster(sim$rasterToMatch, sim$ageMap, stopiffalse = FALSE,
                     extent = TRUE, rowcol = TRUE, crs = TRUE, res = TRUE)) {
    ## ensure ageMap matches rasterToMatch
    useTerra <- getOption("reproducible.useTerra") ## TODO: reproducible#242
    options(reproducible.useTerra = FALSE) ## TODO: reproducible#242
    sim$ageMap <- postProcess(sim$ageMap, rasterToMatch = sim$rasterToMatch)
    options(reproducible.useTerra = useTerra) ## TODO: reproducible#242
  }

  ## we will use our colour choices, not whatever may have come with the loaded map.
  setColors(sim$ageMap, n = 10, colorRampPalette(c("LightGreen", "DarkGreen"))(10))

  return(invisible(sim))
}

Age <- function(sim) {
  newAges <- pmin(P(sim)$maxAge, getValues(sim$ageMap) + P(sim)$returnInterval)
  sim$ageMap <- setValues(sim$ageMap, newAges)

  if (!is.null(sim$rstCurrentBurn)) {
    compareRaster(sim$rasterToMatch, sim$ageMap, sim$rstCurrentBurn,
                  extent = TRUE, rowcol = TRUE, crs = TRUE, res = TRUE)
    burn <- getValues(sim$rstCurrentBurn)
    sim$ageMap[!is.na(burn) & burn == 1] <- 0
  }

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(inputPath(sim), 1)

  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using random polygon in Alberta")

    studyArea <- pemisc::randomStudyArea(size = 2000000000, seed = 23657)
    sim$studyArea <- studyArea
  }

  if (!suppliedElsewhere("ageMap", sim)) {
    sim$ageMap <- LandR::prepInputsStandAgeMap(
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      destinationPath = dPath,
      startTime  = start(sim)
    )
  }

  if (!compareRaster(sim$rasterToMatch, sim$ageMap, stopiffalse = FALSE,
                     extent = TRUE, rowcol = TRUE, crs = TRUE, res = TRUE)) {
    ## ensure ageMap matches rasterToMatch
    useTerra <- getOption("reproducible.useTerra") ## TODO: reproducible#242
    options(reproducible.useTerra = FALSE) ## TODO: reproducible#242
    sim$ageMap <- postProcess(sim$ageMap, rasterToMatch = sim$rasterToMatch)
    options(reproducible.useTerra = useTerra) ## TODO: reproducible#242
  }

  if (!suppliedElsewhere("rstCurrentBurn", sim)) {
    stop("Please supply rstCurrentBurn by running scfmSpread or another fire model")
  }

  return(invisible(sim))
}
