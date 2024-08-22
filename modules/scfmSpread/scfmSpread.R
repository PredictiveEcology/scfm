defineModule(sim, list(
  name = "scfmSpread",
  description = "model fire spread",
  keywords = c("fire", "spread", "scfm"),
  authors = c(
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = list(scfmSpread = "2.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "scfmSpread.Rmd"),
  loadOrder = list(after = c("scfmLandcoverInit", "scfmRegime", "scfmDriver", "scfmIgnition", "scfmEscape"),
                   before = c("Biomass_regeneration", "Biomass_regenerationPM")),
  reqdPkgs = list("data.table", "fpCompare", "sf", "terra", "viridis",
                  "PredictiveEcology/LandR (>= 1.1.1)",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/SpaDES.tools (>= 2.0.7)",
                  "PredictiveEcology/scfmutils@development (>= 2.0.1)"),
  parameters = rbind(
    defineParameter("dataYear", "numeric", 2011, 1985, 2020,
                    desc = paste("used to select the year of landcover data used to create",
                                 "`flammableMap` if the obejct is unsupplied")),
    defineParameter("neighbours", "numeric", 8, NA, NA,
                    desc = "Number of immediate cell neighbours"),
    defineParameter("pSpread", "numeric", 0.23, 0, 1,
                    desc = "default spread probability if `fireRegimePolys` is missing attribute `pSpread`."),
    defineParameter("returnInterval", "numeric", 1.0, NA, NA,
                    desc = "Time interval between burn events"),
    defineParameter("startTime", "numeric", start(sim), NA, NA,
                    desc = "Simulation time at which to initiate burning"),
    defineParameter(".plotInterval", "numeric", 10, NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plots", "character", c("screen"), NA, NA,
                    desc = "Used by Plots function, which can be optionally used here"),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    desc = "Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bindrows(
    expectsInput("fireRegimePolys", "sf",
                 desc = "`fireRegimePolys` with fire attributes appended."),
    expectsInput("fireRegimeRas", "SpatRaster",
                 desc = "raster with fire regimes from `fireRegimePolys`."),
    expectsInput("flammableMap", "SpatRaster",
                 desc = "binary map of landscape flammability"),
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = "template raster for raster GIS operations. Must be supplied by user."),
    expectsInput("spreadState", "data.table", desc = "see `SpaDES.tools::spread2`"),
    expectsInput("studyArea", "sf",
                 desc = "Polygon to use as the simulation study area."),
    expectsInput("studyAreaReporting", "sf",
                 desc = paste("multipolygon (typically smaller/unbuffered than `studyArea`)",
                              "to use for plotting/reporting."))
  ),
  outputObjects = bindrows(
    createsOutput("burnDT", "data.table", desc = "data table with pixel IDs of most recent burn"),
    createsOutput("burnMap", "SpatRaster", desc = "cumulative burn map"),
    createsOutput("burnSummary", "data.table", desc = "describes details of all burned pixels"),
    createsOutput("pSpread", "SpatRaster", desc = "spread probability applied to flammability map"),
    createsOutput("rstCurrentBurn", "SpatRaster", desc = "annual burn map")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmSpread = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmSpread", "burn", 7.5)

      if (anyPlotting(P(sim)$.plots)) {
        sim <- scheduleEvent(sim, P(sim)$startTime, "scfmSpread", "plot", 7.5)
      }
    },
    burn = {
      if (!is.null(sim$spreadState)) {

        if (NROW(sim$spreadState[state == "activeSource"]) > 0) {
          sim <- Burnemup(sim) ## fire sizes recorded in  sim$burnSummary
        } else {
          ## make sure to record fires that did not escape/spread
          tempDT <- countBurnedPixelsInSAR(sim$spreadState)
          tempDT$year <- time(sim)
          tempDT[, areaBurned := N * unique(sim$fireRegimePolys$cellSize)]
          tempDT$PolyID <- if (length(tempDT$initialPixels) > 0) sim$fireRegimeRas[tempDT$initialPixels] else NA_integer_
          setnames(tempDT, c("initialPixels"), c("igLoc"))
          sim$burnSummary <- rbind(sim$burnSummary, tempDT)
        }
      }
      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "scfmSpread", "burn", eventPriority = 7.5)
    },
    plot = {
      if (!is.null(sim$rstCurrentBurn)) {
        Plots(sim$rstCurrentBurn, fn = scfmutils::plot_burnMap, type = P(sim)$.plots,
              filename = paste0("currentBurnMap_year_", time(sim)),
              title = paste0("Annual Burn: year ", time(sim)))
        Plots(sim$burnMap, fn = scfmutils::plot_burnMap, type = P(sim)$.plots,
              filename = paste0("cumulativeBurnMap_year_", time(sim)),
              title = paste0("Cumulative Burn: year ", time(sim)))
      }
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmSpread", "plot", eventPriority = 8)
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  compareGeom(sim$rasterToMatch, sim$fireRegimeRas)
  compareGeom(sim$fireRegimeRas, sim$flammableMap)
  tmpRas <- sim$rasterToMatch
  values(tmpRas) <- 1:ncell(tmpRas)
  tmpRas <- postProcess(tmpRas, studyArea = sim$studyAreaReporting)
  mod$pixInSAR <- na.omit(as.vector(tmpRas))

  ## better to use fireRegimeRas than flammableMap, or burnMap inherits attributes
  sim$burnMap <- rast(sim$fireRegimeRas)
  sim$burnMap[!is.na(sim$flammableMap[])] <- 0
  sim$burnMap[sim$flammableMap[] %==% 0] <- NA

  if (!is.null(sim$fireRegimePolys$pSpread)) {
    sprValues <- data.table(PolyID = sim$fireRegimePolys$PolyID,
                            pSpread = sim$fireRegimePolys$pSpread)
    sprRas <- data.table(PolyID = as.vector(sim$fireRegimeRas),
                         flam = as.vector(sim$flammableMap))
    sprValues <- sprValues[sprRas, on = c("PolyID")]
    sprValues[flam != 1, pSpread := NA]
    sim$pSpread <- rast(sim$fireRegimeRas)
    sim$pSpread <- setValues(sim$pSpread, sprValues$pSpread)
  } else {
    warning("using default pSpread as no pSpread column in fireRegimePolys")
    sim$pSpread <- P(sim)$pSpread * sim$flammableMap
  }
  #Create empty data table to store each year's burn data
  sim$burnSummary <- data.table(igLoc = integer(0),
                                grp = integer(0),
                                N = numeric(0),
                                year = numeric(0),
                                areaBurned = numeric(0),
                                PolyID = integer(0))

  sim$rstCurrentBurn <- rast(sim$fireRegimeRas)
  sim$rstCurrentBurn[sim$flammableMap[] %==% 1] <- 0 # reset annual burn
  sim$rstCurrentBurn[sim$flammableMap[] %==% 0] <- NA # might have to ignore warnings

  return(invisible(sim))
}

countBurnedPixelsInSAR <- function(burnDT) {
  tempDT <- copy(burnDT)

  ## grp 1: pixels from fires ignited in SAR & spread in SAR
  ## grp 2: pixels from fires ignited in SAR & spread outside SAR
  ## grp 3: pixels from fires ignited outside SAR & spread in SAR
  ## grp 4: pixels from fires ignited outside SAR & spread outside SAR

  tempDT[(initialPixels %in% mod$pixInSAR) & (pixels %in% mod$pixInSAR), grp := 1L]
  tempDT[(initialPixels %in% mod$pixInSAR) & !(pixels %in% mod$pixInSAR), grp := 2L]
  tempDT[!(initialPixels %in% mod$pixInSAR) & (pixels %in% mod$pixInSAR), grp := 3L]
  tempDT[!(initialPixels %in% mod$pixInSAR) & !(pixels %in% mod$pixInSAR), grp := 4L]
  tempDT <- tempDT[, .N, by = c("initialPixels", "grp")]

  return(tempDT)
}

## name 'Burnemup' is a homage to Walters and Hillborne
Burnemup <- function(sim) {

  threadsDT <- data.table::getDTthreads()
  setDTthreads(1)
  on.exit(setDTthreads(threadsDT), add = TRUE)

  if (!is.null(sim$fireRegimePolys$maxBurnCells)){
    maxSize <- data.table(PolyID = sim$fireRegimePolys$PolyID,
                          maxBurnCells = sim$fireRegimePolys$maxBurnCells)
    burning <- as.vector(sim$fireRegimeRas)[sim$spreadState[state == "inactive"]$initialPixels]
    maxSizes <- maxSize$maxBurnCells[match(burning, maxSize$polyID)]
  } else {
    maxSizes <- sum(as.vector(sim$flammableMap), na.rm = TRUE)
  }

  #TODO: check maxSizes
  sim$burnDT <- SpaDES.tools::spread2(sim$flammableMap,
                                      start = sim$spreadState,
                                      spreadProb = sim$pSpread,
                                      #spreadState = sim$spreadState,
                                      directions = P(sim)$neighbours,
                                      maxSize = maxSizes,
                                      asRaster = FALSE)

  sim$rstCurrentBurn <- rast(sim$fireRegimeRas) ## use fireRegimeRas as template
  sim$rstCurrentBurn[sim$flammableMap[] %==% 1] <- 0 ## reset annual burn
  sim$rstCurrentBurn[sim$flammableMap[] %==% 0] <- NA ## might have to ignore warnings.

  sim$rstCurrentBurn[sim$burnDT$pixels] <- 1 ## update annual burn
  # sim$rstCurrentBurn@data@attributes <- list("Year" == time(sim))

  sim$burnMap[sim$burnDT$pixels] <- sim$burnMap[sim$burnDT$pixels] + 1 ## update cumulative burn

  ## get fire year, pixels burned, area burned, poly ID of all burned pixels in studyAreaReporting
  tempDT <- countBurnedPixelsInSAR(sim$burnDT)
  tempDT$year <- time(sim)
  tempDT$areaBurned <- tempDT$N * unique(sim$fireRegimePolys$cellSize)
  tempDT$PolyID <- if (length(tempDT$initialPixels) > 0) sim$fireRegimeRas[tempDT$initialPixels] else NA_integer_
  setnames(tempDT, c("initialPixels"), c("igLoc"))
  sim$burnSummary <- rbind(sim$burnSummary, tempDT)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  if (!suppliedElsewhere("fireRegimeRas", sim)) {
    stop("you should run scfmLandCoverInit")
  }

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
  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaReporting <- sim$studyArea
  }
  return(sim)
}
