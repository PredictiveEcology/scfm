defineModule(sim, list(
  name = "scfmDriver",
  description = "generate parameters for the generic percolation model",
  keywords = c("fire"),
  authors = c(
    person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(),
  version = numeric_version("0.1.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmDriver.Rmd"),
  loadOrder = list(after = c("scfmLandcoverInit", "scfmRegime"),
                   before = c("scfmIgnition")),
  reqdPkgs = list("fasterize", "parallel", "sf", "spatialEco", "stats",
                  "PredictiveEcology/LandR@development",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/scfmutils (>= 0.0.7.9001)",
                  "PredictiveEcology/SpaDES.tools (>= 1.0.2.9001)"),
  parameters = rbind(
    defineParameter("buffDist", "numeric", 5e3, 0, 1e5,
                    "Buffer width for fire landscape calibration"),
    defineParameter("bufferLCCYear", "numeric", 2010, NA, 2010,
                    paste("If relying on default buffered flammable map",
                          "the year of LCC to use for defining flammable classes.")),
    defineParameter("cloudFolderID", "character", NULL, NA, NA, "URL for Google-drive-backed cloud cache"),
    defineParameter("pJmp", "numeric", 0.23, 0.18, 0.25, "default spread prob for degenerate polygons"),
    defineParameter("pMax", "numeric", 0.253, 0.24, 0.26, "maximum spread range for calibration"),
    defineParameter("pMin", "numeric", 0.185, 0.15, 0.225, "minimum spread range for calibration"),
    defineParameter("scamOptimizer", "character", "bfgs", NA, NA,
                    "numerical optimization method used in fitting scam model; see `?scam`."),
    defineParameter("targetN", "numeric", 4000, 1, NA, "target sample size for determining true spread probability"),
    defineParameter(".plotInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".useCache", "character", FALSE, NA, NA,
                    "Can be names of events or the whole module name; these will be cached by SpaDES"),
    defineParameter(".useCloud", "logical", getOption("reproducible.useCloud", FALSE), NA, NA,
                    "should a cloud cache be used for heavy operations"),
    defineParameter(".useParallelFireRegimePolys", "logical", getOption("pemisc.useParallel", FALSE), NA, NA,
                    "should driver use parallel? Alternatively accepts a numeric argument, i.e., how many cores.")
  ),
  inputObjects = bindrows(
    expectsInput("cloudFolderID", "character",
                 paste("URL for Google-drive-backed cloud cache.",
                       "Note: turn `cloudCache` on or off with `options('reproducible.useCloud')`.")),
    expectsInput("fireRegimePolys", "sf",
                 paste("Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
                       "Must have numeric field 'PolyID' or it will be created for individual polygons.")),
    expectsInput("flammableMapLarge", "SpatRaster",
                 paste("a flammable map of study area after buffering by `P(sim)$buffDist`.",
                       "Defaults to LCC2010. Must be supplied by user if `flammableMap` is also supplied.")),
    expectsInput("rasterToMatch", "SpatRaster",
                 "template raster for raster GIS operations. Must be supplied by user.")
  ),
  outputObjects = bindrows(
    createsOutput("fireRegimePolys", "sf", "fireRegimePolys with driver attributes appended")
  )
))

## event types
#   - type `init` is required for initilization

doEvent.scfmDriver = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
    },

    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  if (is(sim$fireRegimePolys, "SpatialPolygonsDataFrame")) {
    sim$fireRegimePolys <- st_as_sf(sim$fireRegimePolys)
  }

  ## Check to see if it is a Cache situation -- if it is, don't make a cl -- on Windows, takes too long
  seeIfItHasRun <- CacheDigest(
    list(
      Map2,
      polygonType = unique(sim$fireRegimePolys$PolyID),
      MoreArgs = list(
        targetN = P(sim)$targetN,
        fireRegimePolys = sim$fireRegimePolys,
        buffDist = P(sim)$buffDist,
        pJmp = P(sim)$pJmp,
        pMin = P(sim)$pMin,
        pMax = P(sim)$pMax,
        flammableMap = sim$flammableMapLarge
      ),
      f = calibrateFireRegimePolys ## scfmutils
    )
  )

  if (NROW(showCache(userTags = seeIfItHasRun$outputHash)) == 0) {
    cl <- pemisc::makeOptimalCluster(
      useParallel = P(sim)$.useParallelFireRegimePolys,
      ## Estimate as the area of polygon * 2 for "extra" / raster resolution + 400 for fixed costs
      MBper = units::drop_units(sf::st_area(sim$fireRegimePolys)) / prod(res(sim$rasterToMatch)) / 1e3 * 2 + 4e2,
      maxNumClusters = length(sim$scfmRegimePars),
      outfile = file.path(outputPath(sim), "log", "scfm.log"),
      objects = c(), envir = environment(),
      libraries = c("scfmutils")
    )

    on.exit({
      if (!is.null(cl)) {
        parallel::stopCluster(cl)
      }
    })
  } else {
    cl <- NULL
  }

  if (!compareGeom(sim$flammableMap, sim$flammableMapLarge, ext = FALSE, rowcol = FALSE, res = TRUE)) {
    stop("mismatch in resolution of buffered flammable map. Please supply this object manually.")
  }

  message("Running calibrateFireRegimePolys()...")

  flammableMapLarge <- terra::wrap(sim$flammableMapLarge)
  scfmDriverPars <- Cache(pemisc::Map2,
                              cl = cl,
                              cloudFolderID = sim$cloudFolderID,
                              #function-level cache is controlled by option("reproducible.useCache")
                              useCloud = P(sim)$.useCloud,
                              omitArgs = c("cl", "cloudFolderID", "plotPath", "useCache", "useCloud"),
                              polygonType = unique(sim$fireRegimePolys$PolyID),
                              MoreArgs = list(targetN = P(sim)$targetN,
                                              fireRegimePolys = sim$fireRegimePolys,
                                              buffDist = P(sim)$buffDist,
                                              pJmp = P(sim)$pJmp,
                                              pMin = P(sim)$pMin,
                                              pMax = P(sim)$pMax,
                                              flammableMap = flammableMapLarge,
                                              plotPath = file.path(outputPath(sim), "figures"),
                                              outputPath = file.path(outputPath(sim)),
                                              optimizer = P(sim)$scamOptimizer
                              ),
                              f = scfmutils::calibrateFireRegimePolys,
                              userTags = c("scfmDriver", "scfmDriverPars"))


  scfmDriverPars <- rbindlist(scfmDriverPars)
  # names(sim$scfmDriverPars) <- names(sim$scfmRegimePars) #replicate the polygon labels
  sim$fireRegimePolys  <- left_join(sim$fireRegimePolys, scfmDriverPars, by = "PolyID")

 return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (!suppliedElsewhere("flammableMapLarge")) {
    bufferedPoly <- st_buffer(sim$fireRegimePolys, (abs(P(sim)$buffDist)))
    bufferedPoly <- fixErrors(bufferedPoly)
    landscapeLCC <- Cache(prepInputsLCC,
                          year = P(sim)$bufferLCCYear,
                          destinationPath = dataPath(sim),
                          studyArea = bufferedPoly,
                          projectTo = sim$rasterToMatch,
                          cropTo = bufferedPoly,
                          maskTo = bufferedPoly,
                          method = "near",
                          omitArgs = "destinationPath")

    landscapeLCC <- setValues(landscapeLCC, as.integer(values(landscapeLCC)))

    if (P(sim)$bufferLCCYear == 2010 | P(sim)$bufferLCCYear == 2015) {
      nonFlamClasses <- c(13L, 16L, 17L, 18L, 19L)
    } else if (P(sim)$bufferLCCYear == 2005) {
      nonFlamClasses <- c(0L, 25L, 30L, 33L, 36L, 37L, 38L, 39L)
    } else {
      stop("invalid bufferLCCYear - please supply flammableMapLarge")
    }
    sim$flammableMapLarge <- defineFlammable(landscapeLCC, nonFlammClasses = nonFlamClasses)
  }

  return(invisible(sim))
}
