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
  version = list(scfmDriver = "2.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmDriver.Rmd"),
  loadOrder = list(after = c("scfmLandcoverInit", "scfmRegime"),
                   before = c("scfmEscape", "scfmIgnition", "scfmSpread")),
  reqdPkgs = list("parallel", "sf", "spatialEco", "stats", "terra",
                  "PredictiveEcology/LandR (>= 1.1.1)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/scfmutils@development (>= 2.0.7)",
                  "PredictiveEcology/SpaDES.tools (>= 1.0.2.9001)"),
  parameters = rbind(
    defineParameter("buffDist", "numeric", 5e3, 0, 1e5,
                    "Buffer width for fire landscape calibration"),
    defineParameter("cloudFolderID", "character", NULL, NA, NA, "URL for Google-drive-backed cloud cache"),
    defineParameter("dataYear", "numeric", 2011, 1985, 2020,
                    desc = paste("used to select the year of landcover data used to create",
                                 "flammableMapCalibration if the object is unsupplied")),
    defineParameter("pJmp", "numeric", 0.23, 0.18, 0.25, "default spread prob for degenerate polygons"),
    defineParameter("pMax", "numeric", 0.253, 0.24, 0.26, "maximum spread range for calibration"),
    defineParameter("pMin", "numeric", 0.185, 0.15, 0.225, "minimum spread range for calibration"),
    defineParameter("scamOptimizer", "character", "bfgs", NA, NA,
                    "numerical optimization method used in fitting scam model; see `?scam`."),
    defineParameter("targetN", "numeric", 4000, 1, NA, "target sample size for determining true spread probability"),
    defineParameter(".plotInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "simulation time at which the first plot event should occur"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
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
    expectsInput("flammableMapCalibration", "SpatRaster",
                 paste("a flammable map of study area after buffering by `P(sim)$buffDist`.",
                       "Must be supplied by user if `flammableMap` is also supplied.")),
    expectsInput("rasterToMatch", "SpatRaster",
                 "template raster for raster GIS operations. Must be supplied by user.")
  ),
  outputObjects = bindrows(
    createsOutput("fireRegimePolys", "sf",
                  "`fireRegimePolys` with driver attributes appended")
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
        flammableMap = sim$flammableMapCalibration
      ),
      f = scfmutils::calibrateFireRegimePolys
    )
  )

  if (NROW(showCache(userTags = seeIfItHasRun$outputHash)) == 0) {
    cl <- pemisc::makeOptimalCluster(
      useParallel = P(sim)$.useParallelFireRegimePolys,
      ## Estimate as the area of polygon * 2 for "extra" / raster resolution + 400 for fixed costs
      MBper = units::drop_units(sf::st_area(sim$fireRegimePolys)) / prod(res(sim$rasterToMatch)) / 1e3 * 2 + 4e2,
      maxNumClusters = length(unique(sim$fireRegimePolys$PolyID)),
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

  if (!compareGeom(sim$flammableMap, sim$flammableMapCalibration, ext = FALSE, rowcol = FALSE, res = TRUE)) {
    stop("mismatch in resolution of buffered flammable map. Please supply this object manually.")
  }

  message("Running calibrateFireRegimePolys()...")

  flammableMapCalibration <- terra::wrap(sim$flammableMapCalibration)
  scfmDriverPars <- Cache(pemisc::Map2,
                          cl = cl,
                          cloudFolderID = sim$cloudFolderID,
                          ## function-level cache is controlled by option("reproducible.useCache")
                          useCloud = P(sim)$.useCloud,
                          omitArgs = c("cl", "cloudFolderID", "plotPath", "useCache", "useCloud"),
                          polygonType = unique(sim$fireRegimePolys$PolyID),
                          MoreArgs = list(targetN = P(sim)$targetN,
                                          fireRegimePolys = sim$fireRegimePolys,
                                          buffDist = P(sim)$buffDist,
                                          pJmp = P(sim)$pJmp,
                                          pMin = P(sim)$pMin,
                                          pMax = P(sim)$pMax,
                                          flammableMap = flammableMapCalibration,
                                          plotPath = figurePath(sim),
                                          outputPath = outputPath(sim),
                                          optimizer = P(sim)$scamOptimizer
                          ),
                          f = scfmutils::calibrateFireRegimePolys,
                          userTags = c("scfmDriver", "scfmDriverPars"))

  scfmDriverPars <- rbindlist(scfmDriverPars)

  ## drop the attributes if they are present
  colsToDrop <- c("pSpread", "p0", "naiveP0", "pIgnition", "maxBurnCells")
  colsToKeep <- setdiff(names(sim$fireRegimePolys), colsToDrop)
  sim$fireRegimePolys <- sim$fireRegimePolys[colsToKeep]

  sim$fireRegimePolys  <- left_join(sim$fireRegimePolys, scfmDriverPars, by = "PolyID")

 return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(inputPath(sim), 1)

  if (!suppliedElsewhere("fireRegimePolys", sim)) {
    ## it is impossible to get this behaviour correct without also testing for
    ## rasterToMatch, studyArea, and then supplying artificial regime and landcover
    ## attributes anyway.
    stop("fireRegimePolys unsupplied - please run scfmLandcoverInit and scfmRegime")
  }

  if (!suppliedElsewhere("flammableMapCalibration", sim)) {
    bufferedPoly <- st_buffer(sim$fireRegimePolys, (abs(P(sim)$buffDist)))
    bufferedPoly <- fixErrors(bufferedPoly)
    landscapeLCC <- prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      projectTo = sim$rasterToMatch,
      cropTo = bufferedPoly,
      maskTo = bufferedPoly)
    if (!identical(res(landscapeLCC), res(sim$rasterToMatch))) {
      #warning is about identical crs
      landscapeLCC <- suppressWarnings(expr = eval(
        #we want the resolution of rasterToMatch, but not the extent
        Cache(project,
              landscapeLCC,
              method = "near",
              res = res(sim$rasterToMatch),
              y = crs(bufferedPoly),
              userTags = c("scfmDriver", "projectBufferedLCC"))
      ))
    }

    landscapeLCC <- LandR::asInt(landscapeLCC)

    sim$flammableMapCalibration <- defineFlammable(landscapeLCC,
                                             nonFlammClasses = c(20, 31, 32, 33))
  }

  return(invisible(sim))
}
