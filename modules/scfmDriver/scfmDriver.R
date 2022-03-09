defineModule(sim, list(
  name = "scfmDriver",
  description = "generate parameters for the generic percolation model",
  keywords = c("fire"),
  authors = c(person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre")),
              person("Ian", "Eddy", email = "ian.eddy@canada.ca", role = c("aut"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmDriver.Rmd"),
  reqdPkgs = list("fasterize", "PredictiveEcology/LandR", "magrittr", "parallel",
                  "PredictiveEcology/pemisc", "reproducible",
                  "scam", "sf", "SpaDES.tools", "stats", "spatialEco"), #scam  (==1.2.3)
  parameters = rbind(
    defineParameter("neighbours", "numeric", 8, 4, 8, "number of cell immediate neighbours"),
    defineParameter("buffDist", "numeric", 5e3, 0, 1e5, "Buffer width for fire landscape calibration"),
    defineParameter("bufferLCCYear", "numeric", 2010, NA, 2010,
                    desc = paste("If relying on default buffered flammable map",
                                 "the year of LCC to use for defining flammable classes")),
    defineParameter("pJmp", "numeric", 0.23, 0.18, 0.25, "default spread prob for degenerate polygons"),
    defineParameter("pMin", "numeric", 0.185, 0.15, 0.225, "minimum spread range for calibration"),
    defineParameter("pMax", "numeric", 0.253, 0.24, 0.26, "maximum spread range for calibration"),
    defineParameter("targetN", "numeric", 4000, 1, NA, "target sample size for determining true spread probability"),
    defineParameter("cloudFolderID", "character", NULL, NA, NA, "URL for Google-drive-backed cloud cache"),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name; these will be cached by SpaDES"),
    defineParameter(".useCloud", "logical", getOption("reproducible.useCloud", FALSE), NA, NA,
                    desc = "should a cloud cache be used for heavy operations"),
    defineParameter(".useParallel", class = "logical",
                    default = getOption("pemisc::useParallel", FALSE), min = NA, max = NA,
                    desc = "should driver use parallel? Alternatively accepts a numeric argument, ie how many cores")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "flammableMapLarge", "RasterLayer",
                 desc = paste("a flammable map of study area after buffering by P(sim)$buffDist.",
                              "Defaults to LCC2010. Must be supplied by user flammableMap is also supplied")),
    expectsInput(objectName = "cloudFolderID", "character",
                 paste("URL for Google-drive-backed cloud cache. ",
                       "Note: turn cloudCache on or off with options('reproducible.useCloud')")),
    expectsInput(objectName = "fireRegimePolys", "list",
                 desc = paste("Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
                              "Must have numeric field 'PolyID' or it will be created for individual polygons")),
    expectsInput(objectName = "rasterToMatch", "RasterLayer",
                 desc = "template raster for raster GIS operations. Must be supplied by user"),
    expectsInput(objectName = "scfmRegimePars", objectClass = "list",
                 desc = "list of fire regime parameters for each polygon"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list",
                 desc = "contains landscape attributes for each polygon")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "scfmDriverPars", objectClass = "list",
                  desc = "burn parameters for each polygon in fireRegimePolys")
  )
))

## event types
#   - type `init` is required for initiliazation

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

# 1 - (1-p0)**N = pEscape
# 1 - pEscape = (1-p0)**N
# (1 - pEscape)**1/N = 1 - p0
# p0 = 1 - (1 - pEscape)**1/N

hatP0 <- function(pEscape, n = 8) {
  1 - (1 - pEscape) ** (1 / n)
}

#a real clever boots would minimise the abs log odds ratio.
#be my guest.

escapeProbDelta <- function(p0, w, hatPE) {
  abs(sum(w*(1 - (1 - p0) ** (0:8))) - hatPE)
}

Init <- function(sim) {
  cellSize <- sim$landscapeAttr[[1]]$cellSize

  # Download 1 canonical version of the LCC, cropped to the sim$fireRegimePolys + buffer,
  #  pass this one into the calibrateFireRegimePolys, avoiding many downloads (esp when
  #  in parallel)

  # Check to see if it is a Cache situation -- if it is, don't make a cl -- on Windows, takes too long
  seeIfItHasRun <- CacheDigest(list(pemisc::Map2,
                                    regime = sim$scfmRegimePars,
                                    polygonType = names(sim$scfmRegimePars),
                                    MoreArgs = list(targetN = P(sim)$targetN,
                                                    landAttr = sim$landscapeAttr,
                                                    cellSize = cellSize,
                                                    fireRegimePolys = sim$fireRegimePolys,
                                                    buffDist = P(sim)$buffDist,
                                                    pJmp = P(sim)$pJmp,
                                                    pMin = P(sim)$pMin,
                                                    pMax = P(sim)$pMax,
                                                    neighbours = P(sim)$neighbours,
                                                    flammableMap = sim$flammableMapLarge
                                    ),
                                    calibrateFireRegimePolys))

  if (NROW(showCache(userTags = seeIfItHasRun$outputHash)) == 0) {
    cl <- pemisc::makeOptimalCluster(
      useParallel = P(sim)$.useParallel,
      # Estimate as the area of polygon * 2 for "extra" / raster resolution + 400 for fixed costs
      MBper = rgeos::gArea(sim$fireRegimePolys)/(prod(res(sim$landscapeLCC)))/1e3 * 2 + 4e2, # in MB
      maxNumClusters = length(sim$scfmRegimePars),
      outfile = "scfmLog",
      objects = c("genSimLand"), envir = environment(),
      libraries = c("rlang", "raster", "rgeos", "reproducible",
                    "LandR", "sf", "fasterize", "data.table")
    )

    on.exit({
      if (!is.null(cl))
        parallel::stopCluster(cl)
    })
  } else {
    cl <- NULL
  }

  if (!identical(res(sim$flammableMap), res(sim$flammableMapLarge))) {
    stop("mismatch in resolution of buffered flammable map. Please supply this object manually")
  }

  message("Running calibrateFireRegimePolys()...")
  sim$scfmDriverPars <- Cache(pemisc::Map2,
                              cl = cl,
                              cloudFolderID = sim$cloudFolderID,
                              useCache = P(sim)$.useCache, #getOption("reproducible.useCache", TRUE),
                              useCloud = P(sim)$.useCloud,
                              omitArgs = c("useCloud", "useCache", "cloudFolderID", "cl"),
                              regime = sim$scfmRegimePars,
                              polygonType = names(sim$scfmRegimePars),
                              MoreArgs = list(targetN = P(sim)$targetN,
                                              landAttr = sim$landscapeAttr,
                                              cellSize = cellSize,
                                              fireRegimePolys = sim$fireRegimePolys,
                                              buffDist = P(sim)$buffDist,
                                              pJmp = P(sim)$pJmp,
                                              pMin = P(sim)$pMin,
                                              pMax = P(sim)$pMax,
                                              neighbours = P(sim)$neighbours,
                                              flammableMap = sim$flammableMapLarge
                              ),
                              calibrateFireRegimePolys,
                              userTags = c("scfmDriver", "scfmDriverPars"))

  names(sim$scfmDriverPars) <- names(sim$scfmRegimePars) #replicate the polygon labels
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

  if (any(!suppliedElsewhere("scfmRegimePars", sim),
          !suppliedElsewhere("landscapeAttr", sim))) {
    stop("this module cannot be run without scfmRegime and scfmLandcoverInit")
  }

  if (!suppliedElsewhere("flammableMapLarge")) {
    bufferedPoly <- buffer(sim$fireRegimePolys, (abs(P(sim)$buffDist)))
    bufferedPoly <- fixErrors(bufferedPoly)
    landscapeLCC <- Cache(prepInputsLCC,
                          year = P(sim)$bufferLCCYear,
                          destinationPath = dataPath(sim),
                          studyArea = bufferedPoly, useSAcrs = TRUE,
                          omitArgs = "destinationPath")
    if (!identical(res(landscapeLCC), res(sim$rasterToMatch))) {
      #warning is about identical crs
     landscapeLCC <- suppressWarnings(expr = eval(
        #we want the resolution of rasterToMatch, but not the extent
        Cache(projectRaster,
              landscapeLCC,
              method = "ngb",
              res = res(sim$rasterToMatch),
              crs = crs(bufferedPoly),
              userTags = c("scfmDriver", "projectBufferedLCC"))
      ))
    }

    landscapeLCC <- setValues(landscapeLCC, as.integer(getValues(landscapeLCC)))

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
