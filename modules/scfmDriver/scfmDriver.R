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
  reqdPkgs = list("fasterize", "PredictiveEcology/LandR@development", "magrittr", "parallel",
                  "PredictiveEcology/pemisc@development", "reproducible", "rgeos",
                  "scam", "sf", "sp", "SpaDES.tools", "stats"),
  parameters = rbind(
    defineParameter("neighbours", "numeric", 8, 4, 8, "number of cell immediate neighbours"),
    defineParameter("buffDist", "numeric", 5e3, 0, 1e5, "Buffer width for fire landscape calibration"),
    defineParameter("pJmp", "numeric", 0.23, 0.18, 0.25, "default spread prob for degenerate polygons"),
    defineParameter("pMin", "numeric", 0.185, 0.15, 0.225, "minimum spread range for calibration"),
    defineParameter("pMax", "numeric", 0.253, 0.24, 0.26, "maximum spread range for calibration"),
    defineParameter("targetN", "numeric", 4000, 1, NA, "target sample size for determining true spread probability"),
    defineParameter("useCloudCache", "logical", getOption("reproducible.useCloud", FALSE), NA, NA,
                    desc = "should a cloud cache be used for heavy operations"),
    defineParameter("cloudFolderID", "character", NULL, NA, NA, "URL for Google-drive-backed cloud cache"),
    defineParameter(".useParallel", class = "logical",
                    default = getOption("pemisc::useParallel", FALSE), min = NA, max = NA,
                    desc = "should driver use parallel? Alternatively accepts a numeric argument, ie how many cores")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "cloudFolderID", "character",
                 paste("URL for Google-drive-backed cloud cache. ",
                       "Note: turn cloudCache on or off with options('reproducible.useCloud')")),
    expectsInput(objectName = "scfmRegimePars", objectClass = "list", desc = ""),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = ""),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame",
                 desc = "shapefile of study area"),
    expectsInput(objectName = "fireRegimePolys", objectClass = "SpatialPolygonsDataFrame",
                 desc = "Areas to calibrate individual fire regime parameters",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "scfmDriverPars", objectClass = "list", desc = "")
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

  if (isTRUE(P(sim)$.useParallel) || P(sim)$.useParallel > 1) { # works if numeric or logical
    maxNumClusters <- length(sim$scfmRegimePars) # this is maximum that is needed
    if (is.numeric(P(sim)$.useParallel)) {
      maxNumClusters <- min(maxNumClusters, P(sim)$.useParallel )  # user may set a smaller maximum passed with P(sim)$.useParallel as a numeric
    }
    cl <- pemisc::makeOptimalCluster(MBper = 5000,
                                     maxNumClusters = maxNumClusters,
                                     outfile = "scfmLog")
    on.exit(try(stopCluster(cl), silent = TRUE))
  } else {
    cl <- NULL
  }

  # Eliot modified this to use cloudCache -- need all arguments named, so Cache works
  sim$scfmDriverPars <- Cache(pemisc::Map2,
                              cl = cl,
                              cloudFolderID = sim$cloudFolderID,
                              useCache = getOption("reproducible.useCache", TRUE),
                              useCloud = getOption("reproducible.useCloud", FALSE),
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
                                              neighbours = P(sim)$neighbours
                              ),
                              f = calibrateFireRegimePolys,
                              userTags = c("mainFunction::Map2", "objectName::scfmDriverPars"))

  names(sim$scfmDriverPars) <- names(sim$scfmRegimePars) #replicate the polygon labels

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using random polygon in Alberta")
    #TODO: remove LandR once this is confirmed working
    studyArea <- LandR::randomStudyArea(size = 1e4*1e6, seed = 23654) #10,000 km * 1000^2m^2
    sim$studyArea <- studyArea
  }

  if (!suppliedElsewhere("fireRegimePolys", sim)) {
    message("fireRegimePolys not supplied. Using default ecoregions of Canada")

    sim$fireRegimePolys <- prepInputs(url = extractURL("fireRegimePolys", sim),
                                      destinationPath = dPath,
                                      studyArea = sim$studyArea,
                                      rasterToMatch = sim$rasterToMatch,
                                      filename2 = TRUE,
                                      overwrite = TRUE,
                                      userTags = c("cacheTags", "fireRegimePolys"))
  }
  return(invisible(sim))
}
