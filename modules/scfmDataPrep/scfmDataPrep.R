defineModule(sim, list(
  name = "scfmDataPrep",
  description = paste("This module first generates some relevant fire regime statistics for each fire regime over
                      `studyAreaCalibration` then filters the fire regime polys to those inside studyArea.",
                      "It will combine fire regime polys with smaller area than that denoted by the `sliverThreshold`",
                      "param before calculating the flammable area, mean fire size, maximum fire size,",
                      "number of flammable neighbouring pixels from 0-8, and lastly, the ignition rate, escape rate,
                      and spread rate in each polygon. By default these estimates are based on lightning-caused fires",
                      "from 1970-2000 in the NFDB dataset. However, these params can be overriden by a user.",
                      "The FRI can be set using the `targetBurnRate` param, in which case the mean fire size, ignition rate",
                      "and escape rate will be incrementally adjusted to match the target FRI. An important limitation",
                      "is that all spatial objects must share the same CRS and resolution, where relevant,",
                      "and they must utilize a crs projected in metres"),
  keywords =  c("fire regime", "fire percolation model", "National Fire Data Base (NFBD)"),
  authors =  c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(scfmDataPrep = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "scfmDataPrep.Rmd"),
  reqdPkgs = list(
    "dplyr", "ggplot2", "parallel",
    "PredictiveEcology/LandR (>= 1.1.1)",
    "PredictiveEcology/pemisc@development",
    "PredictiveEcology/scfmutils@development (>= 2.0.7)",
    "PredictiveEcology/SpaDES.core@development (>= 2.1.5.9002)",
    "PredictiveEcology/SpaDES.tools (>= 1.0.2.9001)",
    "purrr", "reproducible", "sf", "stats", "terra"),
  parameters = rbind(
    defineParameter("buffDist", "numeric", 2e4, 1, 1e5,
                    paste("Buffer width to mitigate edge effects in fire landscape calibration.",
                          "If studyAreaCalibration is not supplied, this parameter will also be",
                          "used to create it via buffering studyArea")),
    defineParameter("cloudFolderID", "character", NULL, NA, NA, "URL for Google-drive-backed cloud cache"),
    defineParameter("dataYear", "numeric", 2011, 1985, 2020,
                    desc = paste("used to select the year of landcover data used to create",
                                 "flammableMap if the obejct is unsupplied")),
    defineParameter("empiricalMaxSizeFactor", "numeric", 1.2, 1, 10,
                    desc = "scale `xMax` by this if HD estimator fails"),
    defineParameter("eventsToPrepare", "character", c("scfmLandcoverInit", "scfmRegime", "scfmDriver"), NA, NA,
                    paste("which of three possible dataPrep steps to run? scfmLandcoverInit will clean GIS",
                          "and generate landscape stats regarding flammability in each fire regime poly;",
                          "scfmRegime will prepare fire regime attributes inc. mean fire size and ignition rate;",
                          "scfmDriver will estimate the spread probability of flammable pixels")),
    defineParameter("fireCause", "character", c("N", "L"), NA_character_, NA_character_,
                    desc = "subset of `c('H', 'H-PB', 'N', 'Re', 'U')`"),
    defineParameter("fireCauseColumnName", "character", "CAUSE", NA, NA,
                    desc = "Name of the column that has fire cause, consistent with `P(sim)$fireCause`."),
    defineParameter("fireEpoch", "numeric", c(1971, 2020), NA, NA, "start of normal period"),
    defineParameter("fireRegimePolysType", "character", "ECOREGION", NA, NA,
                    paste("Polygon type to use for scfm `fireRegimePolys`:",
                          "see `?scfmutils::prepInputsFireRegimePolys` for allowed types.")),
    defineParameter("fireSizeColumnName", "character", "SIZE_HA", NA, NA,
                    desc = "Name of the column that has fire size"),
    defineParameter("fireYearColumnName", "character", "YEAR", NA, NA,
                    desc = "Name of the column that has fire size"),
    defineParameter("flammabilityThreshold", "numeric", 0.25, 0, 1,
                    paste("Minimum proportion of flammable old pixel needed to define a new pixel
                          as flammable when upscaling the default flammable maps`.")),
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of immediate cell neighbours"),
    defineParameter("pJmp", "numeric", 0.23, 0.18, 0.25, "default spread prob for degenerate polygons"),
    defineParameter("pMax", "numeric", 0.253, 0.24, 0.26, "maximum spread range for calibration"),
    defineParameter("pMin", "numeric", 0.185, 0.15, 0.225, "minimum spread range for calibration"),
    defineParameter("scamOptimizer", "character", "bfgs", NA, NA,
                    "numerical optimization method used in fitting scam model; see `?scam`."),
    defineParameter("sliverThreshold", "numeric", 6.25e8, NA, NA,
                    paste("fire regime polygons with area (in m2) less than this number will be merged",
                          "with their closest non-sliver neighbour using `sf::st_nearest_feature`.")),
    defineParameter("targetBurnRate", "numeric", NA, 0, 1,
                    desc = paste("a named vector giving the proportional annual area burned of each fire regime polygon.",
                                 "These override the default estimate of scfm and are used to estimate a new mean",
                                 "fire size and ignition rate. Names should correspond to `PolyID`.",
                                 "A partial set of polygons is allowed - missing polys are estimated from data.")),
    defineParameter("targetMaxFireSize", "numeric", NA, 0, NA,
                    desc = paste("a named vector giving the estimated max fire size (in $ha$) of each fire regime polygon.",
                                 "These will override the default estimate of scfm and will be used to estimate",
                                 "a new spread probability. Names should correspond to `PolyID`.",
                                 "A partial set of polygons is allowed - missing polys are estimated from data.")),
    defineParameter("targetN", "numeric", 4000, 1, NA, "target sample size for determining true spread probability"),

    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "Initial time for plotting"),
    defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, "Interval between plotting"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by `Plots` function, which can be optionally used here."),
    defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, "Initial time for saving"),
    defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, "Interval between save events"),
    defineParameter(".useCache", "character", ".inputObjects", NA, NA,
                    "Use caching of events - not recommended as of 10/05/2023"),
    defineParameter(".useCloud", "logical", getOption("reproducible.useCloud", FALSE), NA, NA,
                    "should a cloud cache be used for heavy operations"),
    defineParameter(".useParallelFireRegimePolys", "logical", getOption("pemisc.useParallel", FALSE), NA, NA,
                    "should driver use parallel? Alternatively accepts a numeric argument, i.e., how many cores.")
  ),
  inputObjects = bindrows(
    expectsInput("cloudFolderID", "character",
                 paste("URL for Google-drive-backed cloud cache.",
                       "Note: turn `cloudCache` on or off with `options('reproducible.useCloud')`.")),
    expectsInput("firePoints", "sf",
                 desc = paste0("Historical fire data in point form. Must contain fields 'CAUSE',
                               'YEAR', and 'SIZE_HA', or pass the parameters to identify those."),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput("fireRegimePolys", "sf",
                 desc = paste("Areas to calibrate individual fire regime parameters.",
                              "Defaults to ecozones of Canada.",
                              "Must have numeric field 'PolyID' or it will be created for individual polygons.")),
    expectsInput("fireRegimePolysCalibration", "sf",
                 desc = paste("if `studyAreaCalibration` is supplied, the corresponding fire regime areas.",
                              "Requires integer field `PolyID` if supplied. Uses same defaults as `fireRegimePolys`.")),
    expectsInput("flammableMap", "SpatRaster",
                 desc = "binary flammability map - defaults to using LandR::prepInputsLCC"),
    expectsInput("flammableMapCalibration", "SpatRaster",
                 desc = paste("binary flammability map corresponding to `rasterToMatchCalibration`.",
                              "It should extent from studyArea by >= scfmDriver's `P(sim)$buffDist`.",
                              "and if unsupplied, will be created using `LandR::prepInputs_NTEMS_LCC_FAO")),
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = "template raster for raster GIS operations. Must be supplied by user"),
    expectsInput("rasterToMatchCalibration", "SpatRaster",
                 desc = paste("Template raster for studyAreaCalibration - will be created based on rasterToMatch if unsupplied")),
    expectsInput("studyArea", "sf", desc = "Polygon to use as the simulation study area (typically buffered)."),
    expectsInput("studyAreaCalibration", "sf", desc = "optional larger study area used for parameterization only")
  ),
  outputObjects = bindrows(
    createsOutput("fireRegimePoints", "sf",
                  desc = "Fire locations. These are filtered according to criteria set in params (i.e. epoch, cause)"),
    createsOutput("fireRegimePolys", "sf",
                  desc = "`fireRegimePolys` with fire attributes appended."),
    createsOutput("fireRegimePolysCalibration", "sf",
                  desc = "`fireRegimePolysCalibration` with attributes appended"),
    createsOutput("fireRegimeRas", "SpatRaster",
                  desc = "Rasterized version of fireRegimePolys with values representing polygon ID")
  )
))

doEvent.scfmDataPrep = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {

      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmDataPrep", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "scfmDataPrep", "save")
    },
    #
    plot = {
      ## NOTE: these objects don't change during sim, so only need to be plotted once
      flamRegime <- sim$fireRegimeRas
      flamRegime[sim$flammableMap[] == 0] <- NA
      Plots(flamRegime, fn = scfmutils::plot_fireRegimeRas, type = P(sim)$.plots,
            filename = paste0("flam_fireRegimeRas"),
            title = paste0("Fire regimes"))
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

### template initialization
Init <- function(sim) {

  if ("scfmLandcoverInit" %in% P(sim)$eventsToPrepare){
    sim <- prepare_scfmLandcoverInit(sim)
  }
  if ("scfmRegime" %in% P(sim)$eventsToPrepare){
    sim <- prepare_scfmRegime(sim)
  }
  if ("scfmDriver" %in% P(sim)$eventsToPrepare){
    sim <- prepare_scfmDriver(sim)
  }
  return(invisible(sim))
}


prepare_scfmLandcoverInit <- function(sim) {

  ## ensure flammability maps are integer ('binary') maps
  if (!LandR::isInt(sim$flammableMap)) {
    sim$flammableMap <- LandR::asInt(sim$flammableMap)
  }

  if (!is.integer(sim$flammableMapCalibration[])) {
    sim$flammableMapCalibration <- LandR::asInt(sim$flammableMapCalibration)
  }

  stopifnot(
    all(unique(sim$flammableMap[]) %in% c(NA_integer_, 0L, 1L)),
    all(unique(sim$flammableMapCalibration[]) %in% c(NA_integer_, 0L, 1L))
  )

  message("checking sim$fireRegimePolys for sliver polygons...")
  # this only needs to be done on the larger area, if it is provided
  # doing so on larger and smaller has the potential to
  # mismatch slivers between calibration/simulation

  sim$fireRegimePolysCalibration <- checkForIssues(
    fireRegimePolys = sim$fireRegimePolysCalibration,
    studyArea = sim$studyAreaCalibration,
    rasterToMatch = sim$rasterToMatchCalibration,
    flammableMap = sim$flammableMapCalibration,
    sliverThresh = P(sim)$sliverThreshold,
    cacheTag = c("scfmLandcoverInit", "fireRegimePolysCalibration")
  )

  ## now that slivers are removed, remake frp from the larger object
  sim$fireRegimePolys <- postProcess(sim$fireRegimePolysCalibration,
                                     studyArea = sim$studyArea)
  ## for now - GIS operations with sf objects are causing sliver polygons (area < 0.001 m2)

  if (is(st_geometry(sim$fireRegimePolys), "sfc_GEOMETRY")) {
    # this object may have empty geometries, which can occur when SAC and SA are both subsets
    # of the same file. the empty geometries will cause an error.
    sim$fireRegimePolys <- sim$fireRegimePolys[as.numeric(st_area(sim$fireRegimePolys)) > 0, ]
    #in the event this results in LINESTRING or POINT objects,remove them to prevent error
    sim$fireRegimePolys <- st_collection_extract(sim$fireRegimePolys, "POLYGON")
    sim$fireRegimePolys <- st_cast(sim$fireRegimePolys, "MULTIPOLYGON")
  }

  temp <- sim$fireRegimePolysCalibration[order(sim$fireRegimePolysCalibration$PolyID), ]
  sim$fireRegimePolysCalibration <- temp #to fit on two lines easily
  sim$fireRegimePolysCalibration <- Cache(genFireMapAttr,
                                          flammableMap = sim$flammableMapCalibration,
                                          fireRegimePolys = sim$fireRegimePolysCalibration,
                                          neighbours = P(sim)$neighbours,
                                          userTags = c(currentModule(sim),
                                                       "genFireMapAttr",
                                                       "studyAreaCalibration")
  )


  sim$fireRegimePolys <- checkForIssues(
    fireRegimePolys = sim$fireRegimePolys,
    studyArea = sim$studyArea,
    rasterToMatch = sim$rasterToMatch,
    flammableMap = sim$flammableMap,
    sliverThresh = P(sim)$sliverThreshold,
    cacheTag = c("scfmLandcoverInit", "fireRegimePolys")
  )
  sim$fireRegimePolys <- sim$fireRegimePolys[order(sim$fireRegimePolys$PolyID),]

  sim$fireRegimePolys <- Cache(genFireMapAttr,
                               flammableMap = sim$flammableMap,
                               fireRegimePolys = sim$fireRegimePolys,
                               neighbours = P(sim)$neighbours,
                               userTags = c(currentModule(sim), "genFireMapAttr", "studyArea")
  )

  ## doing this prevents fireRegimeRas from inheriting colormaps
  sim$fireRegimeRas <- rasterize(sim$fireRegimePolys, sim$rasterToMatch, fun = "max", field = "PolyID")
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

prepare_scfmRegime <- function(sim) {

  tmp <- sim$firePoints

  ## extract and validate fireCause spec
  fc <- P(sim)$fireCause

  ## review that sf can be used like this.
  ## should verify CAUSE is a column in the table...
  if (is.factor(tmp[[P(sim)$fireCauseColumnName]])) {
    causeSet <- levels(tmp[[P(sim)$fireCauseColumnName]])
  } else {
    causeSet <- unique(tmp[[P(sim)$fireCauseColumnName]])
  }

  if ("N" %in% fc & "L" %in% causeSet) fc[fc == "N"] <- "L"
  if ("L" %in% fc & "N" %in% causeSet) fc[fc == "L"] <- "N"

  if (all(!(fc %in% causeSet))) {
    notPresent <- fc[!fc %in% causeSet]
    warning(paste0("This firecause is not present: ", notPresent,
                   " The following are the fire causes: ",
                   paste(causeSet, collapse = ", "),
                   ". Original cause will be replaced by ",
                   paste(causeSet, collapse = ", ")), immediate. = TRUE)
    fc <- causeSet
  }

  tmp <- subset(tmp, get(P(sim)$fireCauseColumnName) %in% fc)

  #extract and validate fireEpoch
  epoch <- P(sim)$fireEpoch
  if (length(epoch) != 2 || !is.numeric(epoch) || any(!is.finite(epoch)) || epoch[1] > epoch[2]) {
    stop("illegal fireEpoch: ", epoch)
  }

  quotes <- paste0("tmp$", paste(eval(P(sim)$fireYearColumnName)))
  tmp <- subset(tmp, get(P(sim)$fireYearColumnName) >= epoch[1] &
                  get(P(sim)$fireYearColumnName) <= epoch[2])

  epochLength <- as.numeric(epoch[2] - epoch[1] + 1)

  if (sf::st_crs(tmp) != sf::st_crs(sim$fireRegimePolysCalibration)) {
    tmp <- sf::st_transform(tmp, crs = sf::st_crs(sim$fireRegimePolysCalibration))
  }

  tmp <- sf::st_intersection(tmp, sim$fireRegimePolysCalibration) ## gives studyArea colnames to points

  if (any(is.na(tmp$PolyID))) {
    tmp <- tmp[!is.na(tmp$PolyID), ] ## need to remove NA points
  }
  sim$fireRegimePoints <- tmp

  ## this function estimates the ignition probability and escape probability based on NFDB
  scfmRegimePars <- unique(sim$fireRegimePolysCalibration$PolyID) |>
    lapply(
      FUN = calcZonalRegimePars,
      firePolys = sim$fireRegimePolysCalibration,
      firePoints = sim$fireRegimePoints,
      epochLength = epochLength,
      maxSizeFactor = P(sim)$empiricalMaxSizeFactor,
      fireSizeColumnName = P(sim)$fireSizeColumnName,
      targetBurnRate = P(sim)$targetBurnRate,
      targetMaxFireSize = P(sim)$targetMaxFireSize
    ) |>
    rbindlist(fill = TRUE)

  ## drop the attributes if they are present
  colsToDrop <- c("ignitionRate", "pEscape", "xBar", "lxBar",
                  "xMax", "emfs_ha", "empiricalBurnRate")
  colsToKeep <- setdiff(names(sim$fireRegimePolys), colsToDrop)
  sim$fireRegimePolys <- sim$fireRegimePolys[colsToKeep]

  ## only keep the attributes that are in study area
  sim$fireRegimePolys <- left_join(sim$fireRegimePolys, scfmRegimePars, by = "PolyID")

  return(invisible(sim))
}

prepare_scfmDriver <- function(sim) {

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

  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)

  # object check for SA/FRP/FRPC/SAC - better to be strict with stops
  hasSA <- suppliedElsewhere("studyArea", sim)
  hasSAC <- suppliedElsewhere("studyAreaCalibration", sim)
  hasFRP <- suppliedElsewhere("fireRegimePolys", sim)
  hasFRPC <- suppliedElsewhere("fireRegimePolysCalibration", sim)
  hasRTM <- suppliedElsewhere("rasterToMatch", sim)
  hasRTMC <- suppliedElsewhere("rasterToMatchCalibration", sim)
  hasFM <- suppliedElsewhere("flammableMap", sim)
  hasFMC <- suppliedElsewhere("flammableMapCalibration", sim)

  if (c(hasFRP & !hasFRPC) | c(hasFM & !hasFMC)) {
    stop("if supplying flammableMap or fireRegimePolys",
         "the equivalent calibration-sized object must also be provided")
  }

  # supply objects
  if (!hasSA) {
    message("study area not supplied. Using random polygon in Alberta")
    sim$studyArea <- LandR::randomStudyArea(size = 1500000* 1000, seed = 23654)
    sim$studyArea <- terra::project(sim$studyArea, y = "EPSG:3348")
    #this is 1,500,000 km2 - somewhere in eastern Rockies
    #the crs is Canada equal alberts - unfortunately there is no way to set

  }

  if (!hasSAC){
    # buffDist is necessary only to ensure fires aren't extinguished from edges
    # during the spread calibration - whereas the buffer distance here is to establish
    # studyAreaCalibration, whihc is intended to provide additional fire data for
    # fire regime polygons that are otherwise too small after intersecting with studyArea.
    # however - this distance must logically exceed P(sim)$buffDist
    #ideally it is larger than the sqrt(max(sim$firePoints$SIZE_HA))
    sim$studyAreaCalibration <- buffer(sim$studyArea, P(sim)$buffDist * 2)
  }
  if (hasRTM & !hasRTMC) {
    warning("rasterToMatchCalibration not supplied")
    sim$rasterToMatchCalibration <- terra::extend(sim$rasterToMatch,
                                                  sim$studyAreaCalibration,
                                                  fill = 1L)
  }

  if (!hasRTM & !hasRTMC) {
    warning(paste(
      "rasterToMatch not supplied. generating from NTEMS LCC",
      " - It is strongly recommended to supply a rasterToMatch"
    ))
    sim$rasterToMatchCalibration <- LandR::prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      cropTo = sim$studyAreaCalibration,
      projectTo = sim$studyAreaCalibration,
      maskTo= sim$studyAreaCalibration,
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatchCalibration")
    )
    sim$rasterToMatch <- postProcess(sim$rasterToMatchCalibration,
                                     to = sim$studyArea)
  }

  if (!hasRTM & hasRTMC) {
    sim$rasterToMatch <- postProcess(sim$rasterToMatchCalibration,
                                     to = sim$studyArea)
  }

  if (!hasFMC) {
    fmc <- prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      maskTo = sim$studyAreaCalibration,
      cropTo = sim$rasterToMatchCalibration,
      #projectTo = sim$rasterToMatchCalibration, #should be done after defineFlammable
      userTags = c("prepInputs_NTEMS_LCC_FAO", "studyArea")
    )
    fmc[] <- asInteger(fmc[])
    fmc <- defineFlammable(fmc,
                           nonFlammClasses = c(20, 31, 32, 33))

    sim$flammableMapCalibration <- postProcess(fmc,
                                               to = sim$rasterToMatchCalibration,
                                               method = "mode")
  }

  if (!hasFM) {
    sim$flammableMap <- postProcess(sim$flammableMapCalibration, to = sim$rasterToMatch,
                                    method = "near")
  }

  ## this is TRUE unless fireRegimePolysCalibration is supplied, in which case we drop that object
  if (!hasFRPC) {

    message("fireRegimePolys not supplied. Using default ecoregions of Canada")
    # cannot use prepInputs with a vector for prepInputs - unreliable w/ GDAL

    sim$fireRegimePolysCalibration <- Cache(prepInputsFireRegimePolys, url = NULL,
                                            destinationPath = dPath,
                                            studyArea = sim$studyAreaCalibration,
                                            type = P(sim)$fireRegimePolysType)
  }

  if (!hasFRP) {
    #avoid GIS issue with sf
    sim$fireRegimePolys <- postProcess(terra::vect(sim$fireRegimePolysCalibration),
                                       to = sim$studyArea) |>
      sf::st_as_sf()
  }

  if (!suppliedElsewhere("firePoints", sim)) {
    ## NOTE: do not use fireSenseUtils - it removes the cause column...among other issues
    sim$firePoints <- getFirePoints_NFDB_scfm(
      studyArea = sim$fireRegimePolysCalibration,
      NFDB_pointPath = checkPath(file.path(dPath, "NFDB_point"), create = TRUE)
    )
    #TODO: why is this necessary?
    sim$firePoints <- postProcess(sim$firePoints, studyArea = sim$fireRegimePolysCalibration)
  }


  return(invisible(sim))
}
