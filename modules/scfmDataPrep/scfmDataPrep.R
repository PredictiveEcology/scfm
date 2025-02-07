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
  authors = authors = c(
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
    "ggplot2", "parallel",
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
    defineParameter("fireCause", "character", c("N"), NA_character_, NA_character_,
                    desc = "subset of `c('H', 'H-PB', 'N', 'Re', 'U')`"),
    defineParameter("fireCauseColumnName", "character", "CAUSE", NA, NA,
                    desc = "Name of the column that has fire cause, consistent with `P(sim)$fireCause`."),
    defineParameter("fireEpoch", "numeric", c(1971, 2000), NA, NA, "start of normal period"),
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
  return(invisible(sim))
}


prepare_scfmlandcoverInit <- function(sim) {

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
    vegMap <- prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      maskTo = sim$studyAreaCalibration,
      cropTo = sim$rasterToMatchCalibration,
      #projectTo = sim$rasterToMatchCalibration, #should be done after defineFlammable
      userTags = c("prepInputs_NTEMS_LCC_FAO", "studyArea")
    )
    vegMap[] <- asInteger(vegMap[])
    fmc <- defineFlammable(vegMap,
                           nonFlammClasses = c(20, 31, 32, 33))

    sim$flammableMapCalibration <- postProcess(sim$flammableMapCalibration,
                                               to = sim$rasterToMatchCalibration,
                                               method = "mode")
  }

  if (!hasFM) {
    sim$flammableMap <- postProcess(sim$flammableMapCalibration, to = sim$rasterToMatch)
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
    sim$fireRegimePolys <- postProcess(terra::vect(fireRegimePolys), #avoid GIS issue with sf
                                       to = sim$studyArea) |>
      sf::st_as_sf()
  }

  return(invisible(sim))
}
