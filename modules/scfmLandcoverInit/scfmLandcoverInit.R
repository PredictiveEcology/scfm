defineModule(sim, list(
  name = "scfmLandcoverInit",
  description = paste(
    "Generates some relevant statistics for each fire regime over a studyArea.",
    "if scfm is being parameterized over a larger area (`studyAreaLarge`), then the",
    "following objects must be supplied with identical CRS and resolution, where applicable:",
    "`studyArea`, `studyAreaLarge`, `rasterToMatch`, `rasterToMatchLarge.`",
    "The extent should differ between objects and their 'large' counterparts."
  ),
  keywords = c("fire", "LCC2010", "land cover classification 2010", "BEACONs"),
  childModules = character(),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c("2005-01-01", NA)),
  documentation = list("README.md", "scfmLandcoverInit.Rmd"), # same file
  timeunit = "year",
  citation = list(),
  reqdPkgs = list(
    "fasterize", "purrr", "raster", "sf",
    "PredictiveEcology/LandR@development",
    "PredictiveEcology/reproducible@development",
    "PredictiveEcology/scfmutils (>= 0.0.0.9005)"
  ),
  parameters = rbind(
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of immediate cell neighbours"),
    defineParameter("sliverThreshold", "numeric", 1e8, NA, NA,
                    paste("fire regime polygons with area less than this number will be merged",
                          "with their closest non-sliver neighbour using sf::st_nearest_feature.")),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "Initial time for plotting"),
    defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, "Interval between plotting"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, "Initial time for saving"),
    defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, "Interval between save events"),
    defineParameter(".useCache", "logical", TRUE, NA, NA, "Use cache")
  ),
  inputObjects = bindrows(
    expectsInput(
      "studyArea", "SpatialPolygonsDataFrame",
      desc = "", ## TODO
      sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"
    ),
    expectsInput(
      "studyAreaLarge", "SpatialPolygonsDataFrame",
      desc = "optional larger study area used for parameterization but not simulation"
    ),
    expectsInput(
      "flammableMap", "RasterLayer",
      desc = "binary flammability map - defaults to using LandR::prepInputsLCC"
    ),
    expectsInput(
      "flammableMapLarge", "RasterLayer",
      desc = paste(
        "binary flammability map - defaults to using LandR::prepInputsLCC.",
        "This is only necessary if passing studyAreaLarge OR running scfmDriver.",
        "It should match the extent of studyAreaLarge, and if running scfmDriver,",
        "it should extend by >= scfmDriver's P(sim)$buffDist."
      )
    ),
    expectsInput(
      "rasterToMatch", "RasterLayer",
      desc = "template raster for raster GIS operations. Must be supplied by user"
    ),
    expectsInput(
      "rasterToMatchLarge", "RasterLayer",
      desc = paste(
        "Template raster for raster GIS operations. Only necessary if SAL is passed.",
        "Must be supplied by user."
      )
    ),
    expectsInput(
      "fireRegimePolys", "sf",
      desc = paste(
        "Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
        "Must have numeric field 'PolyID' or it will be created for individual polygons"
      ),
      sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip"
    ),
    expectsInput(
      "fireRegimePolysLarge", "sf",
      desc = paste(
        "if StudyAreaLarge is supplied, the corresponding fire regime areas. Must have",
        "numeric field 'PolyID' if supplied, and uses same defaults as fireRegimePolys"
      ),
      sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip"
    )
  ),
  outputObjects = bindrows(
    createsOutput(
      "cellsByZone", "data.frame",
      desc = "explains which raster cells are in which polygon"
    ),
    createsOutput(
      "landscapeAttr", "list",
      desc = "list of polygon attributes inc. area"),
    createsOutput(
      "landscapeAttrLarge", "list",
      desc = paste(
        "if SAL is passed, this object will supersede landscapeAttr in scfmRegmie, so that",
        "estimates of mean fire size, max fire size, ignition prob, and escape prob",
        "are based on fireRegimePolysLarge. Allows for calibration over larger area."
      )
    ),
    createsOutput(
      "fireRegimePolys", "SpatialPolygonsDataFrame",
      desc = paste(
        "areas to calibrate individual fire regime parameters. If supplied, it must",
        "have a field called PolyID that defines unique regimes. Defaults to ecozones"
      )
    ),
    createsOutput(
      "fireRegimePolysLarge", "SpatialPolygonsDataFrame",
      desc = paste(
        "areas to calibrate individual fire regime parameters if studyAreaLarge is passed.",
        "If supplied, it MUST have a field PolyID used to define unique fire regimes"
      )
    ),
    createsOutput(
      "fireRegimeRas", "RasterLayer",
      desc = "Rasterized version of fireRegimePolys with values representing polygon ID"
    )
  )
))

doEvent.scfmLandcoverInit <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
    init = {
      sim <- Init(sim)

      if ("screen" %in% P(sim)$.plots) {
        sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmLandcoverInit", "plot")
      }

      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "scfmLandcoverInit", "save")
    },
    plot = {
      if (time(sim) > start(sim)) {
        Plot(sim$fireRegimeRas, title = c("fire regimes"), new = TRUE)
        Plot(sim$flammableMap, legend = FALSE)
      }
      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmLandcoverInit", "plot")
    },
    save = {
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "scfmLandcoverInit", "save")
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
      "' in module '", events(sim)[1, "moduleName", with = FALSE], "'",
      sep = ""
    ))
  )

  return(invisible(sim))
}

Init <- function(sim) {
  if (is(sim$fireRegimePolys, "SpatialPolygonsDataFrame")) {
    sim$fireRegimePolys <- sf::st_as_sf(sim$fireRegimePolys)
  }

  message("checking sim$fireRegimePolys for sliver polygons...")

  # this only needs to be done on the larger area, if it is provided
  # doing so on larger and smaller has the potential to mismatch slivers between calibration/simulation
  if (!is.null(sim$fireRegimePolysLarge)) {
    if (is(sim$fireRegimePolysLarge, "SpatialPolygonsDataFrame")) {
      sim$fireRegimePolysLarge <- sf::st_as_sf(sim$fireRegimePolysLarge)
    }

    sim$fireRegimePolysLarge <- checkForIssues(
      fireRegimePolys = sim$fireRegimePolysLarge,
      studyArea = sim$studyAreaLarge,
      rasterToMatch = sim$rasterToMatchLarge,
      flammableMap = sim$flammableMapLarge,
      sliverThresh = P(sim)$sliverThreshold,
      cacheTag = c("scfmLandcoverInit", "fireRegimePolysLarge")
    )

    # now that slivers are removed, remake frp from the larger object
    useTerra <- getOption("reproducible.useTerra") ## TODO: reproducible#242
    options(reproducible.useTerra = FALSE) ## TODO: reproducible#242
    sim$fireRegimePolys <- postProcess(sim$fireRegimePolysLarge, studyArea = sim$studyArea)
    options(reproducible.useTerra = useTerra) ## TODO: reproducible#242

    # remnant slivers will be whole in the larger object
    if (is(sim$fireRegimePolys$geometry, "sfc_GEOMETRY")) {
      sim$fireRegimePolys <- st_cast(sim$fireRegimePolys, "MULTIPOLYGON")
    }

    # This makes sim$landscapeAttr & sim$cellsByZone
    sim$landscapeAttrLarge <- Cache(genFireMapAttr,
      flammableMap = sim$flammableMapLarge,
      fireRegimePolys = sim$fireRegimePolysLarge,
      neighbours = P(sim)$neighbours,
      userTags = c(currentModule(sim), "genFireMapAttr", "studyAreaLarge")
    )
  } else {
    sim$fireRegimePolys <- checkForIssues(
      fireRegimePolys = sim$fireRegimePolys,
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      flammableMap = sim$flammableMap,
      sliverThresh = P(sim)$sliverThreshold,
      cacheTag = c("scfmLandcoverInit", "fireRegimePolys")
    )
  }
  sim$landscapeAttr <- Cache(genFireMapAttr,
    flammableMap = sim$flammableMap,
    fireRegimePolys = sim$fireRegimePolys,
    neighbours = P(sim)$neighbours,
    userTags = c(currentModule(sim), "genFireMapAttr", "studyArea")
  )

  ##############
  # ONLY FOR SA

  # fireRegimeRas is handy for post-simulation analyses
  fireRegimeRas <- fasterize::fasterize(sim$fireRegimePolys, sim$rasterToMatch, field = "PolyID")

  # doing this prevents fireRegimeRas from inheriting colormaps
  sim$fireRegimeRas <- raster(fireRegimeRas)
  sim$fireRegimeRas <- setValues(sim$fireRegimeRas, getValues(fireRegimeRas))

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  # object check for SA/FRP/FRPL/SAL - better to be strict with stops
  hasSA <- suppliedElsewhere("studyArea", sim)
  hasSAL <- suppliedElsewhere("studyAreaLarge", sim)
  hasFRP <- suppliedElsewhere("fireRegimePolys", sim)
  hasFRPL <- suppliedElsewhere("fireRegimePolysLarge", sim)

  # supply objects
  if (!hasSA & !hasSAL) {
    message("study area not supplied. Using random polygon in Alberta")
    studyArea <- LandR::randomStudyArea(size = 15000000000, seed = 23654)
    sim$studyArea <- studyArea
    sim$studyAreaLarge <- studyArea
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message(paste(
      "rasterToMatch not supplied. generating from LCC2010 using studyArea CRS",
      " - It is strongly recommended to supply a rasterToMatch"
    ))
    sim$rasterToMatch <- LandR::prepInputsLCC(
      year = 2010,
      destinationPath = dPath,
      studyArea = sim$studyArea,
      useSAcrs = TRUE,
      filename2 = NULL,
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatch")
    )
  }

  if (hasSAL & !suppliedElsewhere(sim$rasterToMatchLarge, sim)) {
    message(paste(
      "rasterToMatch not supplied. generating from LCC2010 using studyArea CRS",
      " - It is strongly recommended to supply a rasterToMatch"
    ))
    sim$rasterToMatchLarge <- LandR::prepInputsLCC(
      year = 2010,
      destinationPath = dPath,
      studyArea = sim$studyAreaLarge,
      useSAcrs = TRUE,
      filename2 = NULL,
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatchLarge")
    )
  }

  if (!suppliedElsewhere("flammableMapLarge") & hasSAL) {
    vegMap <- prepInputsLCC(
      year = 2010,
      destinationPath = dPath,
      studyArea = sim$studyAreaLarge,
      rasterToMatch = sim$rasterToMatchLarge,
      userTags = c("prepInputsLCC", "studyAreaLarge")
    )
    vegMap[] <- asInteger(vegMap[])
    sim$flammableMapLarge <- defineFlammable(vegMap,
      mask = sim$rasterToMatchLarge,
      nonFlammClasses = c(13L, 16L:19L)
    )
  }

  if (!suppliedElsewhere("flammableMap", sim)) {
    if (hasSAL) {
      useTerra <- getOption("reproducible.useTerra") ## TODO: reproducible#242
      options(reproducible.useTerra = FALSE) ## TODO: reproducible#242
      sim$flammableMap <- postProcess(sim$flammableMapLarge, rasterToMatch = sim$rasterToMatch)
      options(reproducible.useTerra = useTerra) ## TODO: reproducible#242
    } else {
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
  }

  ## this is TRUE unless fireRegimePolysLarge is supplied, in which case we drop that object
  if (!hasFRP & !hasFRPL) {
    sa <- if (hasSAL) {
      sim$studyAreaLarge
    } else {
      sim$studyArea
    }
    message("fireRegimePolys not supplied. Using default ecoregions of Canada")
    # cannot use prepInputs with a vector for prepInputs - unreliable w/ GDAL
    fireRegimePolys <- prepInputs(
      url = extractURL("fireRegimePolys", sim),
      destinationPath = dPath,
      studyArea = sa,
      fun = st_read,
      useSAcrs = TRUE,
      filename2 = NULL,
      userTags = c(cacheTags, "fireRegimePolys")
    )
    fireRegimePolys <- st_transform(fireRegimePolys, crs = crs(sim$rasterToMatch))

    # this should preserve ecoregions in row-names -
    fireRegimePolys$PolyID <- fireRegimePolys$ECOREGION
    if (hasSAL) {
      sim$fireRegimePolysLarge <- fireRegimePolys
      sim$fireRegimePolys <- postProcess(fireRegimePolys,
        studyArea = sim$studyArea
      )
    } else {
      sim$fireRegimePolys <- fireRegimePolys
    }
  }

  return(invisible(sim))
}
