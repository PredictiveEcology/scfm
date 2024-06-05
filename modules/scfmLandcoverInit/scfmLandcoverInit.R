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
  timeframe = as.POSIXlt(c("2005-01-01", NA)),
  documentation = list("README.md", "scfmLandcoverInit.Rmd"), # same file
  loadOrder = list(after = c("Biomass_speciesData", "Biomass_borealDataPrep"),
                   before = c("scfmRegime", "scfmDriver", "scfmEscape", "scfmIgnition", "scfmSpread")),
  timeunit = "year",
  citation = list(),
  reqdPkgs = list(
    "fasterize", ## TODO: use terra::rasterize()
    "PredictiveEcology/LandR@development",
    "purrr",
    "reproducible",
    "PredictiveEcology/scfmutils (>= 0.0.7.9001)",
    "sf", "terra"
  ),
  parameters = rbind(
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of immediate cell neighbours"),
    defineParameter("sliverThreshold", "numeric", 1e8, NA, NA,
                    paste("fire regime polygons with area less than this number will be merged",
                          "with their closest non-sliver neighbour using `sf::st_nearest_feature`.")),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "Initial time for plotting"),
    defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, "Interval between plotting"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by `Plots` function, which can be optionally used here."),
    defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, "Initial time for saving"),
    defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, "Interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Use caching of events - not recommended as of 10/05/2023")
  ),
  inputObjects = bindrows(
    expectsInput(
      "fireRegimePolys", "sf",
      desc = paste(
        "Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
        "Must have numeric field 'PolyID' or it will be created for individual polygons"
      )
    ),
    expectsInput(
      "fireRegimePolysLarge", "sf",
      desc = paste(
        "if `studyAreaLarge` is supplied, the corresponding fire regime areas.",
        "Must have integer field `PolyID` if supplied, and uses same defaults as `fireRegimePolys`."
      )
    ),
    expectsInput(
      "flammableMap", "SpatRaster",
      desc = "binary flammability map - defaults to using LandR::prepInputsLCC"
    ),
    expectsInput(
      "flammableMapLarge", "SpatRaster",
      desc = paste(
        "binary flammability map - defaults to using `LandR::prepInputsLCC`.",
        "This is only necessary if passing `studyAreaLarge` OR running `scfmDriver`.",
        "It should match the extent of `studyAreaLarge`, and if running `scfmDriver`,",
        "it should extend by at least `scfmDriver`'s `buffDist`."
      )
    ),
    expectsInput(
      "rasterToMatch", "SpatRaster",
      desc = "template raster for raster GIS operations. Must be supplied by user"
    ),
    expectsInput(
      "rasterToMatchLarge", "SpatRaster",
      desc = paste(
        "Template raster for raster GIS operations. Only necessary if `studyAreaLarge` is passed.",
        "Must be supplied by user."
      )
    ),
    expectsInput(
      "studyArea", "SpatialPolygonsDataFrame", ## TODO: should be sf?
      desc = "Polygon to use as the simulation study area (typically buffered)."
    ),
    expectsInput(
      "studyAreaLarge", "SpatialPolygonsDataFrame", ## TODO: should be sf?
      desc = "optional larger study area used for parameterization but not simulation"
    )
  ),
  outputObjects = bindrows(
    createsOutput(
      "cellsByZone", "data.frame",
      desc = "explains which raster cells are in which polygon"
    ),
    createsOutput(
      "landscapeAttr", "list", ## TODO: use sf object (#32)
      desc = "list of polygon attributes inc. area"),
    createsOutput(
      "landscapeAttrLarge", "list", ## TODO: use sf object (#32)
      desc = paste(
        "if `studyAreaLarge` is passed, this object will supersede `landscapeAttr` in `scfmRegime`,",
        "so that estimates of mean fire size, max fire size, ignition prob, and escape prob",
        "are based on `fireRegimePolysLarge`. Allows for calibration over larger area."
      )
    ),
    createsOutput(
      "fireRegimePolys", "SpatialPolygonsDataFrame",
      desc = paste(
        "areas to calibrate individual fire regime parameters. If supplied, it must",
        "have a field called `PolyID` that defines unique regimes. Defaults to ecozones."
      )
    ),
    createsOutput(
      "fireRegimePolysLarge", "SpatialPolygonsDataFrame",
      desc = paste(
        "areas to calibrate individual fire regime parameters if `studyAreaLarge` is passed.",
        "If supplied, it must have a field `PolyID` used to define unique fire regimes."
      )
    ),
    createsOutput(
      "fireRegimeRas", "SpatRaster",
      desc = "Rasterized version of `fireRegimePolys` with values representing polygon ID."
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

  ## ensure flammability maps are integer ('binary') maps
  stopifnot(
    all(unique(sim$flammableMap[]) %in% c(NA_integer_, 0L, 1L)),
    all(unique(sim$flammableMapLarge[]) %in% c(NA_integer_, 0L, 1L))
  )

  if (!is.integer(sim$flammableMap[])) {
    sim$flammableMap[] <- as.integer(sim$flammableMap[])
  }

  if (!is.integer(sim$flammableMapLarge[])) {
    sim$flammableMapLarge <- setValues(sim$flammableMapLarge,
                                       as.integer(values(sim$flammableMapLarge,
                                                         mat = FALSE)))
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
    sim$fireRegimePolys <- postProcessTerra(sim$fireRegimePolysLarge, studyArea = sim$studyArea)
    #for now - GIS operations with sf objects are causing sliver polygons (area < 0.001 m2)

    if (is(st_geometry(sim$fireRegimePolys), "sfc_GEOMETRY")) {
      #this object may have empty geometries, which can occur when SAL and SA are both subsets of the same file
      #the empty geometries will cause an error
      sim$fireRegimePolys <- sim$fireRegimePolys[as.numeric(st_area(sim$fireRegimePolys)) > 0,]
      sim$fireRegimePolys <- st_cast(sim$fireRegimePolys, "MULTIPOLYGON")
    }

    # This makes sim$landscapeAttr & sim$cellsByZone
    sim$fireRegimePolysLarge <- sim$fireRegimePolysLarge[order(sim$fireRegimePolysLarge$PolyID),]

    sim$landscapeAttrLarge <- Cache(genFireMapAttr,
      flammableMap = sim$flammableMapLarge,
      fireRegimePolys = sim$fireRegimePolysLarge,
      neighbours = P(sim)$neighbours,
      userTags = c(currentModule(sim), "genFireMapAttr", "studyAreaLarge")
    )
  }

  sim$fireRegimePolys <- checkForIssues(
    fireRegimePolys = sim$fireRegimePolys,
    studyArea = sim$studyArea,
    rasterToMatch = sim$rasterToMatch,
    flammableMap = sim$flammableMap,
    sliverThresh = P(sim)$sliverThreshold,
    cacheTag = c("scfmLandcoverInit", "fireRegimePolys")
  )
  sim$fireRegimePolys <- sim$fireRegimePolys[order(sim$fireRegimePolys$PolyID),]

  sim$landscapeAttr <- Cache(genFireMapAttr,
    flammableMap = sim$flammableMap,
    fireRegimePolys = sim$fireRegimePolys,
    neighbours = P(sim)$neighbours,
    userTags = c(currentModule(sim), "genFireMapAttr", "studyArea")
  )

  ##############
  # ONLY FOR SA

  # doing this prevents fireRegimeRas from inheriting colormaps
  sim$fireRegimeRas <- rasterize(sim$fireRegimePolys, sim$rasterToMatch, fun = "max", field = "PolyID")

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)

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

  if (hasSAL & !suppliedElsewhere("rasterToMatchLarge", sim)) {
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

  if (!suppliedElsewhere("flammableMapLarge", sim) & hasSAL) {
    if (!is.null(sim$flammableMap)) {
      stop("flammableMap was supplied but not flammableMapLarge. Please supply neither or both")
    }

    vegMap <- prepInputsLCC(
      year = 2010,
      destinationPath = dPath,
      studyArea = sim$studyAreaLarge,
      rasterToMatch = sim$rasterToMatchLarge,
      userTags = c("prepInputsLCC", "studyAreaLarge")
    )
    if (!is.integer(vegMap[])) {
      vegMap <- setValues(vegMap, as.integer(values(vegMap)))
    }

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

    fireRegimePolys <- Cache(prepInputsFireRegimePolys, url = NULL, destinationPath = dPath,
                             studyArea = sa, type = "ECOREGION") %>%
      st_transform(., st_crs(sa))

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
