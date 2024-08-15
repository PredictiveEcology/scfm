defineModule(sim, list(
  name = "scfmLandcoverInit",
  description = paste(
    "Generates some relevant statistics for each fire regime over a `studyArea`.",
    "If scfm is being parameterized over a larger area (`studyAreaCalibration`), then the",
    "following objects must be supplied with identical CRS and resolution, where applicable:",
    "`studyArea`, `studyAreaCalibration`, `rasterToMatch`, `rasterToMatchCalibration`.",
    "The extent should differ between objects and their 'large' counterparts."
  ),
  keywords = c("fire", "land cover classification"),
  childModules = character(),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  version = list(scfmLandcoverInit = "2.0.0"),
  timeframe = as.POSIXlt(c("2005-01-01", NA)),
  documentation = list("README.md", "scfmLandcoverInit.Rmd"), # same file
  loadOrder = list(after = c("Biomass_speciesData", "Biomass_borealDataPrep"),
                   before = c("scfmRegime", "scfmDriver", "scfmEscape", "scfmIgnition", "scfmSpread")),
  timeunit = "year",
  citation = list(),
  reqdPkgs = list(
    "PredictiveEcology/LandR (>= 1.1.1)",
    "purrr",
    "PredictiveEcology/scfmutils (>= 2.0.1)",
    "reproducible", "sf", "terra"
  ),
  parameters = rbind(
    defineParameter("dataYear", "numeric", 2011, 1985, 2020,
                    desc = paste("used to select the year of landcover data used to create",
                                 "flammableMap if the obejct is unsupplied")),
    defineParameter("fireRegimePolysType", "character", "ECOREGION", NA, NA,
                    paste("Polygon type to use for scfm `fireRegimePolys`:",
                          "see `?scfmutils::prepInputsFireRegimePolys` for allowed types.")),
    defineParameter("neighbours", "numeric", 8, NA, NA, "Number of immediate cell neighbours"),
    defineParameter("sliverThreshold", "numeric", 6.25e8, NA, NA,
                    paste("fire regime polygons with area (in m2) less than this number will be merged",
                          "with their closest non-sliver neighbour using `sf::st_nearest_feature`.")),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "Initial time for plotting"),
    defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, "Interval between plotting"),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by `Plots` function, which can be optionally used here."),
    defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, "Initial time for saving"),
    defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, "Interval between save events"),
    defineParameter(".useCache", "character", ".inputObjects", NA, NA,
                    "Use caching of events - not recommended as of 10/05/2023")
  ),
  inputObjects = bindrows(
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
                 desc = paste("binary flammability map - defaults to using `LandR::prepInputsLCC`.",
                              "This is only necessary if passing `studyAreaCalibration` OR running `scfmDriver`.",
                              "It should match the extent of `studyAreaCalibration`, and if running `scfmDriver`,",
                              "it should extend by >= scfmDriver's `P(sim)$buffDist`.")),
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = "template raster for raster GIS operations. Must be supplied by user"),
    expectsInput("rasterToMatchCalibration", "SpatRaster",
                 desc = paste("Template raster for raster GIS operations. Only necessary if `studyAreaCalibration` is passed.",
                              "Must be supplied by user.")),
    expectsInput("studyArea", "sf", desc = "Polygon to use as the simulation study area (typically buffered)."),
    expectsInput("studyAreaCalibration", "sf", desc = "optional larger study area used for parameterization only")
  ),
  outputObjects = bindrows(
    createsOutput("fireRegimePolys", "sf",
                  desc = "`fireRegimePolys` with landcover attributes appended"),
    createsOutput("fireRegimePolysCalibration", "sf",
                  desc = "`fireRegimePolysCalibration` with landcover attributes appended"),
    createsOutput("fireRegimeRas", "SpatRaster",
                  desc = "Rasterized version of fireRegimePolys with values representing polygon ID")
  )
))

doEvent.scfmLandcoverInit <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
    init = {
      sim <- Init(sim)

      if (anyPlotting(P(sim)$.plots)) {
        sim <- scheduleEvent(sim, start(sim), "scfmLandcoverInit", "plot")
      }

      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "scfmLandcoverInit", "save")
    },
    plot = {
      ## NOTE: these objects don't change during sim, so only need to be plotted once
      Plots(sim$fireRegimeRas, fn = scfmutils::plot_fireRegimeRas, type = P(sim)$.plots,
            filename = paste0("fireRegimeRas"),
            title = paste0("Fire regimes"))

      Plots(sim$flammableMap, fn = scfmutils::plot_flammableMap, type = P(sim)$.plots,
            filename = paste0("flammableMap"),
            title = paste0("landscape flammability map"))
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
  # doing so on larger and smaller has the potential to mismatch slivers between calibration/simulation
  if (!is.null(sim$fireRegimePolysCalibration)) {
    if (is(sim$fireRegimePolysCalibration, "SpatialPolygonsDataFrame")) {
      sim$fireRegimePolysCalibration <- sf::st_as_sf(sim$fireRegimePolysCalibration)
    }

    sim$fireRegimePolysCalibration <- checkForIssues(
      fireRegimePolys = sim$fireRegimePolysCalibration,
      studyArea = sim$studyAreaCalibration,
      rasterToMatch = sim$rasterToMatchCalibration,
      flammableMap = sim$flammableMapCalibration,
      sliverThresh = P(sim)$sliverThreshold,
      cacheTag = c("scfmLandcoverInit", "fireRegimePolysCalibration")
    )

    ## now that slivers are removed, remake frp from the larger object
    sim$fireRegimePolys <- postProcessTerra(sim$fireRegimePolysCalibration, studyArea = sim$studyArea)
    ## for now - GIS operations with sf objects are causing sliver polygons (area < 0.001 m2)

    if (is(st_geometry(sim$fireRegimePolys), "sfc_GEOMETRY")) {
      ## this object may have empty geometries, which can occur when SAC and SA are both subsets
      ## of the same file. the empty geometries will cause an error.
      sim$fireRegimePolys <- sim$fireRegimePolys[as.numeric(st_area(sim$fireRegimePolys)) > 0, ]
      sim$fireRegimePolys <- st_cast(sim$fireRegimePolys, "MULTIPOLYGON")
    }

    sim$fireRegimePolysCalibration <- sim$fireRegimePolysCalibration[order(sim$fireRegimePolysCalibration$PolyID), ]

    sim$fireRegimePolysCalibration <- Cache(genFireMapAttr,
      flammableMap = sim$flammableMapCalibration,
      fireRegimePolys = sim$fireRegimePolysCalibration,
      neighbours = P(sim)$neighbours,
      userTags = c(currentModule(sim), "genFireMapAttr", "studyAreaCalibration")
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

  sim$fireRegimePolys <- Cache(genFireMapAttr,
    flammableMap = sim$flammableMap,
    fireRegimePolys = sim$fireRegimePolys,
    neighbours = P(sim)$neighbours,
    userTags = c(currentModule(sim), "genFireMapAttr", "studyArea")
  )

  ## doing this prevents fireRegimeRas from inheriting colormaps
  sim$fireRegimeRas <- rasterize(sim$fireRegimePolys, sim$rasterToMatch, fun = "max", field = "PolyID")

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

  # supply objects
  if (!hasSA & !hasSAC) {
    message("study area not supplied. Using random polygon in Alberta")
    studyArea <- LandR::randomStudyArea(size = 15000000000, seed = 23654)
    sim$studyArea <- studyArea
    sim$studyAreaCalibration <- studyArea
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message(paste(
      "rasterToMatch not supplied. generating from LCC2010 using studyArea CRS",
      " - It is strongly recommended to supply a rasterToMatch"
    ))
    sim$rasterToMatch <- LandR::prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      cropTo = sim$studyArea,
      maskTo= sim$studyArea,
      projectTo = sim$studyArea,
      filename2 = NULL,
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatch")
    )
  }

  if (hasSAC & !suppliedElsewhere("rasterToMatchCalibration", sim)) {
    message(paste(
      "rasterToMatch not supplied. generating from NTEMS LCC using studyArea CRS",
      " - It is strongly recommended to supply a rasterToMatch"
    ))
    sim$rasterToMatchCalibration <- LandR::prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      maskTo = sim$studyAreaCalibration,
      projectTo = sim$studyAreaCalibration,
      cropTo = sim$studyAreaCalibration,
      filename2 = NULL,
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatchCalibration")
    )
  }

  if (!suppliedElsewhere("flammableMapCalibration", sim) & hasSAC) {
    if (!is.null(sim$flammableMap)) {
      stop("flammableMap was supplied but not flammableMapCalibration. Please supply neither or both")
    }

    vegMap <- prepInputs_NTEMS_LCC_FAO(
      year = P(sim)$dataYear,
      destinationPath = dPath,
      maskTo = sim$studyAreaCalibration,
      cropTo = sim$rasterToMatchCalibration,
      projectTo = sim$rasterToMatchCalibration,
      userTags = c("prepInputs_NTEMS_LCC_FAO", "studyArea")
    )
    vegMap[] <- asInteger(vegMap[])
    sim$flammableMapCalibration <- defineFlammable(vegMap,
                                             mask = sim$rasterToMatchCalibration,
                                             nonFlammClasses = c(20, 31, 32, 33)
    )
  }

  if (!suppliedElsewhere("flammableMap", sim)) {
    if (hasSAC) {
      useTerra <- getOption("reproducible.useTerra") ## TODO: reproducible#242
      options(reproducible.useTerra = FALSE) ## TODO: reproducible#242
      sim$flammableMap <- postProcess(sim$flammableMapCalibration, rasterToMatch = sim$rasterToMatch)
      options(reproducible.useTerra = useTerra) ## TODO: reproducible#242
    } else {
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
  }

  ## this is TRUE unless fireRegimePolysCalibration is supplied, in which case we drop that object
  if (!hasFRP & !hasFRPC) {
    sa <- if (hasSAC) {
      sim$studyAreaCalibration
    } else {
      sim$studyArea
    }
    message("fireRegimePolys not supplied. Using default ecoregions of Canada")
    # cannot use prepInputs with a vector for prepInputs - unreliable w/ GDAL

    fireRegimePolys <- Cache(prepInputsFireRegimePolys, url = NULL, destinationPath = dPath,
                             studyArea = sa, type = P(sim)$fireRegimePolysType) %>%
      st_transform(., st_crs(sa))

    if (hasSAC) {
      sim$fireRegimePolysCalibration <- fireRegimePolys
      sim$fireRegimePolys <- postProcess(fireRegimePolys,
                                         studyArea = sim$studyArea
      )
    } else {
      sim$fireRegimePolys <- fireRegimePolys
    }
  }

  return(invisible(sim))
}
