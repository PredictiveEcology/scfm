defineModule(sim, list(
  name = "scfmRegime",
  description = "estimates fire regime parameters for BEACONs a la Steve's method",
  keywords = c("fire regime", "BEACONs"),
  authors = c(
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(),
  version = numeric_version("0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.md", "scfmRegime.Rmd"), ## same file
  loadOrder = list(after = c("scfmLandcoverInit"),
                   before = c("scfmDriver", "scfmIgnition", "scfmEscape", "scfmSpread")),
  reqdPkgs = list("raster", "reproducible", "PredictiveEcology/scfmutils (>= 0.0.13.9001)"),
  parameters = rbind(
    defineParameter("empiricalMaxSizeFactor", "numeric", 1.2, 1, 10, "scale xMax by this is HD estimator fails "),
    defineParameter("fireCause", "character", c("N", "L"), NA_character_, NA_character_,
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
    defineParameter("targetBurnRate", "numeric", NA, 0, 1,
                    desc = paste("a named vector giving the proportional annual area burned of each fire regime polygon.",
                                 "These override the default estimate of scfm and are used to estimate a new mean",
                                 "fire size and ignition rate. Names should correspond to `PolyID`.",
                                 "A partial set of polygons is allowed - missing polys are estimated from data")),
    defineParameter("targetMaxFireSize", "numeric", NA, 0, NA,
                    desc = paste("a named vector giving the estimated max fire size in ha of each fire regime polygon.",
                                 "These will override the default estimate of scfm and will be used to estimate",
                                 "a new spread probability. Names should correspond to `PolyID`.",
                                 "A partial set of polygons is allowed - missing polys are estimated from data")),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name to be cached by SpaDES.")
  ),
  inputObjects = bindrows(
    expectsInput("firePoints", "SpatialPointsDataFrame",
                 desc = paste0("Historical fire data in point form. Must contain fields 'CAUSE',
                               'YEAR', and 'SIZE_HA', or pass the parameters to identify those."),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput("fireRegimePolys", "sf",
                 desc = paste("Areas to calibrate individual fire regime parameters. Defaults to ecoregions.",
                              "Must have numeric field 'PolyID' or it will be created for individual polygons")),
    expectsInput("fireRegimePolysLarge", "sf",
                 desc = paste("A polygons file with field 'PolyID' describing unique fire regimes in a larger",
                              "study area. Not required - but useful if the parameterization region is different",
                              "from the simulation region.")),
    expectsInput("landscapeAttr", "list", ## TODO: use sf object (#32)
                 desc = "list of landscape attributes for each polygon"),
    expectsInput("landscapeAttrLarge", "list", ## TODO: use sf object (#32)
                 desc = paste("list of landscape attributes for larger study area - if supplied, the module",
                              "will generate fire regime parameters for the polygons in `landscapeAttr`",
                              "using the attributes from `landscapeAttrLarge`.")),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = paste("template raster for raster GIS operations.",
                              "Must be supplied by user with same CRS as `studyArea`.")),
    expectsInput("rasterToMatchLarge", "RasterLayer",
                 desc = paste("large template raster for raster GIS operations.",
                              "Must be supplied by user with same CRS as `studyAreaLarge`.")),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = "Polygon to use as the simulation study area.",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the parametrisation study area.",
                              "Note that `studyAreaLarge` is only used for parameter estimation, and",
                              "can be larger than the actual study area used for simulations."),
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip")
  ),
  outputObjects = bindrows(
    createsOutput("fireRegimePoints", "SpatialPointsDataFrame",
                  desc = "Fire locations. Points outside studyArea are removed"),
    createsOutput("scfmRegimePars", "list", ## TODO: use sf object (#32)
                  desc =  "list of fire regime parameters for each polygon")
  )
))

doEvent.scfmRegime = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- Init(sim)
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

Init <- function(sim) {
  tmp <- sim$firePoints
  ## extract and validate fireCause spec

  fc <- P(sim)$fireCause

  ## review that sf can be used like this.
  ## should verify CAUSE is a column in the table...
  if (!P(sim)$fireCauseColumnName %in% names(tmp)) {
    stop("The column ", P(sim)$fireCauseColumnName, " does not exist in the fire database used. ",
         "Please pass the correct column name for the fire cause.")
  }
  if (is.factor(tmp[[P(sim)$fireCauseColumnName]])) {
    causeSet <- levels(tmp[[P(sim)$fireCauseColumnName]])
  } else {
    causeSet <- unique(tmp[[P(sim)$fireCauseColumnName]])
  }

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
  if (length(epoch) != 2 || !is.numeric(epoch) || any(!is.finite(epoch)) || epoch[1] > epoch[2])
    stop("illegal fireEpoch: ", epoch)

  quotes <- paste0("tmp$", paste(eval(P(sim)$fireYearColumnName)))
  tmp <- subset(tmp, get(P(sim)$fireYearColumnName) >= epoch[1] &
                  get(P(sim)$fireYearColumnName) <= epoch[2])

  epochLength <- as.numeric(epoch[2] - epoch[1] + 1)

  if (sf::st_crs(tmp) != sf::st_crs(sim$fireRegimePolysLarge)) {
    tmp <- sf::st_transform(tmp, crs = sf::st_crs(sim$fireRegimePolysLarge))
  }

  tmp <- sf::st_intersection(tmp, sim$fireRegimePolysLarge) ## gives studyArea colnames to points

  if (any(is.na(tmp$PolyID))) {
    tmp <- tmp[!is.na(tmp$PolyID), ] ## need to remove NA points
  }
  sim$fireRegimePoints <- tmp

  ## this function estimates the ignition probability and escape probability based on NFDB
  scfmRegimePars <- lapply(names(sim$landscapeAttrLarge),
                           FUN = calcZonalRegimePars,
                           firePolys = sim$fireRegimePolysLarge,
                           landscapeAttr = sim$landscapeAttrLarge,
                           firePoints = sim$fireRegimePoints,
                           epochLength = epochLength,
                           maxSizeFactor = P(sim)$empiricalMaxSizeFactor,
                           fireSizeColumnName = P(sim)$fireSizeColumnName,
                           targetBurnRate = P(sim)$targetBurnRate,
                           targetMaxFireSize = P(sim)$targetMaxFireSize)

  names(scfmRegimePars) <- names(sim$landscapeAttrLarge)
  ## only keep the attributes that are in study area
  scfmRegimePars <- scfmRegimePars[names(scfmRegimePars) %in% names(sim$landscapeAttr)]

  nullIdx <- sapply(scfmRegimePars, is.null)
  if (any(nullIdx)) {
    scfmRegimePars <- scfmRegimePars[-which(nullIdx)]
  }
  sim$scfmRegimePars <- scfmRegimePars

  stopifnot(
    identical(unique(names(sim$scfmRegimePars)), names(sim$scfmRegimePars))
  )

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    sim$studyAreaLarge <- sim$studyArea
  }

  if (!suppliedElsewhere("rasterToMatchLarge", sim)) {
    sim$rasterToMatchLarge <- sim$rasterToMatch
  }

  if (!suppliedElsewhere("fireRegimePolys", sim)) {
    message("fireRegimePolys not supplied. Using default ", P(sim)$fireRegimePolysType, " of Canada.")

    sim$fireRegimePolys <- Cache(
      scfmutils::prepInputsFireRegimePolys,
      url = extractURL("fireRegimePolys", sim),
      destinationPath = dPath,
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      type = P(sim)$fireRegimePolysType,
      userTags = c(cacheTags, "fireRegimePolys")
    )
  }

  if (!suppliedElsewhere("fireRegimePolysLarge", sim) & !is.null(sim$studyAreaLarge)) {
    message("fireRegimePolys not supplied. Using default ", P(sim)$fireRegimePolysType, " of Canada.")

    sim$fireRegimePolysLarge <- Cache(
      scfmutils::prepInputsFireRegimePolys,
      url = NULL,
      destinationPath = dPath,
      studyArea = sim$studyAreaLarge,
      rasterToMatch = sim$rasterToMatchLarge,
      type = P(sim)$fireRegimePolysType,
      userTags = c(cacheTags, "fireRegimePolys")
    )
  } else {
    sim$fireRegimePolysLarge <- sim$fireRegimePolys
  }

  if (!suppliedElsewhere("firePoints", sim)) {
    ## NOTE: do not use fireSenseUtils - it removes the cause column...among other issues
    sim$firePoints <- getFirePoints_NFDB_scfm(
      studyArea = sim$fireRegimePolysLarge,
      NFDB_pointPath = checkPath(file.path(dPath, "NFDB_point"), create = TRUE)
    )
    sim$firePoints <- postProcess(sim$firePoints, studyArea = sim$fireRegimePolysLarge)
  }

  return(invisible(sim))
}
