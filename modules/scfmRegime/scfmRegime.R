defineModule(sim, list(
  name = "scfmRegime",
  description = "estimates fire regime parameters for BEACONs a la Steve's method",
  keywords = c("fire regime", "BEACONs"),
  authors = c(person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
              person("Ian", "Eddy", email = "ian.eddy@canada.ca", role = c("aut"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmRegime.Rmd"),
  reqdPkgs = list("raster", "reproducible", "sp"),
  parameters = rbind(
    defineParameter("empiricalMaxSizeFactor", "numeric", 1.2, 1, 10, "scale xMax by this is HD estimator fails "),
    defineParameter("fireCause", "character", c("L"), NA_character_, NA_character_, "subset of c(H,H-PB,L,Re,U)"),
    defineParameter("fireCauseColumnName", "character", "CAUSE", NA, NA,
                    desc = "Name of the column that has fire cause, consistent with P(sim)$fireCause"),
    defineParameter("fireEpoch", "numeric", c(1971,2000), NA, NA, "start of normal period"),
    defineParameter("fireSizeColumnName", "character", "SIZE_HA", NA, NA,
                    desc = "Name of the column that has fire size"),
    defineParameter("fireYearColumnName", "character", "YEAR", NA, NA,
                    desc = "Name of the column that has fire size"),
    defineParameter("targetBurnRate", "numeric", NA, 0, 1,
                    desc = paste("a named vector giving the proportional annual area burned of each fire regime polygon.",
                                 "These override the default estimate of scfm and are used to estimate a new mean",
                                 "fire size and ignition rate. Names should correspond to polyID.",
                                 "A partial set of polygons is allowed - missing polys are estimated from data")),
    defineParameter("targetMaxFireSize", "numeric", NA, 0, NA,
                    desc = paste("a named vector giving the estimated max fire size in ha of each fire regime polygon.",
                                 "These will override the default estimate of scfm and will be used to estimate",
                                 "a new spread probability. Names should correspond to polyID.",
                                 "A partial set of polygons is allowed - missing polys are estimated from data")),
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "firePoints", objectClass = "SpatialPointsDataFrame",
                 desc = paste0("Historical fire data in point form. Must contain fields 'CAUSE',
                               'YEAR', and 'SIZE_HA', or pass the parameters to identify those"),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list",
                 desc = "list of landscape attributes for each polygon"),
    expectsInput(objectName = "landscapeAttrLarge", objectClass = "list",
                 desc = paste("list of landscape attributes for larger study area - if supplied, the module",
                              "will generate fire regime parameters for the polygons in landscapeAttr",
                              "using the attributes from landscapeAttrLarge.")),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer",
                 desc = "template raster for raster GIS operations. Must be supplied by user with same CRS as studyArea"),
    expectsInput(objectName = "fireRegimePolys", objectClass = "SpatialPolygonsDataFrame",
                 desc = paste("Areas to calibrate individual fire regime parameters. Defaults to ecoregions.",
                              "Must have numeric field 'PolyID' or it will be created for individual polygons"),
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip"),
    expectsInput(objectName = "fireRegimePolysLarge", objectClass = "SpatialPolygonsDataFrame",
                 desc = paste("A polygons file with field 'PolyID' describing unique fire regimes in a larger",
                              "study area. Not required - but useful if the parameterization region is different",
                              "from the simulation region."))
  ),
  outputObjects = bindrows(
   createsOutput(objectName = "scfmRegimePars", objectClass = "list",
                 desc =  "list of fire regime parameters for each polygon"),
   createsOutput(objectName = "fireRegimePoints", objectClass = "SpatialPointsDataFrame",
                 desc = "Fire locations. Points outside studyArea are removed")
  )
))


## event types
#   - type `init` is required for initiliazation

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
  #extract and validate fireCause spec

  fc <- P(sim)$fireCause

  #review that sf can be used like this.
  #should verify CAUSE is a column in the table...
  if (!P(sim)$fireCauseColumnName %in% names(tmp))
    stop(paste0("The column ", P(sim)$fireCauseColumnName, " does not exist",
                " in the fire database used. Please pass the correct column name ",
                "for the fire cause."))
  if (is.factor(tmp[[P(sim)$fireCauseColumnName]])){
    causeSet <- levels(tmp[[P(sim)$fireCauseColumnName]])}
  else {
    causeSet <- unique(tmp[[P(sim)$fireCauseColumnName]])
  }
  if (any(!(fc %in% causeSet))) {
    notPresent <- fc[!fc %in% causeSet]
    warning(paste0("This firecause is not present: ", notPresent,
                   " The following are the fire causes: ",
                   paste(causeSet, collapse = ", "),
                   ". Original cause will be replaced by ",
                   paste(causeSet, collapse = ", ")), immediate. = TRUE)
    fc <- causeSet
  }
  tmp <- subset(tmp,  get(P(sim)$fireCauseColumnName) %in% fc)

  #extract and validate fireEpoch
  epoch <- P(sim)$fireEpoch
  if (length(epoch) != 2 ||
      !is.numeric(epoch) || any(!is.finite(epoch)) || epoch[1] > epoch[2])
    stop("illegal fireEpoch: ", epoch)

  quotes <- paste0("tmp$", paste(eval(P(sim)$fireYearColumnName)))
  tmp <- subset(tmp, get(P(sim)$fireYearColumnName) >= epoch[1] &
                  get(P(sim)$fireYearColumnName) <= epoch[2])

  epochLength <- as.numeric(epoch[2] - epoch[1] + 1)

  if (!is.null(sim$fireRegimePolysLarge) & !is.null(sim$landscapeAttrLarge)) {
    #sp over uses identical CRS

    if (st_crs(tmp) != st_crs(sim$fireRegimePolysLarge)) {
      tmp <- st_transform(tmp, crs = st_crs(sim$fireRegimePolysLarge))
    }

    tmp <- sf::st_intersection(tmp, sim$fireRegimePolysLarge) #gives studyArea row name to point

    if (any(is.na(tmp$PolyID))) {
      tmp <- tmp[!is.na(tmp$PolyID),] #have to remove NA points
    }
    sim$fireRegimePoints <- tmp
    #this function estimates the ignition probability and escape probability based on NFDB
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
    #only keep the attribtues that are in study area
    scfmRegimePars <- scfmRegimePars[names(scfmRegimePars) %in% names(sim$landscapeAttr)]

  } else {

    tmp$PolyID <- sp::over(tmp, sim$fireRegimePolys)$PolyID #gives studyArea row name to point

    if (any(is.na(tmp$PolyID))) {
      tmp <- tmp[!is.na(tmp$PolyID),] #have to remove NA points
    }
    sim$fireRegimePoints <- tmp

    #this function estimates the ignition probability and escape probability based on NFDB
    scfmRegimePars <- lapply(names(sim$landscapeAttr),
                             FUN = calcZonalRegimePars,
                             firePolys = sim$fireRegimePolys,
                             landscapeAttr = sim$landscapeAttr,
                             firePoints = sim$fireRegimePoints,
                             epochLength = epochLength,
                             maxSizeFactor = P(sim)$empiricalMaxSizeFactor,
                             fireSizeColumnName = P(sim)$fireSizeColumnName,
                             targetBurnRate = P(sim)$targetBurnRate,
                             targetMaxFireSize = P(sim)$targetMaxFireSize)

    names(scfmRegimePars) <- names(sim$landscapeAttr)
  }

  nullIdx <- sapply(scfmRegimePars, is.null)
  if (any(nullIdx)){
    scfmRegimePars <- scfmRegimePars[-which(nullIdx)]
  }
  sim$scfmRegimePars <- scfmRegimePars

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("fireRegimePolys", sim)) {
    message("fireRegimePolys not supplied. Using default ecoregions of Canada")

    sim$fireRegimePolys <- prepInputs(url = extractURL("fireRegimePolys", sim),
                                      destinationPath = dPath,
                                      studyArea = sim$studyArea,
                                      rasterToMatch = sim$rasterToMatch,
                                      overwrite = TRUE,
                                      userTags = c("cacheTags", "fireRegimePolys"))
    sim$fireRegimePolys$PolyID <- as.numeric(sim$fireRegimePolys$ECOREGION)
  }
  ## this module has many dependencies that aren't sourced in .inputObjects
  ## this workaround prevents checksums updating due to daily name change of NFDB files
  if (!suppliedElsewhere("firePoints", sim)) {
    if (!is.null(sim$firePolysLarge)) {
      SA <- sim$firePolysRegimeLarge
      RTM <- sim$rasterToMatchLarge
      } else {
      SA <- sim$fireRegimePolys
      RTM <- sim$rasterToMatch
      }

    #do not use fireSenseUtils - it removes the cause column...among other issues
    #this function came first - fireSenseUtils copied the name!
    sim$firePoints <- getFirePoints_NFDB_scfm(studyArea = SA,
                                              rasterToMatch = RTM,
                                              NFDB_pointPath = checkPath(file.path(dataPath(sim), "NFDB_point"),
                                                                         create = TRUE))
  }

  return(invisible(sim))
}
