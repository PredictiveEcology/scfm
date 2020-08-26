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
  reqdPkgs = list("raster", "reproducible", "rgdal", "sp"),
  parameters = rbind(
    defineParameter("empiricalMaxSizeFactor", "numeric", 1.2, 1, 10, "scale xMax by this is HD estimator fails "),
    defineParameter("fireCause", "character", c("L"), NA_character_, NA_character_, "subset of c(H,H-PB,L,Re,U)"),
    defineParameter("fireEpoch", "numeric", c(1971, 2000), NA, NA, "start of normal period"),
    defineParameter("targetBurnRate", "numeric", NA, 0, 1, "proportion anual area burned"),
    defineParameter("fireCauseColumnName", "character", "CAUSE", NA, NA,
                    paste0("Name of the column that has fire cause. ",
                           "In the latest dataset, its FIRECAUS.")),
    defineParameter("fireSizeColumnName", "character", "SIZE_HA", NA, NA,
                    paste0("Name of the column that has fire size. ",
                           "In the latest dataset, its POLY_HA.")),
    defineParameter("fireYearColumnName", "character", "YEAR", NA, NA,
                    paste0("Name of the column that has fire size. ",
                           "In the latest dataset, its YEAR"))
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "firePoints", objectClass = "SpatialPointsDataFrame",
                 desc = paste0("Historical fire data in point form. Must contain fields 'CAUSE',
                               'YEAR', and 'SIZE_HA', or pass the parameters to identify those"),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "binary map of landscape flammbility"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list",
                 desc = "contains landscape attributes for each polygon"),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer",
                 desc = "template raster for raster GIS operations. Must be supplied by user with same CRS as studyArea")
  ),
  outputObjects = bind_rows(
   createsOutput(objectName = "scfmRegimePars", objectClass = "list", desc =  "Fire regime parameters for each polygon"),
   createsOutput(objectName = "firePoints", objectClass = "SpatialPointsDataFrame",
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
  if (is(sim$firePoints, "list")){
    tmp <- do.call(what = bind, args = list(sim$firePoints))
  } else {
    tmp <- sim$firePoints
  }
  if (length(sim$firePoints) == 0) {
    stop("there are no fires in your studyArea. Consider expanding the study Area")
  }
  #extract and validate fireCause spec

  fc <- P(sim)$fireCause

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

  tmp <- subset(tmp, get(P(sim)$fireYearColumnName) >= epoch[1] &
                  get(P(sim)$fireYearColumnName) <= epoch[2])

  epochLength <- as.numeric(epoch[2] - epoch[1] + 1)

  # Assign polygon label to SpatialPoints of fires object
  #should be specify the name of polygon layer? what if it PROVINCE or ECODISTRICT
  #tmp[["ECOREGION"]] <- sp::over(tmp, sim$studyArea[, "ECOREGION"])
  frpl <- sim$fireRegimePolys$PolyID
  tmp$PolyID <- sp::over(tmp, sim$fireRegimePolys) #gives studyArea row name to point
  tmp$PolyID <- tmp$PolyID$PolyID

  if (any(is.na(tmp$PolyID)))
    tmp <- tmp[!is.na(tmp$PolyID),] #have to remove NA points

  # sim$firePoints <- tmp # TM ~ Terrible terrible idea. Other modules might use this object and
  # espect it to be the list of fire points
  # You shouldn't modify it unless its a dataPrep module or smth and even so, its dangerous.
  # Just use your local variable, or save as a new sim.Object if this needs to go across modules.
  # I checked all SCFM modules, and only regime uses firePoints, so I am keeping the local
  # Object

  # firePolys <- unlist(sim$firePoints) # TM ~ Use tmp, as its already unlisted!
  firePolys <- tmp # TM ~ Use tmp, as its already unlisted!

  scfmRegimePars <- lapply(names(sim$landscapeAttr), FUN = calcZonalRegimePars,
                           firePolys = firePolys, landscapeAttr = sim$landscapeAttr,
                           firePoints = firePolys, epochLength = epochLength,
                           maxSizeFactor = P(sim)$empiricalMaxSizeFactor,
                           fireSizeColumnName = P(sim)$fireSizeColumnName)

  names(scfmRegimePars) <- names(sim$landscapeAttr)

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

  if (!suppliedElsewhere("studyArea", sim)) {

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
  ## this module has many dependencies that aren't sourced in .inputObjects
  ## this workaround prevents checksums updating due to daily name change of NFDB files
  if (!suppliedElsewhere("firePoints", sim)) {
    sim$firePoints <- getFirePoints_NFDB(url = extractURL("firePoints", sim),
                                        studyArea = sim$studyArea, rasterToMatch = sim$rasterToMatch,
                                        NFDB_pointPath = checkPath(file.path(dataPath(sim), "NFDB_point"),
                                                                   create = TRUE))
  }
  return(invisible(sim))
}
