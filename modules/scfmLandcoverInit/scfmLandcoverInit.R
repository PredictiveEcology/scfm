stopifnot(packageVersion("SpaDES") >= "0.99.0")

defineModule(sim,list(
    name = "scfmLandcoverInit",
    description = "Takes the LCC05 classification of 39 land cover classes, and reclassifies it to flammable and inflammable [1,0]",
    keywords = c("fire", "LCC05", "land cover classification 2005", "BEACONs"),
    childModules = character(),
    authors = c(
      person(c("Eliot", "J", "B"),"McIntire", email = "Eliot.McIntire@canada.ca", role = c("aut", "cre")),
      person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut"))),
    version = numeric_version("0.1.0"),
    spatialExtent = raster::extent(rep(NA_real_, 4)),
    timeframe = as.POSIXlt(c("2005-01-01", NA)),
    documentation = list("README.txt", "scfmLandcoverInit.Rmd"),
    timeunit = "year",
    citation = list(),
    reqdPkgs = list("raster", "reproducible", "PredictiveEcology/LandR@development", "fasterize", "sf"),
    parameters = rbind(
      defineParameter(".plotInitialTime", "numeric", 0, NA, NA, desc = "Initial time for plotting"),
      defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, desc = "Interval between plotting"),
      defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, desc = "Initial time for saving"),
      defineParameter(".saveIntervalXXX", "numeric", NA_real_, NA, NA, desc = "Interval between save events"),
      defineParameter("useCache", "logical", TRUE, NA, NA, desc = "Use cache"),
      defineParameter("neighbours", "numeric", 8, NA, NA, desc = "Number of immediate cell neighbours")
    ),
    inputObjects = bind_rows(
      expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "",
                   sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
      expectsInput(objectName = "vegMap", objectClass = "RasterLayer", desc = "",
        sourceURL = "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip"),
      expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer", desc = "template raster for raster GIS operations. Must be supplied by user")
    ),
    outputObjects = bind_rows(
      createsOutput(objectName = "cellsByZone", objectClass = "data.frame", desc = ""),
      createsOutput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "map of landscape flammability"),
      createsOutput(objectName = "landscapeAttr", objectClass = "list", desc = "list of polygon attributes inc. area"),
      createsOutput(objectName = "studyAreaRas", objectClass = "RasterLayer", desc = "Rasterized version of study Area with values representing polygon ID")
    )
  )
)

doEvent.scfmLandcoverInit = function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
         init = {

           sim <- Init(sim)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmLandcoverInit", "plot")
           sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "scfmLandcoverInit", "save")
         },
         plot =  {
           Plot(sim$vegMap, new = TRUE)
           Plot(sim$flammableMap, legend = FALSE)
           # schedule future event(s)
           sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmLandcoverInit", "plot")

         },
         save = {

           sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "scfmLandcoverInit", "save")
         },
         warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                       "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = "")))

  return(invisible(sim))
}
Init <- function(sim) {
  if (class(sim$studyArea) == "SpatialPolygons") {
    stop("studyArea must be a SpatialPolygonsDataFrame")
  }
  if (is.null(sim$studyArea$PolyID)) {
    sim$studyArea$PolyID <- row.names(sim$studyArea)
  }

  temp <- sf::st_as_sf(sim$studyArea)
  sim$studyAreaRas <- fasterize(sf = temp, raster = sim$vegMap, field = "PolyID")

  sim$flammableMap <- LandR::defineFlammable(sim$vegMap, filename2 = NULL)
  # setColors(sim$flammableMap, 2) <- colorRampPalette(c("red", "blue"))(2)
  # This makes sim$landscapeAttr & sim$cellsByZone
  outs <- Cache(genFireMapAttr,
                sim$flammableMap,
                sim$studyArea,
                P(sim)$neighbours)
  list2env(outs, envir = envir(sim)) # move 2 objects to sim environment without copy
  return(invisible(sim))
}

genFireMapAttr <- function(flammableMap, studyArea, neighbours) {
  #calculate the cell size, total area, and number of flammable cells, etc.
  #
  #All areas in ha
  #

  cellSize <- prod(res(flammableMap)) / 1e4 # in ha

  if (neighbours == 8)
    w <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3, ncol = 3)
  else if (neighbours == 4)
    w <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
  else
    stop("illegal neighbours specification")

  makeLandscapeAttr <- function(flammableMap, weight, studyArea) {
    neighMap <- Cache(focal, x = flammableMap, w = w, na.rm = TRUE) #default function is sum(...,na.rm)

    # extract table for each polygon
    valsByPoly <- Cache(raster::extract, neighMap, studyArea, cellnumbers = TRUE)
    valsByPoly <- lapply(valsByPoly, na.omit)


    names(valsByPoly) <- studyArea$PolyID
    uniqueZoneNames <- studyArea$PolyID #get unique zones.
    valsByZone <- lapply(uniqueZoneNames, function(ecoName) {
      aa <- valsByPoly[names(valsByPoly) == ecoName]
      if (is.list(aa))
        aa <- do.call(rbind, aa)
      return(aa)
    })
    names(valsByZone) <- uniqueZoneNames

    # Derive frequency tables of number of flammable cells, per polygon type, currently ECOREGION
    nNbrs <- lapply(valsByZone, function(x) {
      nNbrs <- tabulate(x[, 2] + 1, 9)#depends on sfcmLandCoverInit
      names(nNbrs) <- 0:8
      return(nNbrs)
    })


    nFlammable <- lapply(valsByZone, function(x) {

      sum(getValues(flammableMap)[x[, 1]], na.rm = TRUE) #sums flammable pixels in FRI polygons
    })

    landscapeAttr <- purrr::transpose(list(cellSize = rep(list(cellSize), length(nFlammable)),
                                           nFlammable = nFlammable,
                                           nNbrs = nNbrs,
                                           cellsByZone = lapply(valsByZone, function(x) x[, 1])))

      landscapeAttr <- lapply(landscapeAttr, function(x) {
        append(x, list(burnyArea = x$cellSize * x$nFlammable))
      })
      names(landscapeAttr) <- names(valsByZone)

    return(landscapeAttr)
  }

  landscapeAttr <- Cache(makeLandscapeAttr, flammableMap, w, studyArea)

  cellsByZoneFn <- function(flammableMap, landscapeAttr) {

    cellsByZone <- data.frame(cell = 1:ncell(flammableMap), zone = NA_character_, stringsAsFactors = FALSE)

    for (x in names(landscapeAttr)) {
      cellsByZone[landscapeAttr[[x]]$cellsByZone, "zone"] <- x
      }
    return(cellsByZone)
  }

  cellsByZone <- Cache(cellsByZoneFn, flammableMap, landscapeAttr)

  return(invisible(list(landscapeAttr = landscapeAttr, cellsByZone = cellsByZone)))
}

### template initilization

.inputObjects <- function(sim) {

  dPath <- dataPath(sim) #where files will be downloaded
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using random polygon in Alberta")
    #TODO: remove LandR once this is confirmed working
    studyArea <- LandR::randomStudyArea(size = 15000000000, seed = 23654)
    sim$studyArea <- studyArea
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message("rasterToMatch not supplied. generating from LCC2005 using studyArea CRS")

    rasterToMatch <- LandR::prepInputsLCC(year = 2005,
                                               destinationPath = dPath,
                                               studyArea = sim$studyArea,
                                               useSAcrs = TRUE,
                                               filename2 = TRUE,
                                               overwrite = TRUE,
                                               userTags = c("cacheTags", "rasterToMatch"))
    sim$rasterToMatch <- rasterToMatch

  }


  if (!suppliedElsewhere("vegMap", sim)) {
    message("vegMap not supplied. Using default LandCover of Canada 2005 V1_4a")

    sim$vegMap <- prepInputsLCC(year = 2005,
                                destinationPath = dPath,
                                studyArea = sim$studyArea,
                                rasterToMatch = sim$rasterToMatch,
                                filename2 = TRUE,
                                overwrite = TRUE,
                                userTags = c("cacheTags", "vegMap"))
    }

  return(invisible(sim))
}
