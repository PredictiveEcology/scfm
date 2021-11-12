defineModule(sim,list(
    name = "scfmLandcoverInit",
    description = "Takes the LCC2010 classification of land cover classes, and reclassifies it to flammable and inflammable [1,0]",
    keywords = c("fire", "LCC2010", "land cover classification 2010", "BEACONs"),
    childModules = character(),
    authors = c(
      person(c("Eliot", "J", "B"),"McIntire", email = "Eliot.McIntire@canada.ca", role = c("aut", "cre")),
      person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
      person("Ian", "Eddy", email = "ian.eddy@canada.ca", role = c("aut"))),
    version = numeric_version("0.1.0"),
    spatialExtent = raster::extent(rep(NA_real_, 4)),
    timeframe = as.POSIXlt(c("2005-01-01", NA)),
    documentation = list("README.txt", "scfmLandcoverInit.Rmd"),
    timeunit = "year",
    citation = list(),
    reqdPkgs = list("fasterize", "purrr", "raster", "rgeos", "sf",
                    "PredictiveEcology/LandR",
                    "PredictiveEcology/reproducible"),
    parameters = rbind(
      defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, desc = "Initial time for plotting"),
      defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, desc = "Interval between plotting"),
      defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, desc = "Initial time for saving"),
      defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, desc = "Interval between save events"),
      defineParameter("removeSlivers", "logical", TRUE, NA, NA,
                      desc = "automatically merge polygons smaller than sliverThreshold param"),
      defineParameter("useCache", "logical", TRUE, NA, NA, desc = "Use cache"),
      defineParameter("neighbours", "numeric", 8, NA, NA, desc = "Number of immediate cell neighbours"),
      defineParameter("sliverThreshold", "numeric", NA, NA, NA,
                      desc = "fire regime polygons with area less than this number will be merged")
    ),
    inputObjects = bindrows(
      expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "",
                   sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
      expectsInput(objectName = "vegMap", objectClass = "RasterLayer", desc = "Landcover to build flammability map",
                   sourceURL = NA), ## uses default url specified in LandR::prepInputsLCC()
      expectsInput(objectName = "flammableMap", objectClass = "RasterLayer",
                   desc = "binary flammability map - defaults to using LandR::prepInputsLCC"),
      expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer",
                   desc = "template raster for raster GIS operations. Must be supplied by user"),
      expectsInput(objectName = "fireRegimePolys", objectClass = "SpatialPolygonsDataFrame",
                   desc = paste("Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
                                "Must have numeric field 'PolyID' or it will be created for individual polygons"),
                   sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip")
    ),
    outputObjects = bindrows(
      createsOutput(objectName = "cellsByZone", objectClass = "data.frame",
                    desc = "explains which raster cells are in which polygon"),
      createsOutput(objectName = "landscapeAttr", objectClass = "list", desc = "list of polygon attributes inc. area"),
      createsOutput(objectName = 'fireRegimePolys', objectClass = "SpatialPolygonsDataFrame",
                    desc = paste("areas to calibrate individual fire regime parameters",
                                 "modified by removing slivers and adding PolyID field where applicable")),
      createsOutput(objectName = "fireRegimeRas", objectClass = "RasterLayer",
                    desc = "Rasterized version of fireRegimePolys with values representing polygon ID")
    )
))

doEvent.scfmLandcoverInit = function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
         init = {
           sim <- Init(sim)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmLandcoverInit", "plot")
           sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "scfmLandcoverInit", "save")
         },
         plot =  {
           Plot(sim$fireRegimeRas, title = c("fire regimes"), new = TRUE)
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
  message("checking sim$fireRegimePolys for sliver polygons...")

  if (sf::st_is_longlat(sim$fireRegimePolys)) {
    warning("fireRegimePolys has a long/lat projection. This is untested")
    if (is.na(P(sim)$sliverThreshold)) {
      stop("You must supply P(sim)$sliverThreshold for fireRegimePolys with ")
    }
  }
  sim$fireRegimePolys$trueArea <- round(rgeos::gArea(sim$fireRegimePolys, byid = TRUE), digits = 0)
  if (is.na(P(sim)$sliverThreshold)) {
    sim@params[[currentModule(sim)]]$sliverThreshold <- 1e4 * 1e4 #100km2
  }
  if (P(sim)$removeSlivers) {
    if (any(sim$fireRegimePolys$trueArea < P(sim)$sliverThreshold)) {
      message("sliver polygon(s) detected. Merging to their nearest valid neighbour")
      sim$fireRegimePolys <- Cache(deSliver, sim$fireRegimePolys, threshold = P(sim)$sliverThreshold,
                                   userTags = c("deSliver", currentModule(sim)))
    }
  }
  if (is.null(sim$fireRegimePolys$PolyID)) {
    if (is.null(sim$fireRegimePolys$ECOREGION)) {
      warning(paste("no PolyID variable in fireRegimePolys - this is used to determine unique fire regimes.",
                    "Using row.names as default"))
      sim$fireRegimePolys$PolyID <- row.names(sim$fireRegimePolys)
    } else {
      #default fireRegimePolys were likely used
      sim$fireRegimePolys$PolyID <- as.numeric(sim$fireRegimePolys$ECOREGION)
    }
  }

  if (length(unique(sim$fireRegimePolys$PolyID)) != length(sim$fireRegimePolys)) {
    stop("mismatch between PolyID and fireRegimePolys. Must be 1 PolyID value per multipolygon object")
  }

  temp <- sf::st_as_sf(sim$fireRegimePolys)
  temp$PolyID <- as.numeric(temp$PolyID) #fasterize needs numeric; row names must stay char

  fireRegimeRas <- fasterize::fasterize(sf = temp, raster = sim$vegMap, field = "PolyID")
 #fireRegimeRas is inheriting the color table of vegMap. workaround below
  sim$fireRegimeRas <- raster(fireRegimeRas)
  sim$fireRegimeRas <- setValues(sim$fireRegimeRas, getValues(fireRegimeRas))

  # This makes sim$landscapeAttr & sim$cellsByZone
  outs <- Cache(genFireMapAttr,
                flammableMap = sim$flammableMap,
                fireRegimePolys = sim$fireRegimePolys,
                neighbours = P(sim)$neighbours,
                userTags = c(currentModule(sim), "genFireMapAttr"))

  sim$landscapeAttr <- outs$landscapeAttr
  sim$cellsByZone <- outs$cellsByZone
  return(invisible(sim))
}

genFireMapAttr <- function(flammableMap, fireRegimePolys, neighbours) {
  #calculate the cell size, total area, and number of flammable cells, etc.
  #All areas in ha
  cellSize <- prod(res(flammableMap)) / 1e4 # in ha

  if (neighbours == 8)
    w <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3, ncol = 3)
  else if (neighbours == 4)
    w <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
  else
    stop("illegal neighbours specification")

  makeLandscapeAttr <- function(flammableMap, weight, fireRegimePolys) {
    neighMap <- focal(x = flammableMap, w = w, na.rm = TRUE) #default function is sum(...,na.rm)

    # extract table for each polygon
    valsByPoly <- raster::extract(neighMap, fireRegimePolys, cellnumbers = TRUE)
    valsByPoly <- lapply(valsByPoly, na.omit)
    names(valsByPoly) <- fireRegimePolys$PolyID
    uniqueZoneNames <- fireRegimePolys$PolyID #get unique zones.
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

  landscapeAttr <- makeLandscapeAttr(flammableMap, w, fireRegimePolys)

  cellsByZoneFn <- function(flammableMap, landscapeAttr) {

    cellsByZone <- data.frame(cell = 1:ncell(flammableMap), zone = NA_character_, stringsAsFactors = FALSE)

    for (x in names(landscapeAttr)) {
      cellsByZone[landscapeAttr[[x]]$cellsByZone, "zone"] <- x
      }
    return(cellsByZone)
  }

  cellsByZone <- cellsByZoneFn(flammableMap, landscapeAttr)

  return(invisible(list(landscapeAttr = landscapeAttr, cellsByZone = cellsByZone)))
}

### template initilization

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using random polygon in Alberta")
    #TODO: remove LandR once this is confirmed working
    studyArea <- LandR::randomStudyArea(size = 15000000000, seed = 23654)
    sim$studyArea <- studyArea
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message(paste("rasterToMatch not supplied. generating from LCC2010 using studyArea CRS",
                  " - It is strongly recommended to supply a rasterToMatch"))

    sim$rasterToMatch <- LandR::prepInputsLCC(year = 2010,
                                   destinationPath = dPath,
                                   studyArea = sim$studyArea,
                                   useSAcrs = TRUE,
                                   filename2 = NULL,
                                   overwrite = TRUE,
                                   userTags = c(cacheTags, "rasterToMatch"))
  }

  vegMapSupplied <- TRUE
  if (!suppliedElsewhere("vegMap", sim)) {
    vegMapSupplied <- FALSE
    message("vegMap not supplied. Using default 2010 LandCover of Canada")
    sim$vegMap <- LandR::prepInputsLCC(year = 2010,
                                       destinationPath = dPath,
                                       studyArea = sim$studyArea,
                                       rasterToMatch = sim$rasterToMatch,
                                       filename2 = NULL,
                                       useCache = TRUE,
                                       userTags = c(cacheTags, "vegMap"))
  }

  if (!suppliedElsewhere("flammableMap", sim)) {
    if (vegMapSupplied) {
      stop("vegMap supplied but flammableMap is not. Please provide both or neither")
    }
    message("flammableMap not supplied. veg map to create flammableMap")
    flammableMap <- defineFlammable(sim$vegMap,
                                    mask = sim$rasterToMatch,
                                    nonFlammClasses = c(13, 16, 17, 18, 19),
                                    filename2 = NULL)
    sim$flammableMap <- setValues(raster(flammableMap), flammableMap[]) #this removes colour assignment
  }

  if (!suppliedElsewhere("fireRegimePolys", sim)) {
    message("fireRegimePolys not supplied. Using default ecoregions of Canada")
    #cannot use prepInputs with a vector for prepInputs - unreliable w/ GDAL
    sim$fireRegimePolys <- prepInputs(url = extractURL("fireRegimePolys", sim),
                                      destinationPath = dPath,
                                      studyArea = sim$studyArea,
                                      # rasterToMatch = sim$rasterToMatch,
                                      filename2 = NULL,
                                      userTags = c(cacheTags, "fireRegimePolys"))
    sim$fireRegimePolys <- spTransform(sim$fireRegimePolys, CRSobj = crs(sim$rasterToMatch))

    #this is now necessary due to the gridded nature of ecoregion file
    sim$fireRegimePolys <- rgeos::gUnaryUnion(spgeom = sim$fireRegimePolys, id = sim$fireRegimePolys$ECOREGION)
    sim$fireRegimePolys$ECOREGION <- row.names(sim$fireRegimePolys)
  }

  return(invisible(sim))
}
