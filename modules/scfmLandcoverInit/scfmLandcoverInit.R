stopifnot(packageVersion("SpaDES") >= "0.99.0")

defineModule(sim,list(
    name = "scfmLandcoverInit",
    description = "Takes the LCC05 classification of 39 land cover classes, and reclassifies it to flammable and inflammable [1,0]",
    keywords = c("fire", "LCC05", "land cover classification 2005", "BEACONs"),
    childModules = character(),
    authors = c(
      person(c("Eliot", "J", "B"),"McIntire", email = "Eliot.McIntire@NRCan.gc.ca", role = c("aut", "cre")),
      person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut"))),
    version = numeric_version("0.1.0"),
    spatialExtent = raster::extent(rep(NA_real_, 4)),
    timeframe = as.POSIXlt(c("2005-01-01", NA)),
    documentation = list("README.txt", "scfmLandCoverInit.Rmd"),
    timeunit = "year",
    citation = list(),
    reqdPkgs = list("raster", "reproducible"),
    parameters = rbind(
      defineParameter(".plotInitialTime", "numeric", 0, NA, NA, desc = "Initial time for plotting"),
      defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, desc = "Interval between plotting"),
      defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, desc = "Initial time for saving"),
      defineParameter(".saveIntervalXXX", "numeric", NA_real_, NA, NA, desc = "Interval between save events"),
      defineParameter("useCache", "logical", TRUE, NA, NA, desc = "Use cache"),
      defineParameter("neighbours", "numeric", 8, NA, NA, desc = "Number of immediate cell neighbours"),
      defineParameter(".crsUsed", "CRS", raster::crs(paste(
          "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")),
          NA, NA, desc = "CRS to be used. Defaults to the default vegMap projection")
    ),
    inputObjects = bind_rows(
      expectsInput(objectName = "studyArea0", objectClass = "SpatialPolygonsDataFrame", desc = "",
                   sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
      expectsInput(objectName = "vegMap", objectClass = "RasterLayer", desc = "",
        sourceURL = "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip")
    ),
    outputObjects = bind_rows(
      createsOutput(objectName = "cellsByZone", objectClass = "data.frame", desc = ""),
      createsOutput(objectName = "flammableMap", objectClass = "RasterLayer", desc = ""),
      createsOutput(objectName = "landscapeAttr", objectClass = "list", desc = ""),
      createsOutput(objectName = "studyArea", objectClass = "SpatialPolygonsDataLayer", desc = "")
    )
  )
)

doEvent.scfmLandcoverInit = function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
         init = {
           #sim <- scfmLandcoverInitCacheFunctions(sim)
           sim <- Init(sim)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmLandcoverInit", "plot")
           sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "scfmLandcoverInit", "save")
         },
         plot =  {
           Plot(sim$vegMap, new = TRUE)
           Plot(sim$flammableMap, legend = FALSE) # this is failing probably due to a bug in Plot
           # EJM is working on it 20160224
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
    neighMap <- Cache(focal, x = 1 - flammableMap, w = w, na.rm = TRUE) #default function is sum(...,na.rm)

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
      sum(1 - getValues(flammableMap)[x[, 1]], na.rm = TRUE) #sums flammable pixels in FRI polygons
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
Init <- function(sim) {

  if (is.null(sim$studyArea)) {
    sim$studyArea <- sim$studyArea0
  }

  if (!identical(crs(sim$studyArea), P(sim)$.crsUsed)) {

   sim$studyArea <- spTransform(sim$studyArea, CRSobj = P(sim)$.crsUsed)
  }

  if (is.null(sim$studyArea$PolyID)) {
    sim$studyArea$PolyID <- row.names(sim$studyArea)
  }

  nonFlammClasses <- c(36, 37, 38, 39)
  oldClass <- 0:39
  newClass <- ifelse(oldClass %in% nonFlammClasses, 1, 0)   #1 codes for non flammable
  #see mask argument for SpaDES::spread()
  flammableTable <- cbind(oldClass, newClass)

  sim$flammableMap <- makeFlammableMap(sim$vegMap, flammableTable, ls(sim))

  # This makes sim$landscapeAttr & sim$cellsByZone
  outs <- Cache(genFireMapAttr,
                sim$flammableMap,
                sim$studyArea,
                P(sim)$neighbours)
  list2env(outs, envir = envir(sim)) # move 2 objects to sim environment without copy
  return(invisible(sim))
}

testFun <- function(x) {
  sum(na.omit(x) == 1)
}

makeFlammableMap <- function(vegMap, flammableTable, lsSimObjs) {

  flammableMap <- ratify(reclassify(vegMap, flammableTable, count = TRUE))

  if ("Mask" %in% lsSimObjs) {
    flammableMap <- flammableMap * sim$Mask # don't pass in sim$Mask explicitly to fn so not assessed in Cache
  }
  #the count options should cause that "a column with frequencies is added.

  #setColors(sim$flammableMap, n=2) <- c("blue","red")
  setColors(flammableMap, 2) <- colorRampPalette(c("blue", "red"))(2)
  #flammableMap <- writeRaster(flammableMap, filename = "flammableMap.tif",
  #                            datatype = "INT2U", overwrite=TRUE)
  return(flammableMap)
}

.inputObjects <- function(sim) {

  dPath <- dataPath(sim) #where files will be downloaded
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyArea0", sim)) {
    message("study area not supplied. Using Ecodistrict 348")

    #source shapefile from ecodistict in input folder. Use ecodistrict 348
    studyAreaFilename <- file.path(dPath, "ecodistricts.shp")

    SA <- Cache(prepInputs,
                targetFile  = studyAreaFilename,
                fun = "raster::shapefile",
                url = extractURL(objectName = "studyArea0"),
                archive = "ecodistrict_shp.zip",
                filename2 = TRUE,
                userTags = c(cacheTags, "studyArea0"),
                destinationPath = file.path(dPath, "ecodistricts_shp", "Ecodistricts"))

    SA <- SA[SA$ECODISTRIC == 348, ]
    sim$studyArea0 <- SA
    sim$studyArea <- SA

  }

  sim$studyArea0 <- spTransform(sim$studyArea0, CRSobj = P(sim)$.crsUsed) #must happen regardless if supplied
  sim$studyArea <- sim$studyArea0

  if (!suppliedElsewhere("vegMap", sim)) {
    message("vegMap not supplied. Using default LandCover of Canada 2005 V1_4a")

    vegMapFilename <- file.path(dPath, "LCC2005_V1_4a.tif")

    vegMap <- Cache(prepInputs, targetFile = vegMapFilename,
                    url = extractURL(objectName = "vegMap"),
                    archive = "LandCoverOfCanada2005_V1_4.zip",
                    destinationPath = dPath,
                    studyArea = sim$studyArea0,
                    overwrite = TRUE,
                    filename2 = TRUE,
                    userTags = c(cacheTags, "vegMap"),
                    targetCRS = P(sim)$.crsUsed)

    if (!identicalCRS(crs(vegMap), P(sim)$.crsUsed)) {
      vegMap <- raster::projectRaster(vegMap, crs = P(sim)$.crsUsed) #This can screw up resolution. Need to revisit when rtm is available
    }
    sim$vegMap <- vegMap

  }
  return(invisible(sim))
}
