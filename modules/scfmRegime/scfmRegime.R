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
  reqdPkgs = list("rgdal"),
  parameters = rbind(
    defineParameter("empiricalMaxSizeFactor", "numeric", 1.2, 1, 10, "scale xMax by this is HD estimator fails "),
    defineParameter("fireCause", "character", c("L"), NA_character_, NA_character_, "subset of c(H,H-PB,L,Re,U)"),
    defineParameter("fireEpoch", "numeric", c(1971,2000), NA, NA, "start of normal period")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "firePoints", objectClass = "SpatialPointsDataFrame", desc = "Historical fire data in point form. Must contain fields 'CAUSE', 'YEAR', and 'SIZE_HA'",
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "binary map of landscape flammbility"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = "contains landscape attributes for each polygon"),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer", desc = "template raster for raster GIS operations. Must be supplied by user with same CRS as studyArea")
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
    Init(sim)
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
    }
  return(invisible(sim))
}

Init <- function(sim) {

  tmp <- sim$firePoints
  if (length(sim$firePoints) == 0) {
    stop("there are no fires in your studyArea. Consider expanding the study Area")
  }
  #extract and validate fireCause spec

  fc <- P(sim)$fireCause
  #should verify CAUSE is a column in the table...
  causeSet <- if (is.factor(tmp$CAUSE)) levels(tmp$CAUSE) else unique(tmp$CAUSE)
  if (any(!(fc %in% causeSet)))
    stop("illegal fireCause: ", fc)
  tmp <- subset(tmp, CAUSE %in% fc)

  #extract and validate fireEpoch
  epoch <- P(sim)$fireEpoch
  if (length(epoch) != 2 ||
      !is.numeric(epoch) || any(!is.finite(epoch)) || epoch[1] > epoch[2])
    stop("illegal fireEpoch: ", epoch)
  tmp <- subset(tmp, YEAR >= epoch[1] & YEAR <= epoch[2])

  epochLength <- as.numeric(epoch[2] - epoch[1] + 1)

  # Assign polygon label to SpatialPoints of fires object
  #should be specify the name of polygon layer? what if it PROVINCE or ECODISTRICT
  #tmp[["ECOREGION"]] <- sp::over(tmp, sim$studyArea[, "ECOREGION"])
  tmp <- sim$firePoints

  #extract and validate fireCause spec
  fc <- P(sim)$fireCause
  #should verify CAUSE is a column in the table...
  causeSet <- if (is.factor(tmp$CAUSE)) levels(tmp$CAUSE) else unique(tmp$CAUSE)

  if (any(!(fc %in% causeSet)))
    stop("illegal fireCause: ", fc)
  tmp <- subset(tmp,CAUSE %in% fc)

  #extract and validate fireEpoch
  epoch <- P(sim)$fireEpoch
  if (length(epoch) != 2 || !is.numeric(epoch) || any(!is.finite(epoch)) || epoch[1] > epoch[2])
    stop("illegal fireEpoch: ", epoch)
  tmp <- subset(tmp, YEAR >= epoch[1] & YEAR <= epoch[2])

  epochLength <- as.numeric(epoch[2] - epoch[1] + 1)

  # Assign polygon label to SpatialPoints of fires object
  #should be specify the name of polygon layer? what if it PROVINCE or ECODISTRICT
  #tmp[["ECOREGION"]] <- sp::over(tmp, sim$studyArea[, "ECOREGION"])

  frpl <- sim$studyArea$PolyID
  tmp$PolyID <- sp::over(tmp, sim$studyArea[sim$studyArea$PolyID,]) #gives studyArea row name to point
  tmp$PolyID <- tmp$PolyID$PolyID

  tmp <- tmp[!is.na(tmp$PolyID),] #have to remove NA points
  sim$firePoints <- tmp

  firePolys <- unlist(sim$firePoints)

  scfmRegimePars <- lapply(names(sim$landscapeAttr), FUN = calcZonalRegimePars,
                               firePolys = firePolys, landscapeAttr = sim$landscapeAttr,
                               firePoints = sim$firePoints, epochLength = epochLength,
                               maxSizeFactor = P(sim)$empiricalMaxSizeFactor)

  names(scfmRegimePars) <- names(sim$landscapeAttr)

  nullIdx <- sapply(scfmRegimePars, is.null)
  if (any(nullIdx)){
    scfmRegimePars <- scfmRegimePars[-which(nullIdx)]
  }
  sim$scfmRegimePars <- scfmRegimePars

  return(invisible(sim))
}

calcZonalRegimePars <- function(polygonID, firePolys = firePolys, landscapeAttr = sim$landscapeAttr,
                                firePoints = sim$firePoints, epochLength = epochLength, maxSizeFactor) {

  idx <- firePolys$PolyID == polygonID
  tmpA <- firePoints[idx, ]
  landAttr <- landscapeAttr[[polygonID]]
  cellSize = landAttr$cellSize
  nFires <- dim(tmpA)[1]
  if (nFires == 0) {
    return(NULL)
  }
  rate <- nFires / (epochLength * landAttr$burnyArea)   # fires per ha per yr

  pEscape <- 0
  xBar <- 0
  xMax <- 0
   #NA might be better, but would take more downstream work SGC 15.10.2018
  maxFireSize <- lxBar <- NA
  maxFireSize <- cellSize
  xVec <- numeric(0)

  if (nFires > 0) {
    #calculate escaped fires
    #careful to subtract cellSize where appropriate
    xVec <- tmpA$SIZE_HA[tmpA$SIZE_HA > cellSize]

    if (length(xVec) > 0) {
      pEscape <- length(xVec) / nFires
      xBar <- mean(xVec)
      lxBar <- mean(log(xVec))
      xMax <- max(xVec)

      zVec <- log(xVec / cellSize)
      if (length(zVec) < 50)
        warning(
          sprintf(
            "Less than 50 \"large\" fires in zone %s. T estimates may be unstable.\n\tConsider using a larger area and/or longer epoch.",
            polygonID
          )
        )
      hdList <- HannonDayiha(zVec)  #defined in sourced TEutilsNew.R
      That <- hdList$That
      if (That == -1) {
        warning(
          sprintf(
            "Hannon-Dahiya convergence failure in zone %s.\n\tUsing sample maximum fire size",
            polygonID
          )
        )
        #browser()
        maxFireSize <- xMax * maxSizeFactor  #just to be safe, respecify here
      }
      else {
        maxFireSize <- exp(That) * cellSize
        if (!(maxFireSize > xMax)) {
          warning(
            sprintf("Dodgy maxFireSize estimate in zone %s.\n\tUsing sample maximum fire size.",polygonID)
          )
          maxFireSize <- xMax * maxSizeFactor
        }
        #missing BEACONS CBFA truncated at 2*xMax. Their reasons don't apply here.
      }
    }
  } else {
    message(paste("Insufficient data for polygon ", polygonID, ". Default values used."))
    #return(NULL)
  }

  #verify estimation results are reasonable. That=-1 indicates convergence failure.
  #need to addd a name or code for basic verification by Driver module, and time field
  #to allow for dynamic regeneration of disturbanceDriver pars.
  #browser()

  return(list(ignitionRate = rate,
              pEscape = pEscape,
              xBar = xBar,
              #mean fire size
              lxBar = lxBar,
              #mean log(fire size)
              xMax = xMax,
              #maximum observed size
              emfs = maxFireSize  #Estimated Maximum Fire Size in ha
              )
          )
}

.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using Ecodistrict 348")

    #source shapefile from ecodistict in input folder. Use ecodistrict 348
    studyAreaFilename <- file.path(dPath, "ecodistricts.shp")
    SA <- Cache(prepInputs,
                targetFile  = studyAreaFilename,
                fun = "raster::shapefile",
                url = extractURL(objectName = "studyArea"),
                archive = "ecodistrict_shp.zip",
                filename2 = TRUE,
                userTags = c(cacheTags, "studyArea"),
                destinationPath = file.path(dPath, "ecodistricts_shp", "Ecodistricts"))

    SA <- SA[SA$ECODISTRIC == 348, ]
    sim$studyArea <- SA
  }
  #this module has many dependencies that aren't sourced in .inputObjects
  if (!suppliedElsewhere("firePoints", sim)) {
    if (!dir.exists(file.path(dPath, "NFDB_point"))) {

      download.file(
        url = extractURL(objectName = "firePoints"),
        destfile = file.path(dPath, "NFDB_point.zip")
      )
      unzip(
        zipfile = file.path(dPath, "NFDB_point.zip"),
        exdir = file.path(dPath, "NFDB_point")
      )
    }

    #fire points file name changes daily so must be grepped
    zipContents <- list.files(file.path(dPath, "NFDB_point"),
                              all.files = TRUE,
                              full.names = TRUE)

    outFile <- grep(pattern = "*.shp$",
                    x = zipContents,
                    value = TRUE)

    #Reading this shapefile takes forever even when cached so combining these calls.
    #PrepInputs doesn't work here because we don't know the targetFile (name changes daily)
    #And we don't want the file written to the archive because that screws up the grep
    fireDownload <- function(SA, file = outFile) {
      firePoints <- raster::shapefile(file) %>%
      sp::spTransform(CRSobj = crs(SA))
      firePoints <- postProcess(firePoints, studyArea = SA, rasterToMatch = sim$rasterToMatch,
                                filename2 = file.path(dPath, "firePoints_SA.shp"),
                                overwrite = TRUE)
      return(firePoints)
    }

    sim$firePoints <- Cache(fireDownload, SA = sim$studyArea, file = outFile)
  }

  return(invisible(sim))
}
