defineModule(sim, list(
  name = "scfmRegime",
  description = "estimates fire regime parameters for BEACONs a la Steve's method",
  keywords = c("fire regime", "BEACONs"),
  authors = c(person(c("Steven", "G."), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmRegime.Rmd"),
  reqdPkgs = list("rgdal"),
  parameters = rbind(
    defineParameter("fireCause", "character", c("L"), NA_character_, NA_character_, "subset of c(H,H-PB,L,Re,U)"),
    defineParameter("fireEpoch", "numeric", c(1971,2000), NA, NA, "start of normal period"),
    defineParameter("empiricalMaxSizeFactor", "numeric", 1.2,1, 10, "scale xMax by this is HD estimator fails ")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "firePoints", objectClass = "SpatialPointsDataFrame", desc = "",
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = ""),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = "")
  ),
  outputObjects = bind_rows(
   createsOutput(objectName = "scfmRegimePars", objectClass = "list", desc =  "")
  )
))


## event types
#   - type `init` is required for initiliazation

doEvent.scfmRegime = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType == "init") {
    Init(sim)
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
    }
  return(invisible(sim))
}

Init <- function(sim) {
  # verify studyArea is in same projection
  if (!identicalCRS(sim$studyArea, sim$vegMap)) {
    sim$studyArea <- Cache(spTransform, sim$studyArea, CRSobj = crs(sim$vegMap))
  }

  return(invisible(sim))
}

calcZonalRegimePars <- function(polygonID, firePolys = firePolys, landscapeAttr = sim$landscapeAttr,
                                firePoints = sim$firePoints, epochLength = epochLength) {

  idx <- firePolys$PolyID == polygonID
  tmpA <- firePoints[idx, ] #STOP HERE
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
            "Less than 50 \"large\" fires in zone %s. T estimates may be unstable.\
            \n\tConsider using a larger area and/or longer epoch.",
            polygonID
          )
        )
      hdList <- HannonDayiha(zVec)  #defined in sourced TEutilsNew.R
      That <- hdList$That
      if (That == -1) {
        warning(
          sprintf(
            "Hannon-Dahiya convergence failure in zone %s.\n
            \tUsing sample maximum fire size",
            polygonID
          )
        )
        maxFireSize <- xMax * P(sim)$empiricalMaxSizeFactor  #just to be safe, respecify here
      }
      else {
        maxFireSize <- exp(That) * cellSize
        if (!(maxFireSize > xMax)) {
          warning(
            sprintf(
              "Dodgy maxSize estimate in zone %s.\n\tUsing sample maximum fire size.",
              polygonID
            )
          )
          maxFireSize <- xMax * P(sim)$empiricalMaxSizeFactor
        }
        #missing BEACONS CBFA truncated at 2*xMax. Their reasons don't apply here.
      }
    }
  } else {
    return(NULL)
  }

  #verify estimation results are reasonable. That=-1 indicates convergence failure.
  #
  #need to addd a name or code for basic verification by Driver module, and time field
  #to allow for dynamic regeneration of disturbanceDriver pars.
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

    zipContents <- list.files(file.path(dPath, "NFDB_point"),
                              all.files = TRUE,
                              full.names = TRUE)

    outFile <- grep(pattern = "*.shp$",
                    x = zipContents,
                    value = TRUE)

    #browser()
    firePoints <- Cache(rgdal::readOGR, outFile)
    firePoints <- spTransform(firePoints, CRSobj = crs(sim$studyArea))
    firePoints <- firePoints[sim$studyArea, ]

    sim$firePoints <- firePoints
  }

  if (!suppliedElsewhere("scfmRegimePars", sim)) {

    tmp <- sim$firePoints

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
    # browser() #this was line 90 in the master branch
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
    sim$firePoints <- tmp

    # Hack to make a study area level cellSize ... TODO -- this should be removed from landscapeAttr
    cellSize <- sim$landscapeAttr[[1]]$cellSize

    firePolys <- unlist(sim$firePoints)


    sim$scfmRegimePars <- lapply(names(sim$landscapeAttr), FUN = calcZonalRegimePars,
                                 firePolys = firePolys, landscapeAttr = sim$landscapeAttr,
                                 firePoints = sim$firePoints, epochLength = epochLength)

    names(sim$scfmRegimePars) <- names(sim$landscapeAttr)

    sim$scfmRegimePars <- scfmRegimePars[-which(sapply(scfmRegimePars, is.null))]
  }

  return(invisible(sim))
}
