calcZonalRegimePars <- function(polygonID, firePolys,
                                landscapeAttr,
                                firePoints,
                                epochLength,
                                maxSizeFactor,
                                fireSizeColumnName,
                                targetBurnRate = NULL,
                                targetMaxFireSize = NULL) {

  idx <- firePolys$PolyID == polygonID
  tmpA <- firePoints[firePoints$PolyID == as.numeric(polygonID),]
  landAttr <- landscapeAttr[[polygonID]]
  cellSize = landAttr$cellSize
  nFires <- dim(tmpA)[1]
  if (nFires == 0) {
    return(NULL)
  }
  rate <- nFires / (epochLength * landAttr$burnyArea)   # fires per ha per yr
  pEscape <- 0
  xBar <- 0 # mean fire size
  xMax <- 0
  lxBar <- NA
  maxFireSize <- cellSize   #note that maxFireSize has unit of ha NOT cells!!!
  xVec <- numeric(0)

  #check for user supplied defaults
  targetBurnRate <- targetBurnRate[polygonID]
  targetMaxFireSize <- targetMaxFireSize[polygonID]

  if (nFires > 0) {
    #calculate escaped fires
    #careful to subtract cellSize where appropriate

    xVec <- tmpA[[fireSizeColumnName]][tmpA[[fireSizeColumnName]] > cellSize]
    if (length(xVec) > 0) {
      pEscape <- length(xVec) / nFires
      xBar <- mean(xVec)
      lxBar <- mean(log(xVec))
      xMax <- max(xVec)
      xFireSize <- xBar
      zVec <- log(xVec / cellSize)
      if (length(zVec) < 25)
        warning(paste("Less than 25 \"large\" fires in zone", polygonID, ".",
                      "Estimates may be unstable.\n",
                      "\tConsider using a larger area and/or longer epoch.\n"))
      hdList <- HannonDayiha(zVec)  #defined in sourced TEutilsNew.R
      That <- hdList$That
      if (That == -1) {
        warning(
          sprintf(
            "Hannon-Dahiya convergence failure in zone %s.\n\tUsing sample maximum fire size",
            polygonID
          )
        )
        maxFireSize <- xMax * maxSizeFactor  #just to be safe, re-specify here
      } else {
        maxFireSize <- exp(That) * cellSize
        if (!(maxFireSize > xMax)) {
          warning(
            sprintf("Dodgy maxFireSize estimate in zone %s.\n\tUsing sample maximum fire size.",
                    polygonID)
          )
          maxFireSize <- xMax * maxSizeFactor
        }
        #missing BEACONS CBFA truncated at 2*xMax. Their reasons don't apply here.
      }
    } else {
      # there should be a way to pass non-zero defaults but I'm not sure whether we would specify by polygon
      # and if so, how, given the initial polygons may be modified during sliver removal
      message(paste("no fires larger than cellsize in ", polygonID, "."))
    }
  } else {
    message(paste("Insufficient data for polygon ", polygonID, ". Default values used."))
  }

  #verify estimation results are reasonable. That=-1 indicates convergence failure.
  #need to add a name or code for basic verification by Driver module, and time field
  #to allow for dynamic regeneration of disturbanceDriver pars.
  #browser()
  if (maxFireSize < 1) {
    warning("this can't happen")
    maxFireSize = cellSize
  }

  empiricalBurnRate <- sum(tmpA[[fireSizeColumnName]]) / (epochLength * landAttr$burnyArea)

  if (is.na(targetBurnRate) | is.null(targetBurnRate)) {
    ratio <- 1
  }

  if (!is.na(targetBurnRate) | is.null(targetBurnRate)) {
    ratio <-  targetBurnRate / empiricalBurnRate
    if (ratio >= 1) {
      newFireValues <- ratioPartition2(targetBurnRate = targetBurnRate,
                                       empiricalBurnRate = empiricalBurnRate,
                                       pEscape = pEscape,
                                       xBar = xBar,
                                       rate = rate)
      rate <- newFireValues$rate
      pEscape <- newFireValues$pEscape
      xBar <- newFireValues$xBar
    } else {
      warning("ratio cannot be < 1. Please make sure this does not happen")
    }
  }

  ## override maximum fire size if user supplied
  if (!is.na(targetMaxFireSize) | is.null(targetMaxFireSize)) {
    maxFireSize <- targetMaxFireSize
    xMax <- targetMaxFireSize
    ## TODO: add check that max is larger than mean, else stop
  }

  ## max fire size is returned twice - I think this is a backwards compatibility decision
  return(list(ignitionRate = rate,
              pEscape = pEscape,
              xBar = xBar,           ## mean fire size
              lxBar = lxBar,         ## mean log(fire size)
              xMax = xMax,           ## maximum observed size
              emfs_ha = maxFireSize, ## Estimated Maximum Fire Size in ha
              empiricalBurnRate = empiricalBurnRate
  ))
}
