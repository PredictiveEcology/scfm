calcZonalRegimePars <- function(polygonID, firePolys,
                                landscapeAttr,
                                firePoints,
                                epochLength, 
                                maxSizeFactor,
                                fireSizeColumnName,
                                targetBurnRate) {
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
  xBar <- 0 # mean fire size
  xMax <- 0
  lxBar <- NA
  maxFireSize <- cellSize   #note that maxFireSize has unit of ha NOT cells!!!
  xVec <- numeric(0)
  xFireSize <- 0 
  
  if (nFires > 0) {
    #calculate escaped fires
    #careful to subtract cellSize where appropriate
    # xVec <- tmpA$SIZE_HA[tmpA$SIZE_HA > cellSize]fireSizeColumnName # Hardcoded!! Breaks
    # as soon as you use another fire points database
    xVec <- tmpA[[fireSizeColumnName]][tmpA[[fireSizeColumnName]] > cellSize]
    
    if (length(xVec) > 0) {
      pEscape <- length(xVec) / nFires
      xBar <- mean(xVec)
      lxBar <- mean(log(xVec))
      xMax <- max(xVec)
      xFireSize <- mean(xFireSize)
      zVec <- log(xVec / cellSize)
      if (length(zVec) < 30)
        warning(paste("Less than 30 \"large\" fires in zone", polygonID, ".",
                      "T estimates may be unstable.\n",
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
        maxFireSize <- xMax * maxSizeFactor  #just to be safe, respecify here
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
      #TODO Default values need to be used, except they are not being used here! Just a message saying
      # they are. But they are all zeroed, NOT default! This is NOT producing fires!
      message(paste("no fires larger than cellsize in ", polygonID, "."))
    }
  } else {
    message(paste("Insufficient data for polygon ", polygonID, ". Default values used."))
  }
  
  #verify estimation results are reasonable. That=-1 indicates convergence failure.
  #need to addd a name or code for basic verification by Driver module, and time field
  #to allow for dynamic regeneration of disturbanceDriver pars.
  #browser()
  if (maxFireSize < 1){
    warning("this can't happen")
    maxFireSize = cellSize
  }
  
  burnRate <- (nFires * xFireSize) / (epochLength * landAttr$burnyArea) 
  
  if  (!is.na(targetBurnRate)){
    ratio <- targetBurnRate / landAttr$burnyArea
    xFireSize <- xFireSize * ratio
  }
  
  return(list(ignitionRate = rate,
              pEscape = pEscape,
              xBar = xBar,
              #mean fire size
              lxBar = lxBar,
              #mean log(fire size)
              xMax = xMax,
              xFireSize = xFireSize,
              burnRate = burnRate,
              #maximum observed size
              emfs_ha = maxFireSize  #Estimated Maximum Fire Size in ha
  )
  )
}
