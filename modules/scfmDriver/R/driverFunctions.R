#Buffers polygon, generates index raster
genSimLand <- function(coreLand, buffDist) {
  tempDir <- tempdir()
  #Buffer study Area. #rbind had occasional errors before makeUniqueIDs = TRUE
  #TODO: Investigate why some polygons fail
  bfireRegimePoly <- buffer(coreLand, buffDist)
  if (!gIsValid(bfireRegimePoly, byid = FALSE)) {
    bfireRegimePoly <- gBuffer(bfireRegimePoly, width = 0)
  }
  bfireRegimePoly <-  gDifference(bfireRegimePoly, spgeom2 = coreLand, byid = FALSE)
  polyLandscape <- rbind.SpatialPolygons(coreLand, bfireRegimePoly, makeUniqueIDs = TRUE) #
  polyLandscape$zone <- c("core", "buffer")
  polyLandscape$Value <- c(1, 0)

  #Generate flammability raster
  landscapeLCC <- prepInputsLCC(destinationPath = tempDir, studyArea = polyLandscape, useSAcrs = TRUE)
  landscapeFlam <- defineFlammable(landscapeLCC)
  #Generate landscape Index raster
  polySF <- st_as_sf(polyLandscape)
  landscapeIndex <- fasterize(polySF, landscapeLCC, "Value")

  calibrationLandscape <- list(polyLandscape, landscapeIndex, landscapeLCC, landscapeFlam)
  names(calibrationLandscape) <- c("fireRegimePoly", "landscapeIndex", "lcc", "flammableMap")
  return(calibrationLandscape)
}



#make design
#dT <- data.frame("igLoc" = index, p0 = 0.1, p = 0.23)

#this version of makeDesign is the simplest possible...

makeDesign <- function(indices, targetN, pEscape = 0.1, pmin, pmax, q = 1) {
  #TODO: Fix makeDesign to work if polygons have no fires
  sampleSize <- round(targetN / pEscape)
  cellSample <- sample(indices, sampleSize, replace = TRUE)
  pVec <- runif(sampleSize)^q
  pVec <- pVec * (pmax - pmin) + pmin

  #derive p0 from escapeProb
  #steal code from scfmRegime and friends.

  p0 <- 1 - (1 - pEscape)^0.125  #assume 8 neighbours
  #the preceding approximation seems inadequate in practice.
  #when implemented in scfmDriver, make use of correct derivation of p0 from pEscape based on L
  Temp <- data.frame("igLoc" = cellSample, "p0" = p0, "p" = pVec)
  return(Temp)
}

executeDesign <- function(L, dT, maxCells) {
  # extract elements of dT into a three column matrix where column 1,2,3 = igLoc, p0, p

  iter <- 0
  f <- function(x, L, ProbRas) { ## L, P are rasters, passed by reference
    iter <<- iter + 1
    currentTime <- Sys.time()
    diffTime <- currentTime - startTime
    units(diffTime) <- "secs"
    timePer <- as.numeric(diffTime) / iter
    timeLeft <- (NROW(dT) - iter) * timePer
    timeLeft <- round(as.difftime(timeLeft, units = "mins")/60, 1)
    nrowDT <- NROW(dT)
    if (iter %% 200 == 0)
      message("  ", iter, " of ", nrowDT, " total; estimated time remaining: ",
              format(timeLeft, units = "mins"))


    threadsDT <- getDTthreads()
    setDTthreads(1)
    on.exit({setDTthreads(threadsDT)}, add = TRUE)

    i <- x[1]
    p0 <- x[2]
    p <- x[3]

    nbrs <- as.vector(SpaDES.tools::adj(x = L, i, pairs = FALSE, directions = 8))
    #nbrs < nbrs[which(L[nbrs] == 1)] #or this?
    nbrs <- nbrs[L[nbrs] == 1] #only flammable neighbours please. also, verify NAs excluded.
    #nbrs is a vector of flammable neighbours.
    nn <- length(nbrs)
    res <- c(nn, 0, 1)
    if (nn == 0)
      return(res) #really defaults
    #P is still flammableMap.

    ProbRas[nbrs] <- p0
    #Now it is 1, 0, p0, and NA
    spreadState0 <- SpaDES.tools::spread2(landscape = L,
                                          start = i,
                                          iterations = 1,
                                          spreadProb = ProbRas,
                                          asRaster = FALSE)

    tmp <- nrow(spreadState0)
    res[2:3] <- c(tmp - 1,tmp)
    if (tmp == 1) # the fire did not spread
      return(res)
    ProbRas[] <- L[]*p
    spreadState1 <- SpaDES.tools::spread2(landscape = L,
                                          start = spreadState0,
                                          spreadProb = ProbRas,
                                          asRaster = FALSE,
                                          maxSize = maxCells)
    #calculate return data
    res[3] <- nrow(spreadState1)
    return(res)
  }

  probRas <- raster(L)
  probRas[] <- L[]

  startTime <- Sys.time()
  res <- Cache(apply, dT, 1, f, L, ProbRas = probRas) # Parallelizing isn't efficient here. ~TM 15Feb19
  res <- data.frame("nNeighbours" = res[1,], "initSpreadEvents" = res[2,], "finalSize" = res[3,])

  #cbind dT and res, then select the columns we need
  x <- cbind(dT,res)
  x <- x[x$finalSize > 1,]
  return(x)
}

#This is a wrapper on makeDesign and executeDesign
makeAndExecuteDesign <- function(...){

  dots <- list(...)
  designTable <- makeDesign(indices = dots$indices, targetN = dots$targetN,
                            pmin = dots$pmin, pmax = dots$pmax,
                            pEscape = dots$pEscape)

  cD <- executeDesign(L = dots$L,
                      dT = designTable,
                      maxCells = dots$maxCells)

  return(cD)
}
