defineModule(sim, list(
  name = "scfmDriver",
  description = "generate parameters for the generic percolation model",
  keywords = c("fire"),
  authors = c(person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre")),
              person("Ian", "Eddy", email = "ian.eddy@canada.ca", role = c("aut"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmDriver.Rmd"),
  reqdPkgs = list("fasterize", "PredictiveEcology/LandR@development", "magrittr", "parallel",
                  "PredictiveEcology/pemisc@development", "reproducible", "rgeos",
                  "scam", "sf", "sp", "SpaDES.tools", "stats"),
  parameters = rbind(
    defineParameter("neighbours", "numeric", 8, 4, 8, "number of cell immediate neighbours"),
    defineParameter("buffDist", "numeric", 5e3, 0, 1e5, "Buffer width for fire landscape calibration"),
    defineParameter("pJmp", "numeric", 0.23, 0.18, 0.25, "default spread prob for degenerate polygons"),
    defineParameter("targetN", "numeric", 1500, 1, NA, "target sample size for determining true spread probability")
  ),
  inputObjects = bind_rows(
    expectsInput("cloudFolderID", "character",
                 paste("URL for Google-drive-backed cloud cache. ",
                       "Note: turn cloudCache on or off with options('reproducible.useCloud')")),
    expectsInput("scfmRegimePars", "list", desc = ""),
    expectsInput("landscapeAttr", "list", desc = ""),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = "a studyArea where separate polygons denote separate fire regimes")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "scfmDriverPars", objectClass = "list", desc = "")
  )
))

## event types
#   - type `init` is required for initiliazation

doEvent.scfmDriver = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
    },

    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

# 1 - (1-p0)**N = pEscape
# 1 - pEscape = (1-p0)**N
# (1 - pEscape)**1/N = 1 - p0
# p0 = 1 - (1 - pEscape)**1/N

hatP0 <- function(pEscape, n = 8) {
  1 - (1 - pEscape) ** (1 / n)
}

#a real clever boots would minimise the abs log odds ratio.
#be my guest.

escapeProbDelta <- function(p0, w, hatPE) {
  abs(sum(w*(1 - (1 - p0) ** (0:8))) - hatPE)
}

Init <- function(sim) {
  cellSize <- sim$landscapeAttr[[1]]$cellSize

  if (getOption("pemisc.useParallel", FALSE)) {
    #options("pemisc.useParallel" = TRUE)
    cl <- pemisc::makeOptimalCluster(MBper = 5000, maxNumClusters = length(sim$scfmRegimePars))
    on.exit(try(stopCluster(cl), silent = TRUE))
  } else {
    cl <- NULL
  }

  # Eliot modified this to use cloudCache -- need all arguments named, so Cache works
  sim$scfmDriverPars <- cloudCache(
    pemisc::Map2, cl = cl, cloudFolderID = sim$cloudFolderID,
    useCache = getOption("reproducible.useCache", TRUE),
    useCloud = getOption("reproducible.useCloud", FALSE),
    regime = sim$scfmRegimePars, #[[polygonType]]
    landAttr = sim$landscapeAttr, #[[polygonType]]
    omitArgs = c("useCloud", "useCache", "cloudFolderID", "cl"),
    MoreArgs = list(cellSize = cellSize,
                    studyArea = sim$studyArea,
                    buffDist = P(sim)$buffDist,
                    pJmp = P(sim)$pJmp,
                    neighbours = P(sim)$neighbours),
    polygonType = names(sim$scfmRegimePars),
    function(polygonType, targetN = P(sim)$targetN,
             regime = regime, landAttr = landAttr,
             cellSize = cellSize,
             studyArea = sim$studyArea,
             buffDist = buffDist,
             pJmp = pJmp, neighbours = neighbours) {
      #regime <- sim$scfmRegimePars[[polygonType]] #pass as argument
      #landAttr <- sim$landscapeAttr[[polygonType]] #pass as argument
      maxBurnCells <- as.integer(round(regime$emfs_ha / cellSize)) #will return NA if emfs is NA
      if (is.na(maxBurnCells)) {
        warning("This can't happen")
        maxBurnCells = 1
      }

      message("generating buffered landscapes...")
      calibLand <- Cache(genSimLand, studyArea[polygonType,], buffDist = buffDist,
                         userTags = paste("genSimLand ", polygonType))

      #Need a vector of igniteable cells
      #Item 1 = L, the flammable Map
      #Item 2 = B (aka the landscape Index) this denotes buffer
      #Item 3 = igLoc(index of igniteable cells) L[igloc] == 1 &&(B[igLoc]) == 1 (ie within core)
      index <- 1:ncell(calibLand$flammableMap)
      index[calibLand$flammableMap[] != 1 | is.na(calibLand$flammableMap[])] <- NA
      index[calibLand$landscapeIndex[] != 1 | is.na(calibLand$landscapeIndex[])] <- NA
      index <- index[!is.na(index)]
      if (length(index) == 0)
        stop("polygon has no flammable cells!")

      #index is the set of locations where fires may Ignite.

      dT <- Cache(makeDesign, indices = index, targetN = targetN,
                  pEscape = ifelse(regime$pEscape == 0, 0.1, regime$pEscape),
                  userTags = paste("makeDesign", polygonType))

      message(paste0("calibrating for polygon ", polygonType, " (Time: ", Sys.time(), ")"))

      calibData <- Cache(executeDesign,
                         L = calibLand$flammableMap,
                         dT,
                         maxCells = maxBurnCells,
                         userTags = paste("executeDesign", polygonType)
      )

      cD <- calibData[calibData$finalSize > 1,]  #could use [] notation, of course.
      #calibModel <- loess(cD$finalSize ~ cD$p)
      calibModel <- scam::scam(finalSize ~ s(p, bs = "micx", k = 20), data = cD)

      xBar <- regime$xBar / cellSize

      if (xBar > 0) {
        #now for the inverse step.
        Res <- try(stats::uniroot(f <- function(x, cM, xBar) {predict(cM, list("p" = x)) - xBar},
                                  calibModel, xBar, # "..."
                                  interval = c(min(cD$p), max(cD$p)),
                                  extendInt = "no",
                                  tol = 0.00001
        ), silent = TRUE)
        if (class(Res) == "try-error") {
          pJmp <- min(cD$p)
          message("the loess model may underestimate the spread probability for polygon ", polygonType)
        } else {
          pJmp <- Res$root
        }
      } else {
        #pJmp <- P(sim)$pJmp
        calibModel <- "No Model"
        Res <- "No Uniroot result"
      }
      #check convergence, and out of bounds errors etc
      w <- landAttr$nNbrs
      w <- w / sum(w)
      hatPE <- regime$pEscape
      if (hatPE == 0) {
        # no fires in polygon zone escapted
        p0 <- 0
      } else if (hatPE == 1) {
        # all fires in polygon zone escaped
        p0 <- 1
      } else {
        res <- optimise(escapeProbDelta,
                        interval = c(hatP0(hatPE, neighbours),
                                     hatP0(hatPE, floor(sum(w * 0:8)))),
                        tol = 1e-4,
                        w = w,
                        hatPE = hatPE)
        p0 <- res[["minimum"]]
        #It is almost obvious that the true minimum must occurr within the interval specified in the
        #call to optimise, but I have not proved it, nor am I certain that the function being minimised is
        #monotone.
      }
      #don't forget to scale by number of years, as well, if your timestep is ever != 1yr
      rate <- regime$ignitionRate * cellSize #fireRegimeModel and this module must agree on
      #an annual time step. How to test / enforce?
      pIgnition <- rate #approximate Poisson arrivals as a Bernoulli process at cell level.
      #for Poisson rate << 1, the expected values are the same, partially accounting
      #for multiple arrivals within years. Formerly, I used a poorer approximation
      #where 1-p = P[x==0 | lambda=rate] (Armstrong and Cumming 2003).

      return(list(pSpread = pJmp,
                  p0 = p0,
                  naiveP0 = hatP0(regime$pEscape, 8),
                  pIgnition = pIgnition,
                  maxBurnCells = maxBurnCells,
                  calibModel = calibModel,
                  uniroot.Res = Res
      )
      )
    })

  names(sim$scfmDriverPars) <- names(sim$scfmRegimePars) #replicate the polygon labels

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using random polygon in Alberta")
    #TODO: remove LandR once this is confirmed working
    studyArea <- LandR::randomStudyArea(size = 1e4*1e6, seed = 23654) #10,000 km * 1000^2m^2
    sim$studyArea <- studyArea
  }
  return(invisible(sim))
}

#Buffers polygon, generates index raster
genSimLand <- function(coreLand, buffDist) {
  tempDir <- tempdir()
  #Buffer study Area. #rbind had occasional errors before makeUniqueIDs = TRUE
  #TODO: Investigate why some polygons fail
  bStudyArea <- buffer(coreLand, buffDist) %>%
    rgeos::gDifference(., spgeom2 = coreLand, byid = FALSE)
  polyLandscape <- sp::rbind.SpatialPolygons(coreLand, bStudyArea, makeUniqueIDs = TRUE) #
  polyLandscape$zone <- c("core", "buffer")
  polyLandscape$Value <- c(1, 0)

  #Generate flammability raster
  landscapeLCC <- prepInputsLCC(destinationPath = tempDir, studyArea = polyLandscape, useSAcrs = TRUE)
  landscapeFlam <- defineFlammable(landscapeLCC)
  #Generate landscape Index raster
  polySF <- sf::st_as_sf(polyLandscape)
  landscapeIndex <- fasterize::fasterize(polySF, landscapeLCC, "Value")

  calibrationLandscape <- list(polyLandscape, landscapeIndex, landscapeLCC, landscapeFlam)
  names(calibrationLandscape) <- c("studyArea", "landscapeIndex", "lcc", "flammableMap")
  return(calibrationLandscape)
}

#dT <- data.frame("igLoc" = index, p0 = 0.1, p = 0.23)

#this version of makeDesign is the simplest possible...

makeDesign <- function(indices, targetN, pEscape = 0.1, pmin = 0.21, pmax = 0.2525, q = 1) {
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

  return(x)
}
