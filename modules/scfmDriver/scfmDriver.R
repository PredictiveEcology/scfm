defineModule(sim, list(
  name = "scfmDriver",
  description = "generate parameters for the generic percolation model",# spades::spread()",
  keywords = c("fire"),
  authors = c(person(c("Steve", "G"), "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.txt", "scfmDriver.Rmd"),
  reqdPkgs = list("stats", "magrittr", "sf", "rgeos", "fasterize"),
  parameters = rbind(
    defineParameter("neighbours", "numeric", 8, 4, 8, "number of cell immediate neighbours"),
    defineParameter("buffDist", "numeric", 5e3, 0, 1e5, "Buffer width for fire landscape calibration"),
    defineParameter("targetN", "numeric", 1000, 1, NA, "target sample size for determining true spread probability")),
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmRegimePars", objectClass = "list", desc = ""),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = ""),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame",
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
  #
  cellSize <- sim$landscapeAttr[[1]]$cellSize

  sim$scfmDriverPars <- lapply(names(sim$scfmRegimePars), function(polygonType, targetN = P(sim)$sampleSize) {
    regime <- sim$scfmRegimePars[[polygonType]]
    landAttr <- sim$landscapeAttr[[polygonType]]
    maxBurnCells <- as.integer(round(regime$emfs / cellSize))

    #we know this table was produced with MinFireSize=2cells.

    # y <- sim$cTable2$y #What are these supposed to be?
    # x <- sim$cTable2$p
    # m.lw <- lowess(y~x, iter = 2)
    # if (sum(diff(m.lw$y) < 0) > 0)
    #   warning(sprintf("lowess curve non-monotone in zone %s. Proceed with caution", polygonType))
    # targetSize <- regime$xBar/cellSize - 1
    # pJmp <- approx(m.lw$y, m.lw$x, targetSize, rule = 2)$y
    message("generating buffered landscapes...")
    calibLand <- genSimLand(sim$studyArea[polygonType,], buffDist = P(sim)$buffDist)

    #Need a vector of igniteable cells
    #Item 1 = L, the flammable Map
    #Item 2 = B (aka the landscape Index) this denotes buffer
    #Item 3 = igLoc(index of igniteable cells) L[igloc] == 1 &&(B[igLoc]) == 1 (ie within core)
    index <- 1:ncell(calibLand$flammableMap)
    index[calibLand$flammableMap[] != 1 | is.na(calibLand$flammableMap[])] <- NA
    index[calibLand$landscapeIndex[] != 1 | is.na(calibLand$landscapeIndex[])] <- NA
    index <- index[!is.na(index)]
    #index is the set of locations where fires may Ignite.

    dT = makeDesign(indices=index, targetN = P(sim)$targetN, pEscape=regime$pEscape)
    message(paste0("calibrating for polygon ", polygonType))
    calibData <- Cache(executeDesign, L = calibLand$flammableMap, dT,
                       maxCells=maxBurnCells,
                       userTags = paste("executeDesign", polygonType))

    cD <- calibData[calibData$finalSize > 1,]  #could use [] notation, of course.
    calibModel <- loess(cD$finalSize ~ cD$p)
    #now for the inverse step.
    xBar <- regime$xBar / cellSize

    pJmp <- try(stats::uniroot(f <- function(x, cM, xBar) {predict(cM,x) - xBar},
                    calibModel, xBar, # "..."
                    interval=c(min(cD$p), max(cD$p)),
                    extendInt = "no",
                    tol = 0.00001
                    ), silent = TRUE)
    if (class(pJmp) == "try-error") {
      pJmp <- min(calibModel$x)
      message("the loess model may underestimate the spread probability for polygon ", polygonType)
    }
    #check convergence, and out of bounds errors etc
    w <- landAttr$nNbrs
    w <- w/sum(w)
    hatPE <- regime$pEscape
    if (hatPE == 0) {
      # no fires in polygon zone escapted
      p0 <- 0
    } else if (hatPE == 1) {
      # all fires in polygon zone escaped
      p0 <- 1
    } else {
      res <- optimise(escapeProbDelta,
                      interval = c(hatP0(hatPE, P(sim)$neighbours),
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
                maxBurnCells = maxBurnCells
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


genSimLand <- function(coreLand, buffDist){

  tempDir <- tempdir()
  #Buffer study Area
  bStudyArea <- buffer(coreLand, buffDist) %>%
    gDifference(., spgeom2 = coreLand, byid = FALSE)
  polyLandscape <- sp::rbind.SpatialPolygons(coreLand, bStudyArea)
  polyLandscape$zone <- c("core", "buffer")
  polyLandscape$Value <- c(1, 0)

  #Generate flammability raster
  landscapeLCC <- prepInputsLCC(destinationPath = tempDir, studyArea = polyLandscape, useSAcrs = TRUE)
  landscapeFlam <- defineFlammable(landscapeLCC)
  #Generate landscape Index raster
  polySF <- sf::st_as_sf(polyLandscape)
  landscapeIndex <- fasterize(polySF, landscapeLCC, "Value")

  calibrationLandscape <- list(polyLandscape, landscapeIndex, landscapeLCC, landscapeFlam)
  names(calibrationLandscape) <- c("studyArea", "landscapeIndex", "lcc", "flammableMap")
  return(calibrationLandscape)
}


#dT <- data.frame("igLoc" = index, p0 = 0.1, p = 0.23)

#this version of makeDesign is the simplest possible...

makeDesign <- function(indices, targetN, pEscape=0.1, pmin=0.18, pmax=0.26, q=1){
  sampleSize <- round(targetN/pEscape)
  cellSample <- sample(indices, sampleSize, replace = TRUE)
  pVec <- runif(sampleSize)^q
  pVec <- pVec * (pmax-pmin) + pmin

  #derive p0 from escapeProb
  #steal code from scfmRegime and friends.

  p0 <- 1 - (1 - pEscape)^0.125  #assume 8 neighbours
  #the preceding approximation seems inadequate in practice.
  #when implemented in scfmDriver, make use of correct derivation of p0 from pEscape based on L
  T <- data.frame("igLoc" = cellSample, "p0" = p0, "p" = pVec)
  return(T)
}

executeDesign <- function(L, dT, maxCells){

  # extract elements of dT into a three column matrix where column 1,2,3 = igLoc, p0, p

  f <- function(x, L, ProbRas){ #L, P are rasters, passed by reference

    i <- x[1]
    p0 <- x[2]
    p <-x[3]

    nbrs <- raster::adjacent(L, i, pairs=FALSE, directions=8)
    #nbrs < nbrs[which(L[nbrs]==1)] #or this?
    nbrs <- nbrs[L[nbrs]==1] #only flammable neighbours please. also, verify NAs excluded.
    #nbrs is a vector of flammable neighbours.
    nn <- length(nbrs)
    res = c(nn,0,1)
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
    res[2:3] <- c(tmp-1,tmp)
    if (tmp==1) #the fire did not spread.
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


  res <- Cache(apply, dT, 1, f, L, ProbRas = probRas) #f(T[i,], L, P)

  res <- data.frame("nNeighbours" = res[1,], "initSpreadEvents" = res[2,], "finalSize" = res[3,])

  #cbind dT and res, then select the columns we need
  x <- cbind(dT,res)

  return(x)
}



