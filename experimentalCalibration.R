library(SpaDES)
library(LandR)
library(raster)
library(sf)
library(rgeos)
library(magrittr)
library(fasterize)
####Generate a random Polygon somewhere in the Boreal####
#this area varies in space, unlike randomStudyArea, which only changes shape
genRandomBorealArea <- function(size) {
  #Build extent (= to Boreal Plains)
  temp <- matrix(c(-606073.8, -606073.8,
                   1218651, 1218651,
                   5119032, 6379649,
                   6379649, 5119032), ncol = 2) %>%
    Polygon(.)

  temp <- Polygons(srl = list(temp), ID = data.frame("id" = 1))
  tempSP <- SpatialPolygons(Srl = list(temp))
  crs(tempSP) <- crs("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
+ellps=GRS80 +towgs84=0,0,0")
  tempPoint <- sp::spsample(tempSP, n = 1, type = "random") %>%
   SpatialPoints(.)
  crs(tempPoint) <- crs("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
+ellps=GRS80 +towgs84=0,0,0")

  out <- randomPolygon(x = tempPoint, area = size)
  out
}

####Generate the needed rasters####

genSimLand <- function(size = 25000000, buffDist = 5000){

  tempDir <- tempdir()

  #Generate study Area
  studyArea <- genRandomBorealArea(size)

  #Buffer study Area
  bStudyArea <- buffer(studyArea, buffDist) %>%
    gDifference(., spgeom2 = studyArea, byid = FALSE)
  polyLandscape <- sp::rbind.SpatialPolygons(studyArea, bStudyArea)
  polyLandscape$zone <- c("core", "buffer")
  polyLandscape$Value <- c(1, 0)

  #Generate flammability raster
  landscapeLCC <- prepInputsLCC(destinationPath = tempDir, studyArea = polyLandscape, useSAcrs = TRUE)
  landscapeFlam <- defineFlammable(landscapeLCC)
  #Generate landscape Index raster
  polySF <- sf::st_as_sf(polyLandscape)
  landscapeIndex <- fasterize(polySF, landscapeLCC, "Value")

  simLand <- list(polyLandscape, landscapeIndex, landscapeLCC, landscapeFlam)
  names(simLand) <- c("studyArea", "landscapeIndex", "lcc", "flammableMap")
  return(simLand)
  }

sim1 <- genSimLand()

clearPlot()
# Plot(sim1$studyArea)
# Plot(sim1$landscapeIndex)
# Plot(sim1$lcc)
plot(sim1$flammableMap)

#Need a vector of igniteable cells
#Item 1 = L, the flammable Map
#Item 2 = B (aka the landscape Index) this denotes buffer
#Item 3 = igLoc(index of igniteable cells) L[igloc] == 1 &&(B[igLoc]) == 1 (ie within core)
index <- 1:ncell(sim1$flammableMap)
index[sim1$flammableMap[] != 1 | is.na(sim1$flammableMap[])] <- NA
index[sim1$landscapeIndex[] != 1 | is.na(sim1$landscapeIndex[])] <- NA
index <- index[!is.na(index)]

#dT <- data.frame("igLoc" = index, p0 = 0.1, p = 0.23)

#this version of makeDesign is the simplest possible...

makeDesign <- function(indices, targetN=1000, pEscape=0.1, pmin=0.18, pmax=0.26, q=1){
  
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

executeDesign <- function(L, dT){

  # extract elements of dT into a three column matrix where column 1,2,3 = igLoc, p0, p
 #browser()
  f <- function(x, L, P){ #L, P are rasters, passed by reference
    #browser()
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
    P[nbrs] <- p0
    #Now it is 1, 0, p0, and NA
    spreadState0 <- SpaDES.tools::spread2(landscape = L,
                                          start = i,
                                          iterations = 1,
                                          spreadProb = P,
                                          asRaster = FALSE)

    tmp <- nrow(spreadState0)
    res[2:3] <- c(tmp-1,tmp)
    if (tmp==1) #the fire did not spread.
      return(res)
    P[] <- L[]*p
    spreadState1 <- SpaDES.tools::spread2(landscape = L,
                                          start = spreadState0,
                                          spreadProb = P,
                                          asRaster = FALSE)
    #calculate return data
    res[3] <- nrow(spreadState1)
    return(res)
  }

  P <- raster(L) # I think this makes a new raster
  P[] <- L[]
  #@Ian: need to ensure that L and P are passed by reference
  res <- apply(dT, 1, f, L, P) #f(T[i,], L, P)

  return(t(res)) #t() to conform with dT

}

dT = makeDesign(indices=index)
calibData <- executeDesign(L = sim1$flammableMap, dT)

