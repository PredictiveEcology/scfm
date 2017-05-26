library(raster)
library(rgeos)
#### 
# This will load in several maps that are used 

####

#CanadaMap <- readRDS(file.path(paths$inputPath, "CanadaMap.rds"))
landisInputs <- readRDS(file.path(paths$inputPath, "landisInputs.rds"))
spEcoReg <- readRDS(file.path(paths$inputPath, "SpEcoReg.rds"))

seralStageData <- readRDS(file.path(paths$inputPath, "seralStageData.rds"))
vegTypeData <- readRDS(file.path(paths$inputPath, "vegTypeData.rds"))
availableRegions <- unique(vegTypeData$ecoregion)

# Study area original shapefile
loadShpAndMakeValid <- function(file) {
  shapefile(file) %>% gBuffer(byid=TRUE, width=0)
}

shpStudyRegionFull <- SpaDES::Cache(loadShpAndMakeValid, 
                                    file=file.path(paths$inputPath,"shpLandWEB.shp"),
              cacheRepo = paths$cachePath)
shpStudyRegionFull$fireReturnInterval <- shpStudyRegionFull$LTHRC
shpStudyRegionFull@data <- shpStudyRegionFull@data[,!(names(shpStudyRegionFull) %in% "ECODISTRIC")]

crsKNNMaps <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

seedState <- .Random.seed

set.seed(853839)#set.seed(5567913)
if (studyArea != "FULL") {
  if (studyArea == "SMALL") {
    areaKm2 <- 10000#700000#2000#600000#too big for laptop
  } else if (studyArea == "VERYSMALL") {
    areaKm2 <- 1000 #700000#2000#600000#too big for laptop
  } else if (studyArea == "MEDIUM") {
    areaKm2 <- 40000 #700000#2000#600000#too big for laptop
  } else if (studyArea == "LARGE") {
    areaKm2 <- 80000 #700000#2000#600000#too big for laptop
  } 
  
  minY <- 7778877 - 1.6e5
  shpStudyRegionFull <- spTransform(shpStudyRegionFull, crsKNNMaps)
  minX <- -1202250.2
  maxX <- minX + sqrt(areaKm2 * 1e6)
  maxY <- minY + sqrt(areaKm2 * 1e6)
  meanY <- mean(c(minY, maxY))
  
  # Add random noise to polygon
  xAdd <- -3e5#round(runif(1,-5e5, 1.5e6))
  yAdd <- 5e5#round(runif(1, 1e5, 5e5)) - xAdd / 2
  nPoints <- 20
  betaPar <- 0.6
  X <- c(jitter(sort(rbeta(nPoints, betaPar, betaPar) * (maxX - minX) + minX)),
         jitter(sort(rbeta(nPoints, betaPar, betaPar) * (maxX - minX) + minX, decreasing = TRUE)))
  Y <- c(jitter(sort(rbeta(nPoints / 2, betaPar, betaPar) * (maxY - meanY) + meanY)),
         jitter(sort(rbeta(nPoints, betaPar, betaPar) * (maxY - minY) + minY, decreasing = TRUE)),
         jitter(sort(rbeta(nPoints / 2, betaPar, betaPar) * (meanY - minY) + minY)))
  
  Sr1 <- Polygon(cbind(X + xAdd, Y + yAdd))
  Srs1 <- Polygons(list(Sr1), "s1")
  inputMapPolygon <- SpatialPolygons(list(Srs1), 1L)
  crs(inputMapPolygon) <- crsKNNMaps
  shpStudyRegion <- raster::intersect(shpStudyRegionFull, inputMapPolygon)
  
} else {
  shpStudyRegion <- shpStudyRegionFull
}

#ggStudyRegion <- ggvisFireReturnInterval(shpStudyRegion, shpStudyRegionFull)

set.seed(seedState)


