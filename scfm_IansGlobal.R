library(magrittr)
library(raster)
library(SpaDES)

inputDir <- file.path("inputs")
outputDir <- file.path("outputs")

timeunit <- "year"
times <- list(start = 0, end = 30)
mapDim <- 200
defaultInterval <- 1.0
defaultPlotInterval <- 1.0
defaultInitialSaveTime <- NA #don't be saving nuffink

parameters <- list(
  .progress = list(type = "text", interval = 1),
  ageModule = list(
    initialAge = 100,
    maxAge = 200,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = times$start,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmIgnition = list(
    pIgnition = 0.0001,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = NA,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmEscape = list(
    p0 = 0.05,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = NA,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmSpread = list(
    pSpread = 0.235,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = times$start,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmRegime = list(fireCause=c("L", "H")),
  scfmDriver = list(targetN = 1500)
  # andisonDriver =   list(pSpreadOddsRatio = 1,#1.025,
  #                        mfsMaxRatio = 3,
  #                        mfsMultiplier = 3.25),
  # andisonDriver_dataPrep = list(minFRI=0)
)

modules <- list("scfmLandcoverInit","scfmIgnition","scfmDriver",
                "ageModule", "scfmRegime", "scfmEscape", "scfmSpread")
#, "andisonDriver_dataPrep","andisonDriver"

# AndisonFRI <- shapefile("modules/andisonDriver/data/landweb_ltfc_v6.shp")
# AndisonFRI <- raster::aggregate(AndisonFRI[AndisonFRI$LTHFC > 40,],
#                                by = 'LTHFC', dissolve = TRUE)

ecoDistricts <- shapefile("C:/Ian/Data/Canada Ecosystem/ecodistrict_shp/Ecodistricts/ecodistricts.shp")
subEcoDistricts <- ecoDistricts[ecoDistricts$DISTRICT_I %in% c(387, 390, 372),]
#Three small relatively contiguous ecodistricts in the RIA

# Generate rasterToMatch WHICH IS NOW REQUIRED!
# rasterToMatch <- Cache(pemisc::prepInputsLCC,
#                        studyArea = subEcoDistricts, year = 2005,
#                        filename2 = "C:/Ian/PracticeDirectory/scfm/subEcoDistricts.tif", overwrite = TRUE)

rasterToMatch <- raster("C:/Ian/PracticeDirectory/scfm/subEcoDistricts.tif")
subEcoDistricts <- spTransform(x = subEcoDistricts,
                               CRSobj = crs(rasterToMatch))


objects <- list(
  studyArea =  subEcoDistricts,
  rasterToMatch = rasterToMatch
)

setPaths(cachePath = file.path("cache"),
         modulePath = file.path("modules"),
         inputPath = inputDir,
         outputPath = outputDir)

options(spades.moduleCodeChecks = TRUE)
mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = getPaths())
dev()
set.seed(23657)
outSim <- SpaDES.core::spades(mySim, progress = FALSE, debug = TRUE)

###############
###############
report<-function(outSim){
 fri <<- unlist(lapply(outSim$scfmDriverPars,function(x)x$fri))
 area <<- unlist(lapply(outSim$scfmDriverPars,function(x)x$burnyArea))
 meanFri <<- sum(fri * (area/sum(area)))
 burnable <<- unlist(lapply(outSim$landscapeAttr,function(x)length(x$cellsByZone)))
 burned <<- unlist(lapply(outSim$landscapeAttr,function(x,y=outSim$burnMap)sum(y[x$cellsByZone])))
 zonalFri <<- 1/((burned/burnable)/end(outSim))
 return(NULL)
}


#####Test Study Areas
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
#build SA
set.seed(23657)
randomArea1 <- genRandomBorealArea(size = 1e8*8)
randomArea1$polyID <- 1

randomRTM1 <- LandR::prepInputsLCC(studyArea = randomArea1, useSAcrs = TRUE, destinationPath = tempdir())
#set SA
objects <- list(
  studyArea =  randomArea1,
  rasterToMatch = randomRTM1
)
#run SA
mySimR1 <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = getPaths())
dev()
outSimR1 <- SpaDES.core::spades(mySim, progress = FALSE, debug = TRUE)

setseed(2300)
randomArea2 <- genRandomBorealArea(size = 1e7*5)
randomArea2$polyID <- 2
randomRTM2 <- LandR::prepInputsLCC(studyArea = randomArea2, useSAcrs = TRUE, destinationPath = tempdir())
