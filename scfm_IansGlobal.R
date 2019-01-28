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
  scfmDriver = list(targetN = 1000)
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


.objects <- list(
  studyArea =  subEcoDistricts,
  rasterToMatch = rasterToMatch
)

setPaths(cachePath = file.path("cache"),
         modulePath = file.path("modules"),
         inputPath = inputDir,
         outputPath = outputDir)

options(spades.moduleCodeChecks = TRUE)

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = .objects, paths = getPaths())
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
genRandomBorealArea <- function(size, seed) {
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
  set.seed(seed)
  tempPoint <- sp::spsample(tempSP, n = 1, type = "random") %>%
    SpatialPoints(.)
  crs(tempPoint) <- crs("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
                        +ellps=GRS80 +towgs84=0,0,0")

  out <- randomPolygon(x = tempPoint, area = size)
  out
}

#RUN 1
#build SA

# randomArea1 <- genRandomBorealArea(size = 1e9*8)
# randomArea1$polyID <- 1
#
#run SA
# writeOGR(randomArea1, dsn = "C:/Ian/PracticeDirectory/scfm", layer = "randomArea1", driver = "ESRI Shapefile")
randomArea1 <- shapefile("C:/Ian/PracticeDirectory/scfm/randomArea1.shp")
# randomRTM1 <- LandR::prepInputsLCC(studyArea = randomArea1, useSAcrs = TRUE, destinationPath = tempdir())
# #set SA
# objects1 <- list(
#   studyArea =  randomArea1,
#   rasterToMatch = randomRTM1
# )

mySimR1 <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects1, paths = getPaths())
dev()
outSimR1 <- SpaDES.core::spades(mySimR1, progress = FALSE, debug = TRUE)

outSim$firePoints$IntensityCol <- "blue"
outSim$firePoints$IntensityCol[outSim$firePoints$SIZE_HA > 1] <- "yellow"
outSim$firePoints$IntensityCol[outSim$firePoints$SIZE_HA > 100] <- "orange"
outSim$firePoints$IntensityCol[outSim$firePoints$SIZE_HA > 1000] <- "red"
outSim$firePoints$IntensityCol[outSim$firePoints$SIZE_HA > 10000] <- "black"
temp <- outSim$firePoints[outSim$firePoints$SIZE_HA > 1,]
Plot(outSim$studyAreaRas)
Plot(temp, addTo = "outSim$studyAreaRas", cols = temp$IntensityCol)

#RUN2

# randomArea2 <- genRandomBorealArea(size = 1e9*3, seed = 2343)
# randomArea2$polyID <- 2
# writeOGR(randomArea2, dsn = "C:/Ian/PracticeDirectory/scfm", layer = "RandomArea2", driver = 'ESRI Shapefile')
randomArea2 <- shapefile("C:/Ian/PracticeDirectory/scfm/randomArea2.shp")
randomRTM2 <- LandR::prepInputsLCC(studyArea = randomArea2, useSAcrs = TRUE, destinationPath = tempdir())
plot(randomRTM2)
objects2 <- list(
  studyArea =  randomArea2,
  rasterToMatch = randomRTM2
)

#run SA
mySimR2 <- simInit(times = times, params = parameters, modules = modules,
                   objects = objects2, paths = getPaths())
dev()
outSimR2 <- SpaDES.core::spades(mySimR2, progress = FALSE, debug = TRUE)

