###This is a script I am using to test calibration over multiple polygons## IE

library(magrittr)
library(raster)
library(SpaDES)

inputDir <- file.path("inputs")
outputDir <- file.path("outputs")

timeunit <- "year"
times <- list(start = 0, end = 100)
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


###STUDY AREAS###
#The OG####
# ecoDistricts <- shapefile("C:/Ian/Data/Canada Ecosystem/ecodistrict_shp/Ecodistricts/ecodistricts.shp")
# #Three small relatively contiguous ecodistricts in the RIA
# subEcoDistricts <- ecoDistricts[ecoDistricts$DISTRICT_I %in% c(387, 390, 372),]
# subEcoDistricts <- spTransform(x = subEcoDistricts,
#                                CRSobj = crs(rasterToMatch))

#The new BCR6 NWT Ecoregion amalgamation
# ecoregions <- shapefile("C:/Ian/Data/Canada Ecosystem/ecoregion_shp/Ecoregions/ecoregions.shp")
# BCR <- shapefile("C:/Ian/Data/BirdConservationRegions/BCR_Terrestrial_master.shp")
# BCR <- BCR[BCR$BCR == 6,]
# NWT <- shapefile("C:/Ian/Data/Can_pol_boundaries/provinces/Provincial_Boundaries_OpenCanada/gpr_000b11a_e.shp")
# NWT <- NWT[NWT$PRUID == 61,]
#
# SAmask <- rgeos::gIntersection(BCR, NWT)
# ecoregions <- spTransform(ecoregions, CRSobj = crs(NWT))
# studyArea <- crop(ecoregions, SAmask)
#
#
# # Generate rasterToMatch WHICH IS NOW REQUIRED!
rasterToMatch <- Cache(LandR::prepInputsLCC, studyArea = studyArea, year = 2005,
                       destinationPath = "C:/Ian/PracticeDirectory/scfm",
                       filename2 = TRUE, overwrite = TRUE)
# studyArea <- spTransform(x = studyArea, CRSobj = crs(rasterToMatch))
# writeOGR(obj = studyArea, dsn = "C:/Ian/PracticeDirectory/scfm", layer = "BCR6_EcoregionsNWT", driver = "ESRI Shapefile")
studyArea <- shapefile("C:/Ian/PracticeDirectory/scfm/BCR6_EcoregionsNWT.shp")

rasterToMatch <- raster("C:/Ian/PracticeDirectory/scfm/SmallLCC2005_V1_4a.tif")

objects <- list(
  studyArea =  studyArea,
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
clearPlot()
set.seed(23658)

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
