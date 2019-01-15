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
  scfmRegime = list(fireCause=c("L", "H"))
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
#tolkoABN <- shapefile("../LandWeb/inputs/FMA_Boundaries/Tolko/Tolko_AB_S_SR.shp")
# tolkoSK <- shapefile("../LandWeb/inputs/FMA_Boundaries/Tolko/Tolko_SK_SR.shp")
objects <- list(
# studyArea = shapefile("~/GitHub/LandWeb/inputs/FMA_Boundaries/Tolko/Tolko_AB_N_SR.shp")
  studyArea =  subEcoDistricts, #ABN  tolkoSK
  rasterToMatch = rasterToMatch
)

paths <- list(
  cachePath = file.path("cache"),
  modulePath = file.path("modules"),
  inputPath = inputDir,
  outputPath = outputDir
)

setPaths(cachePath = paths$cachePath,
         modulePath = paths$modulePath,
         inputPath = paths$inputPath,
         outputPath = paths$outputPath)

.Options$reproducible.cachePath <- paths$cachePath
options(spades.moduleCodeChecks = TRUE)
mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = paths)
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
