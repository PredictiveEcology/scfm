library(SpaDES)
library(magrittr)
savedir<-file.path("outputs")
inputDir <- file.path("inputs") 
outputDir <- file.path("outputs")

timeunit<-"year"
times <- list(start = 0, end = 10)
mapDim <- 200
defaultInterval <- 1.0
defaultPlotInterval <- 1.0
defaultInitialSaveTime <- NA #don't be saving nuffink
parameters <- list(
  .progress = list(type = "text", interval = 1),
  ageModule = list(
    initialAge=100, 
    maxAge=200, 
    returnInterval = defaultInterval, 
    startTime = times$start,
    .plotInitialTime = times$start,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime, 
    .saveInterval=defaultInterval),
  scfmIgnition = list(
    pIgnition=0.0001,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = NA,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmEscape = list(
    p0=0.05,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = NA,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmSpread = list(
    pSpread=0.235,
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = times$start,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval)
)

modules <- list("scfmLandcoverInit", "scfmIgnition","ageModule","scfmRegime", "scfmEscape", "scfmSpread", "scfmDriver")
#
objects <- list(mapDim = mapDim) #note that these definitions are critical

paths <- list(
  cachePath = file.path(outputDir, "cache"),
  modulePath = file.path("modules"),
  inputPath = inputDir,
  outputPath = outputDir
)

setPaths(cachePath = paths$cachePath,
         modulePath = paths$modulePath,
         inputPath = paths$inputPath,
         outputPath = paths$outputPath)

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = paths)

set.seed(2343)
outSim <- spades(mySim, progress = FALSE)


#### Test with multiple study Areas ####

require(raster)
ecoDistricts <- shapefile(file.path(paths$inputPath, "Ecodistricts", "ecodistricts.shp"))
subEcoDistricts <- ecoDistricts[ecoDistricts$DISTRICT_I %in% c(348,350,347),] #Three semi-adjacent ecoDistricts
objects <- list(studyArea = subEcoDistricts)
#test
mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = paths)

set.seed(2343)
outSim <- spades(mySim, progress = FALSE, debug = TRUE)
