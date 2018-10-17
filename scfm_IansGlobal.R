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
  .globals = list("neighbours"=8), # globals(sim)$neighbours
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

modules <- list("scfmIgnition")
#"scfmLandcoverInit", ,"ageModule","scfmRegime", "scfmEscape", "scfmSpread", "scfmDriver"
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
