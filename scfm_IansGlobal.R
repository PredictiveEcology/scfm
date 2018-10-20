library(magrittr)
library(raster)
library(SpaDES)

inputDir <- file.path("inputs")
outputDir <- file.path("outputs")

timeunit <- "year"
times <- list(start = 0, end = 10)
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
    .saveInterval = defaultInterval)
)

modules <- list("andisonDriver_dataPrep", "andisonDriver", "scfmLandcoverInit",
                "scfmIgnition", "ageModule", "scfmRegime", "scfmEscape", "scfmSpread")

objects <- list(
  mapDim = mapDim, # note that these definitions are critical
  studyArea0 = shapefile("~/GitHub/LandWeb/inputs/FMA_Boundaries/Tolko/Tolko_AB_N_SR.shp")
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

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = paths)

set.seed(2343)
outSim <- spades(mySim, progress = FALSE, debug = TRUE)

#### Test with multiple study Areas ####

require(raster)
#Don't try this it won't work yet
# ecoDistricts <- shapefile("modules/scfmLandcoverInit/data/ecodistricts_shp/Ecodistricts/ecodistricts.shp")
# subEcoDistricts <- ecoDistricts[ecoDistricts$DISTRICT_I %in% c(348,350,347),] #Three semi-adjacent ecoDistricts
# objects <- list(studyArea = subEcoDistricts, mapDim = mapDim)


#testing steve's data. #There is a problem right now with input objects. This only works if studyArea is Dave's shapefile

AndisonFRI <- shapefile("modules/andisonDriver/data/landweb_ltfc_v6.shp")
AndisonFRI <- raster::aggregate(AndisonFRI[AndisonFRI$LTHFC > 40,],
                                by = 'LTHFC', dissolve = TRUE)
ecoDistricts <- shapefile("modules/scfmLandcoverInit/data/ecodistricts_shp/Ecodistricts/ecodistricts.shp")
subEcoDistricts <- ecoDistricts[ecoDistricts$ECOREGION %in% c(198),] #Three semi-adjacent ecoDistricts
plot(subEcoDistricts)

subEcoDistricts <- spTransform(subEcoDistricts, CRSobj = crs(AndisonFRI))
#This is producing NULL files, sometimes single lines. wtf
AndisonFRIm <- crop(AndisonFRI, subEcoDistricts)

plot(AndisonFRIm)
AndisonFRIm$PolyID <- row.names(AndisonFRIm)
objects <- list(studyArea = AndisonFRIm)
plot(AndisonFRIm)

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = paths)

set.seed(2343)
outSim <- spades(mySim, progress = FALSE, debug = TRUE)
