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
  andisonDriver =   list(pSpreadOddsRatio = 1,#1.025,
                         mfsMultiplier=1.55)
)

modules <- list("andisonDriver_dataPrep", "andisonDriver", "scfmLandcoverInit",
                "scfmIgnition", "ageModule", "scfmRegime", "scfmEscape", "scfmSpread")

#AndisonFRI <- shapefile("modules/andisonDriver/data/landweb_ltfc_v6.shp")
#AndisonFRI <- raster::aggregate(AndisonFRI[AndisonFRI$LTHFC > 40,],
#                                by = 'LTHFC', dissolve = TRUE)

#ecoDistricts <- shapefile("modules/scfmLandcoverInit/data/ecodistricts_shp/Ecodistricts/ecodistricts.shp")
#subEcoDistricts <- ecoDistricts[ecoDistricts$ECOREGION %in% c(87, 142),] #Two small relatively contiguous ecoregions

tolkoABN <- shapefile("../LandWeb/inputs/FMA_Boundaries/Tolko/Tolko_AB_N_SR.shp")

objects <- list(
# studyArea0 = shapefile("~/GitHub/LandWeb/inputs/FMA_Boundaries/Tolko/Tolko_AB_N_SR.shp")
  studyArea0 = tolkoABN # subEcoDistricts
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
options(spades.moduleCodeChecks = FALSE)


mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = paths)

set.seed(23435)
#outSim <- spades(mySim, progress = FALSE, debug = TRUE)

