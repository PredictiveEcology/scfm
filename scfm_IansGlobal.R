###This is a script I am using to test calibration over multiple polygons## IE

library(magrittr)
library(raster)
library(SpaDES)
library(LandR)

inputDir <- file.path("inputs")
outputDir <- file.path("outputs")

timeunit <- "year"
times <- list(start = 0, end = 100)
defaultInterval <- 1.0
defaultPlotInterval <- 1.0
defaultInitialSaveTime <- NA

parameters <- list(
  .progress = list(type = "text", interval = 1),
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
    returnInterval = defaultInterval,
    startTime = times$start,
    .plotInitialTime = times$start,
    .plotInterval = defaultPlotInterval,
    .saveInitialTime = defaultInitialSaveTime,
    .saveInterval = defaultInterval),
  scfmRegime = list(fireCause=c("L", "H"))
)

modules <- list("scfmLandcoverInit","scfmIgnition","scfmDriver",
                "ageModule", "scfmRegime", "scfmEscape", "scfmSpread")

opts <- options(
  "LandR.verbose" = 1,
  "map.dataPath" = paths$inputPath,
  "map.overwrite" = TRUE,
  "map.useParallel" = TRUE, #!identical("windows", .Platform$OS.type),
  "reproducible.futurePlan" = FALSE,
  "reproducible.inputPaths" = if (pemisc::user("emcintir")) "~/data" else NULL,
  "reproducible.quick" = FALSE,
  "reproducible.overwrite" = TRUE,
  "reproducible.useMemoise" = TRUE, # Brings cached stuff to memory during the second run
  "reproducible.useNewDigestAlgorithm" = TRUE,  # use the new less strict hashing algo
  "reproducible.useCache" = TRUE,
  "reproducible.cachePath" = paths$cachePath,
  "reproducible.showSimilar" = FALSE,
  "spades.moduleCodeChecks" = FALSE, # Turn off all module's code checking
  "spades.useRequire" = FALSE, # assuming all pkgs installed correctly
  "pemisc.useParallel" = FALSE,
  "reproducible.useCloud" = FALSE
)



#### RIA RUN ####
studyArea <- shapefile("C:/Ian/PracticeDirectory/scfm/RIA_studyArea.shp")
rasterToMatch <- raster("C:/Ian/PracticeDirectory/scfm/RIA_studyAreaRas.tif")
studyArea <- spTransform(studyArea, crs(rasterToMatch))
cloudFolderID <- "https://drive.google.com/open?id=1PoEkOkg_ixnAdDqqTQcun77nUvkEHDc0"
objects <- list(
  "cloudFolderID" = cloudFolderID,
  "studyArea" =  studyArea,
  "rasterToMatch" = rasterToMatch
)

setPaths(cachePath = file.path("cache"),
         modulePath = file.path("modules"),
         inputPath = inputDir,
         outputPath = outputDir)

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects, paths = getPaths())

set.seed(23658)

outSim <- SpaDES.core::spades(mySim, progress = FALSE, debug = TRUE)

# ####Tati's run#######
wd <- tempdir()
cloudFolderID <- "https://drive.google.com/open?id=1PoEkOkg_ixnAdDqqTQcun77nUvkEHDc0"
NWT.url <- "https://drive.google.com/open?id=1LUxoY2-pgkCmmNH5goagBp3IMpj6YrdU"
studyAreaNWT <- cloudCache(prepInputs,
                           url = NWT.url,
                           destinationPath = wd,
                           cloudFolderID = cloudFolderID,
                           omitArgs = c("destinationPath", "cloudFolderID"))

rasterToMatchNWT <- cloudCache(prepInputs, url = "https://drive.google.com/open?id=1fo08FMACr_aTV03lteQ7KsaoN9xGx1Df",
                            studyArea = studyAreaNWT,
                            targetFile = "RTM.tif", destinationPath = wd,
                            useCloud = TRUE,
                            cloudFolderID = cloudFolderID, overwrite = TRUE, filename2 = NULL,
                            omitArgs = c("destinationPath", "cloudFolderID", "useCloud", "overwrite", "filename2"))

newObjs <- list(
  "cloudFolderID" = cloudFolderID,
  "studyArea" =  studyAreaNWT,
  "rasterToMatch" = rasterToMatchNWT
)

mySimT <- simInit(times = times, params = parameters, modules = modules,
                      objects = newObjs, paths = getPaths())
outSimT <- SpaDES.core::spades(mySimT, progress = FALSE, debug = TRUE)

