###This is a script I am using to test calibration over multiple polygons## IE

library(magrittr)
library(raster)
library(SpaDES)
library(LandR)

timeunit <- "year"
times <- list(start = 0, end = 50)
defaultInterval <- 1.0
defaultPlotInterval <- 1.0
defaultInitialSaveTime <- NA

setPaths(cachePath = file.path("cache"),
         modulePath = file.path("modules"),
         inputPath = file.path("inputs"),
         outputPath = file.path("outputs"))
paths <- getPaths()

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

studyArea <- shapefile("../../data/scfmRIA/RIA_studyArea.shp")
rasterToMatch <- raster("../../data/scfmRIA/RIA_studyAreaRas.tif")

# studyArea <- shapefile("C:/Ian/PracticeDirectory/scfm/RIA_studyArea.shp")
# rasterToMatch <- raster("C:/Ian/PracticeDirectory/scfm/RIA_studyAreaRas.tif")
# studyArea <- spTransform(studyArea, crs(rasterToMatch))
cloudFolderID <- "https://drive.google.com/open?id=1PoEkOkg_ixnAdDqqTQcun77nUvkEHDc0"

objects <- list(
  "cloudFolderID" = cloudFolderID,
  "studyArea" =  studyArea,
  "rasterToMatch" = rasterToMatch
)

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


comparePredictions <- function(polyList, simList) {
  out <- lapply(polyList, FUN = function(x, sim = simList) {
    regime <- sim$scfmRegimePars[[x]]
    driver <- sim$scfmDriverPars[[x]]
    landscapeAttr <- sim$landscapeAttr[[x]]
    firePoints <- sim$firePoints[sim$firePoints$PolyID == x,]
    hist_MAAB <- sum(firePoints$SIZE_HA[firePoints$SIZE_HA > landscapeAttr$cellSize])/
      (landscapeAttr$burnyArea*(sim@params$scfmRegime$fireEpoch[2] - sim@params$scfmRegime$fireEpoch[1] + 1)) * 100
    #This is a long way of saying, sum of fires/ (flammable landscape * fire epoch )
    #hist_mfs will be NaN if there were no fires larger than one pixel
    pSpread <- driver$pSpread
    burnSum <- sim$burnSummary[sim$burnSummary$polyID == x,]
    burnSum$N <- as.numeric(burnSum$N)
    burnSum$areaBurned <- as.numeric(burnSum$areaBurned)
    burnSum <- burnSum[burnSum$N > 1]
    mod_MAAB <- sum(burnSum$areaBurned)/(landscapeAttr$burnyArea * (times(sim)$end - times(sim)$start)) * 100

    pred <- data.frame("PolyId" = x, #Polygon ID
                       "prdMeanSize" = regime$xBar, #The predicted (empirical) mean size of fires
                       "modMeanSize" = mean(burnSum$areaBurned), #The modeled mean size of fires
                       "pSpread" = pSpread, # The spread probability estimated from the SCAM model
                       "hist_MAAB" = hist_MAAB,#The empirical mean annual area burned (from NFDB 1970-2000)
                       "mod_MAAB" = mod_MAAB) #The modelled mean annual area burned
    return(pred)
  })
  return(out)
}

df <- comparePredictions(names(outSimT$scfmDriverPars), outSimT) %>%
  rbindlist(.)
#Some useful plots
breaks = c(0.20, 0.21, 0.22, 0.23, 0.24, 0.25)
library(ggplot2)
ggplot(df, aes(x = prdMeanSize, y = modMeanSize, color = pSpread)) +
  geom_point(aes(prdMeanSize, modMeanSize)) +
  scale_color_gradient(low = "yellow", high = "red")+
  labs(y = "predicted size", x = "empirical size") +
  theme_minimal() +
  #geom_text(aes(label = PolyId), hjust = 0, vjust = 0)
  geom_abline(slope = 1)


ggplot(df, aes(x = hist_MAAB, y = mod_MAAB, col = pSpread)) +
  geom_point(aes(hist_MAAB, mod_MAAB)) +
  labs(y = "model mean annual area burned (%)", x = "empirical mean annual area burned (%)") +
  scale_color_gradient(low = "yellow", high = "red")+
  theme_minimal() +  #+
  geom_abline(slope = 1)
# xlim(c(0, 8000)) +
# ylim(c(0, 12000)) # +
