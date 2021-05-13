library(SpaDES.core)
library(reproducible)
library(raster)
library(data.table)
library(ggplot2)
library(LandR)
library(pemisc)

paths <- list(inputPath = 'inputs',
              cachePath = 'cache',
              outputPath = 'outputs',
              modulePath = 'modules')

do.call(setPaths, paths)

studyArea <- prepInputs(url = 'https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip',
                        destinationPath = 'inputs')

studyArea <- spTransform(studyArea, crs('+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'))

rasterToMatch <- LandR::prepInputsLCC(destinationPath = 'inputs')
rasterToMatch <- Cache(crop, rasterToMatch, studyArea, userTags = c("cropRTM"))
rasterToMatch <- Cache(mask, rasterToMatch, studyArea, userTags = c("maskRTM"))

fireRegimePolys <- studyArea

L70to2000parameters <- list(
  scfmLandcoverInit = list(
    removeSlivers = FALSE,
    sliverThreshold = 1e8), #polygons <100 km2 are merged with closest non-sliver neighbour
  scfmRegime = list(
    fireCause=c("L"),
    fireEpoch = c(1970-2000))
)


modules <- list('scfmLandcoverInit', 'scfmRegime')
flammableMap <- Cache(defineFlammable, LandCoverClassifiedMap = rasterToMatch, userTags = c("flammableMap"))

objects <- list(fireRegimePolys = fireRegimePolys,
                studyArea = studyArea,
                flammableMap = flammableMap,
                vegMap = rasterToMatch,
                rasterToMatch = rasterToMatch)

L70To2000 <- simInitAndSpades(modules = modules,
                              paths = paths,
                              times = list(start = 1, end = 2),
                              objects = objects,
                              params = L70To2000parameters)


A70to2000parameters <- list(
  scfmLandcoverInit = list(
    removeSlivers = FALSE,
    sliverThreshold = 1e8), #polygons <100 km2 are merged with closest non-sliver neighbour
  scfmRegime = list(
    fireCause = c("L", "H"),
    fireEpoch = c(1970, 2000))
)

A70To2000 <- simInitAndSpades(modules = modules,
                              paths = paths,
                              times = list(start = 1, end = 2),
                              objects = objects,
                              params = A70to2000parameters)

getStat <- function(regime, stat){ return(unlist(lapply(regime, FUN = function(x){x[stat]})))}
allCauseMax <- getStat(A70To2000$scfmRegimePars, stat = "emfs_ha")
LightningMax <- getStat(L70To2000$scfmRegimePars, stat = "emfs_ha")

allCauseMean <- getStat(A70To2000$scfmRegimePars, stat = "xBar")
LightningMean <- getStat(L70To2000$scfmRegimePars, stat = "xBar")
LightningPEscape <- getStat(L70To2000$scfmRegimePars, stat = "pEscape")
allCauseEscape <- getStat(A70To2000$scfmRegimePars, stat = "pEscape")
LightningIgnitionRate <- getStat(L70To2000$scfmRegimePars, stat = "ignitionRate")
allCauseIgnitionRate <- getStat(A70To2000$scfmRegimePars, stat = "ignitionRate")

lightningStats <- data.table(ecozone = names(L70To2000$scfmRegimePars),
                             maxFireSize = LightningMax,
                             meanFireSize = LightningMean,
                             ignitionRate = LightningIgnitionRate,
                             escapeProb = LightningPEscape,
                             ignitionSource = "lightning")
allStats <- data.table(ecozone = names(A70To2000$scfmRegimePars),
                       maxFireSize = allCauseMax,
                       meanFireSize = allCauseMean,
                       ignitionRate = allCauseIgnitionRate,
                       escapeProb = allCauseEscape,
                       ignitionSource = "human and lightning")

scfmRegimeDefaults <- rbind(lightningStats, allStats)
scfmRegimeDefaults[, era := c("1970 to 2000")]
EcozoneNames <- as.data.table(A70To2000$fireRegimePolys@data)
set(EcozoneNames, NULL, c("Legend_1", "Legend", "REGION", "Desc_", "Desc2",
                          "Label_1", "Nom", "Shape_Leng", "Shape_Area", "Label"), NULL)
EcozoneNames[, areaMha := trueArea/100/1e6]
EcozoneNames[, trueArea := NULL]

#whoops just noticed strait of georgia was included
scfmRegimeDefaults <- EcozoneNames[scfmRegimeDefaults, on = c('PolyID' = "ecozone")]
scfmRegimeDefaults <- scfmRegimeDefaults[!Name %in% c("Strait of Georgia", "Northern Arctic")]

#Some plots
firePlot2 <- ggplot(data = scfmRegimeDefaults, aes(y = maxFireSize, x = meanFireSize, col = Name)) +
  geom_point(aes(shape = ignitionSource), size = 3, alpha = 0.6) +
  scale_color_manual(values = pals::glasbey(n = length(unique(scfmRegimeDefaults$Name)))) +
  labs(x = "mean fire size (ha)", y = "max fire size (ha)") +
  scale_y_log +
  scale_x_log10() +
  guides(col = guide_legend("ecozone"),
         shape = guide_legend("ignition source"))
firePlot2
ggsave(firePlot2, filename = "outputs/meanAndMaxFireByEcozone.png", device = "png")
saveRDS(scfmRegimeDefaults, "inputs/scfmRegimeDefaults_250mRes.rds")
