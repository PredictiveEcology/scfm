library(LandR)
library(raster)
library(sf)
library(rgeos)
library(magrittr)
library(SpaDES)
library(fasterize)
####Generate a random Polygon somewhere in the Boreal####
#this area varies in space, unlike randomStudyArea, which only changes shape
genRandomBorealArea <- function(size) {
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
  tempPoint <- sp::spsample(tempSP, n = 1, type = "random") %>%
   SpatialPoints(.)
  crs(tempPoint) <- crs("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
+ellps=GRS80 +towgs84=0,0,0")

  out <- randomPolygon(x = tempPoint, area = size)
  out
}

####Generate the needed rasters####
tempDir <- tempdir()

studyArea <- genRandomBorealArea(size = 62500 * 400)
bStudyArea <- buffer(studyArea, 5000) %>%
 gDifference(., spgeom2 = studyArea, byid = FALSE)

polyLandscape <- sp::rbind.SpatialPolygons(studyArea, bStudyArea)

polyLandscape$zone <- c("core", "buffer")
polyLandscape$Value <- c(1, 0)
landscapeLCC <- prepInputsLCC(destinationPath = tempDir, studyArea = polyLandscape, useSAcrs = TRUE)

landscapeFlam <- defineFlammable(landscapeLCC)

polySF <- sf::st_as_sf(polyLandscape)
#Generate buffer index raster

landscapeIndex <- fasterize(polySF, landscapeLCC, "Value")

# setColors(landscapeIndex, n = 2) <- colorRampPalette(c("purple", "yellow"))(2)
# getValues(landscapeIndex)
colfun <- colorRampPalette(colors = c("dark green", "green", "brown", "blue"))

out <- colfun(39)

Plot(polyLandscape)
Plot(landscapeIndex)
Plot(landscapeLCC, col = out, legend = FALSE)
Plot(landscapeFlam)
