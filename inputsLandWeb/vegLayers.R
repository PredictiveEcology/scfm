library(raster)
#library(SpaDES)
devtools::load_all("~/GitHub/SpaDES/.")
library(magrittr)
library(rgeos)
library(data.table)
library(gdalUtils)
gdal_setInstallation(ignore.full_scan = TRUE, 
                     verbose = TRUE)
startTime <- Sys.time()
dev()
setPaths()

# Step 1 -- Load LandWeb study area shapefile
loadShpAndMakeValid <- function(file) {
  shapefile(file) %>% gBuffer(byid=TRUE, width=0)
}
oldfilename <- "shpLandWEB.shp"

PaulRawFileName <- "C:/data/LandWeb/LANDWEB_DATA_LAYERS/SPP_1990_FILLED_100m_NAD83_LCC_BYTE_VEG.dat"
PaulOnGoogleDrive <- raster(PaulRawFileName)
LandWebStudyAreaRawPoly <- file.path("C:","Eliot","GitHub","LandWeb","inputs",oldfilename)


PaulOnGoogleDrive <- Cache(writeRaster, PaulOnGoogleDrive, filename = "PaulSppFilled.tif", datatype = "INT2U", 
                           overwrite = TRUE)
shpStudyRegionFull <- Cache(loadShpAndMakeValid, file=LandWebStudyAreaRawPoly)
studyAreaTooBig <- sp::spTransform(shpStudyRegionFull, crs(PaulOnGoogleDrive))
studyAreaTooBig$AREA <- round(studyAreaTooBig$AREA,0)


# Make Study Area Mask polygon
makeStudyAreaMask <- function(ras, maskFilename) {
  ras[] <- ras[]
  ras[ras[]%in% c(255)] <- NA_integer_ # out of study area
  rasterOptions(maxmemory=1e9)
  studyAreaMask <- raster(ras)
  studyAreaMask[] <- !is.na(ras[])
  Cache(writeRaster, studyAreaMask, filename = maskFilename, 
                         datatype = "INT1U", overwrite = TRUE)
}
studyAreaMask <- Cache(makeStudyAreaMask, PaulOnGoogleDrive, 
                       maskFilename ="StudyAreaMask.tif")
  
newfilename <- file.path(getwd(),"StudyAreaMask250.tif")
Cache(gdalwarp, overwrite=TRUE, 
      dstalpha = TRUE,
      s_srs= as.character(crs(studyAreaMask)),
      t_srs= as.character(crs(studyAreaMask)),
      multi=TRUE, of="GTiff", 
      tr=c(250, 250),
      filename(studyAreaMask), ot = "Byte", 
      newfilename)

Cache(system, paste("python",
                    file.path(getOption("gdalUtils_gdalPath")[[1]]$path,"gdal_polygonize.py"), 
             basename(newfilename), "shpLandWeb5.shp -f \"ESRI Shapefile\""))

shpPaul <- shapefile("shpLandWeb5.shp")#[2,]

#shpPaul <- shapefile("shpLandWeb3.shp")#[2,]
shpPaul <- shpPaul[shpPaul$DN==1,]
studyArea <- Cache(raster::intersect, studyAreaTooBig, shpPaul)
# cropMask <- function(ras, poly) {
#   message("starting crop")
#   PaulCropped <- crop(ras, poly)
#   message("starting cluster")
#   raster::beginCluster(8)
#   message("starting fastMask")
#   PaulTrimmed <- amc::fastMask(PaulCropped, poly)
#   message("starting endCluster")
#   raster::endCluster()
#   message("starting writeRaster")
#   PaulTrimmed <- writeRaster(PaulTrimmed, filename = "PaulTrimmed.tif", datatype="INT2U")
# }
# PaulTrimmed <- Cache(cropMask, PaulOnGoogleDrive, studyArea)

newfilename <- gsub(filename(PaulOnGoogleDrive), pattern = "\\.tif", replacement = "250.tif")
# Cache(gdalwarp, overwrite=TRUE, 
#       s_srs= as.character(crs(PaulOnGoogleDrive)),
#       t_srs= as.character(crs(PaulOnGoogleDrive)),
#       multi=TRUE, of="GTiff", tr=c(250, 250),  
#       filename(PaulOnGoogleDrive), ot = "Byte", 
#       newfilename)
# PaulOnGoogleDrive <- raster(newfilename)
#PaulOnGoogleDrive[] <- PaulOnGoogleDrive[]

Cache(shapefile, studyArea, filename = "studyArea.shp", overwrite = TRUE)
newfilename <- file.path(getwd(),"PaulTrimmed.tif")
Cache(gdalwarp, overwrite=TRUE, cutline = file.path(getwd(),"studyArea.shp"),
      dstalpha = TRUE,
      s_srs= as.character(crs(PaulOnGoogleDrive)),
      t_srs= as.character(crs(PaulOnGoogleDrive)),
      multi=TRUE, of="GTiff", 
      crop_to_cutline = TRUE, tr=c(250, 250),
      filename(PaulOnGoogleDrive), ot = "Byte", 
      newfilename)
PaulTrimmed <- raster(newfilename)
PaulTrimmed[] <- PaulTrimmed[]
PaulTrimmed[PaulTrimmed[]==255] <- NA_integer_
#gdalwarp -cutline INPUT.shp -crop_to_cutline -dstalpha INPUT.tif OUTPUT.tif



# make "small" number of pixels, centred in Alberta at -115
#studyArea <- spTransform(studyArea, CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-115 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#shapefile(studyArea, filename = "LandWebStudyArea.shp")



# studyAreaRas <- raster("c:/data/studyAreaRas.tif")
# 
# shpFilename <- file.path(getwd(),"shpLandWeb.shp")
# shpFilename <- file.path(getwd(),"shpLandWeb2.shp")
# rstFilename <- gsub(shpFilename, pattern = "\\.shp", replacement = ".tif")
# Cache(writeRaster, studyAreaMask, filename = rstFilename, datatype = "INT1U", overwrite = TRUE)
# #gdal_polygonizeR(rstFilename, readpoly=FALSE, outshape = shpFilename)
# 
# bb <- shapefile("shpLandWeb3.shp")
# 
# newfilename <- gsub(filename(studyAreaRas), pattern = "\\.tif", replacement = "Good.tif")
# Cache(gdalwarp, overwrite=TRUE, 
#       s_srs= as.character(crs(studyAreaRas)),
#       t_srs= as.character(crs(PaulOnGoogleDrive)),
#       multi=TRUE, of="GTiff", 
#       filename(studyAreaRas), 
#       newfilename)
# studyAreaRas <- raster(newfilename)
# studyAreaMask[studyAreaMask==0] <- NA_integer_
# studyAreaMask <- trim(studyAreaMask)
# sa <- crop(studyAreaRas,studyAreaMask)
# sa2 <- crop(studyAreaMask, sa)
# sa3 <- sa*sa2
# 
# raster::beginCluster(8)
# sa2 <- Cache(amc::fastMask, sa, bb) 
# raster::endCluster()
# 
# bbbb <- extract(studyAreaRas, studyAreaMask)
# 
# studyAreaMask[studyAreaMask==0] <- NA_integer_
# 
# a <- gdal_polygonizeR(studyAreaMask, readpoly=FALSE)


#PaulOnGoogleDrive[PaulOnGoogleDrive[]%in% c(230, 220)] <- NA_integer_ # water, non veg
#PaulOnGoogleDrive[PaulOnGoogleDrive[]==220] <- NA_integer_ # water

# studyAreaRas <- Cache(projectRaster, studyAreaRas, to = studyAreaMask)
# 
# StudyAreaPaul <- raster(PaulOnGoogleDrive)
# #StudyAreaPaul[] <- 
# studyAreaTooBig <- sp::spTransform(shpStudyRegionFull, crs(PaulOnGoogleDrive))
# studyAreaTooBig$AREA <- round(studyAreaTooBig$AREA,0)
# Cache(shapefile, studyAreaTooBig, filename = shpFilename, overwrite = TRUE)
# ccc <- intersect(studyAreaTooBig,bb)
# 
# 
# raster::beginCluster(8)
# rstStudyArea <- Cache(rasterize, studyAreaTooBig, PaulOnGoogleDrive)
# raster::endCluster()
# gdal_rasterize(a= "LTHRC", tr=c( 250.0, 250.0), l= "shpLandWeb", shpFilename, 
#                rstFilename, ot = "UInt16")
# studyAreaMask2 <- raster(rstFilename)



CASFRIRas <- raster("C:/data/LandWeb/Landweb_CASFRI_GIDs.tif")
# Step 1a -- warp CASFRI to 250m x 250m
# gdalwarp(overwrite=TRUE, 
#          t_srs= "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs",
#          multi=TRUE, of="GTiff", tr=c(250, 250),  
#          "C:/data/LandWeb/LANDWEB_DATA_LAYERS/SPP_1990_FILLED_100m_NAD83_LCC_BYTE_VEG.dat", 
# #          "C:/data/SppPaul250m.tif")
# newfilename <- gsub(filename(CASFRIRas), pattern = "\\.tif", replacement = "250.tif")
# Cache(gdalwarp, overwrite=TRUE, 
#       s_srs= as.character(crs(CASFRIRas)),
#       t_srs= as.character(crs(PaulOnGoogleDrive)),
#       multi=TRUE, of="GTiff", tr=c(250, 250),  
#       filename(CASFRIRas), 
#       newfilename)
# 
newfilename1 <- gsub(filename(CASFRIRas), pattern = "\\.tif", replacement = "StudyArea.tif")
Cache(gdalwarp, overwrite=TRUE, cutline = file.path(getwd(),"studyArea.shp"),
      s_srs= as.character(crs(CASFRIRas)),
      t_srs= as.character(crs(PaulOnGoogleDrive)),
      dstalpha = TRUE,tr=c(250, 250),  
      multi=TRUE, of="GTiff", 
      crop_to_cutline = TRUE,
      filename(CASFRIRas), #ot = "Byte", 
      newfilename1)

CASFRIRas <- raster(newfilename1)


# Step 2 -- Load CASFRI data and headers, and convert to long format, and define species groups
loadCASFRI <- function(CASFRIRas, attrFile, headerFile) {
  CASFRIattr <- fread(attrFile)
  CASFRIheader <- fread(headerFile, skip = 14, nrows=25, sep = "\t", header = FALSE)
  setnames(CASFRIattr, CASFRIheader$V1)
  set(CASFRIattr, , grep(CASFRIheader$V1, pattern = "^SPECIES|^GID|^AGE", invert = TRUE), NULL)
  setkey(CASFRIattr, "GID")
  NAVals <- c("XXXX MISS", "UNDEF", "XXXX ERRC")
  for(i in 1:5) {
    set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_",i)]]%in%NAVals), paste0("SPECIES_",i), NA_character_ )
    set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_",i)]]%in%NAVals), paste0("SPECIES_",i), NA_character_ )
  }
  for(i in 1:1) {
    CASFRIattr <- CASFRIattr[which(CASFRIattr[[paste0("SPECIES_PER_",i)]]>15),]
  }
  for(i in 2:5) {
    set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_",i)]]<=15), paste0("SPECIES_",i), NA_character_)
  }
  #CASFRIattr[SPECIES_5>15,.N]
  
  keepSpecies <- whSpecies(CASFRIattr, topN = 16) # 16 most abundant species
  CASFRIattrLong <- melt(CASFRIattr, id.vars = c("GID"),
                         measure.vars = paste0("SPECIES_",1:5))
  CA2 <- melt(CASFRIattr, id.vars = c("GID"),
              measure.vars = c(paste0("SPECIES_PER_",1:5)))
  CASFRIattrLong[,pct:=CA2$value]
  rm(CA2)
  CASFRIattrLong <- na.omit(CASFRIattrLong)
  
  CASFRIdt <- CASFRIRas[] %>% data.table(GID=., rastInd=1:ncell(CASFRIRas))
  CASFRIdt <- CASFRIdt[,isNA:=is.na(GID)]
  CASFRIdt <- CASFRIdt[isNA==FALSE]
  setkey(CASFRIdt, GID)
  set(CASFRIdt, ,"isNA",NULL)
  
  return(list(keepSpecies = keepSpecies, CASFRIattrLong = CASFRIattrLong,
              CASFRIdt = CASFRIdt))
}

whSpecies <- function(CASFRIattr, topN = 16) {
  spAbund <- CASFRIattr[,.N,by="SPECIES_1"] %>% setkeyv("N") #%>% print()
  spAbund2 <- CASFRIattr[,.N,by="SPECIES_2"] %>% setkeyv("N") #%>% print()
  setorder(spAbund, -N)
  setorder(spAbund2, N)
  keepSpecies <- data.table(keepSpecies=spAbund$SPECIES_1[1:topN])
  set(keepSpecies, ,"spGroup", keepSpecies$keepSpecies)
  setkey(keepSpecies, keepSpecies)
  keepSpecies <- keepSpecies[!"Pseu menz"]
  keepSpecies[c("Pice glau", "Pice enge", "Pice hybr", "Pice spp."),spGroup:="Pice_gla"]
  keepSpecies["Pice mari",spGroup:="Pice_mar"]
  keepSpecies["Betu papy",spGroup:="Betu_pap"]
  keepSpecies[c("Abie bals", "Abie lasi"),spGroup:="Abie_sp"]
  keepSpecies[c("Lari lari"),spGroup:="Lari_lar"]
  keepSpecies[c("Pinu cont", "Pinu conl"),spGroup:="Pinu_sp"]
  keepSpecies[c("Pinu bank", "Pinu spp."),spGroup:="Pinu_sp"]
  keepSpecies[c("Popu trem", "Popu balb"),spGroup:="Popu_tre"]
  keepSpecies
}


attrFile <- "C:/data/LandWeb/CASFRI for Landweb/Landweb_CASFRI_GIDs_attributes3.csv"
headerFile <- "C:/data/LandWeb/CASFRI for Landweb/Landweb_CASFRI_GIDs_README.txt"
loadedCASFRI <- Cache(loadCASFRI, CASFRIRas, attrFile, headerFile)

PaulTrimmed[PaulTrimmed[]%in% c(230, 220, 255)] <- NA_integer_ # water, non veg
Paulvals <- sort(unique(PaulTrimmed[]))
PaulStack <- list()
uniqueKeepSp <- unique(loadedCASFRI$keepSpecies$spGroup)
# "Abie_sp"  "Betu_pap" "Lari_lar" "Pice_gla" "Pice_mar" "Pinu_sp" "Popu_tre"
makePaulStack <- function(PaulTrimmed, uniqueKeepSp, ...) {
  rasterOptions(maxmemory = 1e9)
  N <- 1
  message("running ", N)
  PaulStack[[uniqueKeepSp[N]]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[uniqueKeepSp[N]]] <- Cache(writeRaster, PaulStack[[uniqueKeepSp[N]]] , 
                                        filename = paste0("Paul",uniqueKeepSp[N], ".tif"),
                                        overwrite=TRUE, ...)
  N <- 2
  message("running ", N)
  PaulStack[[uniqueKeepSp[N]]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[uniqueKeepSp[N]]] <- Cache(writeRaster, PaulStack[[uniqueKeepSp[N]]] , filename = paste0("Paul",uniqueKeepSp[N], ".tif"),
                                        overwrite=TRUE, ...)
  
  N <- 3
  message("running ", N)
  PaulStack[[uniqueKeepSp[N]]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[uniqueKeepSp[N]]] <- Cache(writeRaster, PaulStack[[uniqueKeepSp[N]]] , filename = paste0("Paul",uniqueKeepSp[N], ".tif"),
                                        overwrite=TRUE, ...)
  
  
  N <- "Pice_gla"
  message("running ", N)
  PaulStack[[N]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[N]][PaulTrimmed[] %in% c(41, 42, 43)] <- 60
  PaulStack[[N]][PaulTrimmed[] %in% c(44)] <- 90
  PaulStack[[N]][PaulTrimmed[] %in% c(14, 34)] <- 40
  PaulStack[[N]] <- Cache(writeRaster, PaulStack[[N]] , filename = paste0("Paul",N, ".tif"),
                          overwrite=TRUE, ...)
  
  # 5
  N <- "Pice_mar"
  message("running ", N)
  PaulStack[[N]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[N]][PaulTrimmed[] %in% c(23, 26)] <- 60
  PaulStack[[N]][PaulTrimmed[] %in% c(22)] <- 90
  PaulStack[[N]][PaulTrimmed[] %in% c(32, 42)] <- 40
  PaulStack[[N]] <- Cache(writeRaster, PaulStack[[N]] , filename = paste0("Paul",N, ".tif"),
                          overwrite=TRUE, ...)
  
  # 6
  N <- "Pinu_sp"
  message("running ", N)
  PaulStack[[N]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[N]][PaulTrimmed[] %in% c(31, 32, 34)] <- 60
  PaulStack[[N]][PaulTrimmed[] %in% c(33)] <- 90
  PaulStack[[N]][PaulTrimmed[] %in% c(23, 43)] <- 40
  PaulStack[[N]] <- Cache(writeRaster, PaulStack[[N]] , filename = paste0("Paul",N, ".tif"),
                          overwrite=TRUE, ...)
  
  
  # 7
  N <- "Popu_tre"
  message("running ", N)
  PaulStack[[N]] <- raster(PaulTrimmed) %>% setValues(NA_integer_)
  PaulStack[[N]][PaulTrimmed[] %in% c(14)] <- 60
  PaulStack[[N]][PaulTrimmed[] %in% c(11)] <- 90
  PaulStack[[N]][PaulTrimmed[] %in% c(31, 41)] <- 40
  PaulStack[[N]] <- Cache(writeRaster, PaulStack[[N]] , filename = paste0("Paul",N, ".tif"),
                          overwrite=TRUE, ...)
  
  stack(PaulStack)
  # PaulStack <- Cache(writeRaster, PaulStack, filename = "PaulSpStack.tif", datatype = "INT2U",
  #                    overwrite=TRUE)
  
}

PaulSpStack <- Cache(makePaulStack, PaulTrimmed, uniqueKeepSp)

# Step 3 -- Load CASFRI data and headers, and make a raster for each relevant species group
CASFRItoSpRasts <- function(spRas, loadedCASFRI) {
  spRasts <- list()
  for(sp in unique(loadedCASFRI$keepSpecies$spGroup)) {
    spRasts[[sp]] <- spRas
    message("starting ", sp)
    aa2 <- loadedCASFRI$CASFRIattrLong[value %in% 
                                         loadedCASFRI$keepSpecies[spGroup==sp,
                                                                  keepSpecies]][
                                                                    ,min(100L, sum(pct)), by = GID]
    
    setkey(aa2, GID)
    cc <- aa2[loadedCASFRI$CASFRIdt] %>% na.omit()
    rm(aa2)
    spRasts[[sp]][cc$rastInd] <- cc$V1
    message("  ", sp, " writing to disk")
    system.time(spRasts[[sp]] <- #Cache(
                  writeRaster(spRasts[[sp]], filename = paste0(sp,".tif"), #objectLength = 1e6,
                              datatype = "INT2U", overwrite = TRUE)
    )
    message("  ", sp, " done")
  }
  
  stack(spRasts)
}

spRas <- raster(CASFRIRas) %>% setValues(.,NA_integer_)
#dd <- CASFRItoSpRasts(spRas, loadedCASFRI)
CASFRISpStack <- Cache(CASFRItoSpRasts, spRas, loadedCASFRI)


overlayStacks <- function(stack1, stack2) {
  for(sp in layerNames(stack1)) {
    message(sp)
    nas <- is.na(stack1[[sp]][])
    stack1[[sp]][nas] <- stack2[[sp]][][nas]
    message("Writing to disk")
    stack1[[sp]] <- Cache(writeRaster, stack1[[sp]], datatype = "INT2U",
                          filename=paste0(sp,"_overlay.tif"), overwrite = TRUE)
  }
  stack1
}

outStack <- Cache(overlayStacks, CASFRISpStack, PaulSpStack)#, notOlderThan = Sys.time())
rm(loadedCASFRI, PaulTrimmed, spRas)
endTime <- Sys.time()
print(endTime - startTime)