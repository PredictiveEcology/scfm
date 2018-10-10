stopifnot(packageVersion("SpaDES") >= "0.99.0")
defineModule(sim, list(
  name="scfmCrop",
  description="A translator module. Crops, reprojects and mask all necessary layers to\ 
      specifed extent. Selects or subsets all spatialised data e.g. fires.\ 
      Use cropped RasterLayer, defined by a spatial polygon that defines the area of interest",
  keywords=c("translator", "lcc05", "Land Cover Classification", "vegetation","Canadian National Fire Database"),
  childModules=character(),
  authors=c(person("Steve", "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut")),
            person(c("Eliot", "J","B"), "McIntire", email="emcintir@nrcan.gc.ca", role=c("aut")),
            person("Pierre", "Vernier", email="pierre.vernier@gmail.com")
            ),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c(NA, NA)),
  timeunit=NA_character_,
  citation=list(),
  reqdPkgs=list("raster","rgeos","sp","archivist"),
  documentation = list("README.txt", "scfmCrop.Rmd"),
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", NA_real_, NA, NA, desc="Initial time for plotting"),
    defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, desc="Interval between plotting"),
    defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, desc="Initial time for saving"),
    defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, desc="Interval between save events"),
    defineParameter("useCache", "logical", TRUE, NA, NA, desc="Use cache or not.")),
    inputObjects=data.frame(objectName =c("studyArea",      "vegInput",    "ageMapInit",  "firePointsInput",  "cacheLoc"),
                            objectClass=c("SpatialPolygons","RasterLayer", "RasterLayer", "SpatialPoints",    "character"),
                            sourceURL="",
                           other=rep(NA_character_, 5L), stringsAsFactors=FALSE),
    outputObjects=data.frame(objectName=c("vegMap", "ageMapInit", "firePoints", "ageMap"), #these are the names that matter for dependencies
                             objectClass=c("RasterLayer", "RasterLayer", "SpatialPoints", "RasterLayer"),
                             other=rep(NA_character_, 4L), stringsAsFactors=FALSE)
))

#Data are in ~stevec/Dropbox/SpaDES/Data/CanadianNationalFireDatabase as of 2015.07.10
#Fire data sources.Rmd  NFDB_point.zip		NFDB_poly.zip

doEvent.scfmCrop = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    # do stuff for this event
    #sim <- scfmCropCacheFunctions(sim)
    sim <- scfmCropInit(sim)
    # schedule future event(s)
    sim <- scheduleEvent(sim, params(sim)$scfmCrop$.plotInitialTime, "scfmCrop", "plot")
    sim <- scheduleEvent(sim, params(sim)$scfmCrop$.saveInitialTime, "scfmCrop", "save") 
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with=FALSE],
                  "' in module '", events(sim)[1, "moduleName", with=FALSE], "'", sep=""))
  }
  return(invisible(sim))
}

### template initilization
scfmCropInit <- function(sim) {

  vegProjection <- crs(sim$vegMapInit)
  
  if (is.na(crs(sim$studyArea)))        #in case it was sampled from the vegmap.
    crs(sim$studyArea) <- vegProjection
  
  simProjection<-crs(sim$studyArea)  #this would probably be set to be the same as the veg map at an earlier stage.
  
  #if(ncell(sim$vegMap)>5e5) beginCluster(min(parallel::detectCores(),6))
  
  #Project the study area into each input raster, then crop and mask; 
  #Then project result intpo sim projection.
  cropMaskReproject <- function(map, studyArea, fromProjection, toProjection) {
    studyAreaTmp <- spTransform(studyArea, CRSobj=fromProjection)
    updatedMap <-  Cache(crop, map, studyAreaTmp)
    crs(updatedMap) <- fromProjection
    updatedMap <- Cache(mask, updatedMap, studyAreaTmp) #
    updatedMap <- Cache(projectRaster, updatedMap,
                    crs=toProjection, method="ngb")
    
  }
  
  sim$vegMap <- Cache(cropMaskReproject, sim$vegMapInit, sim$studyArea,
                                  vegProjection, simProjection)
  
  buildBinaryMask <- function(map) {
    Mask <- map
    Mask[!is.na(map[])] <- 1
    Mask
  }
  sim$Mask <- Cache(buildBinaryMask, sim$vegMap)
  
  tmp <- getColors(sim$vegMapInit)[[1]] # mask removes colors!
  setColors(sim$vegMap, n=length(tmp)-1) <- tmp[-1]
  #setColors(sim$vegMap, n=maxValue(sim$vegMap)) <- tmp[1:maxValue(sim$vegMap)]
  
  ageProjection <- crs(sim$ageMapInit)
  sim$ageMap <- Cache(cropMaskReproject, sim$ageMapInit, sim$studyArea,
                      ageProjection, simProjection)
  
  # studyAreaTmp <- Cache(spTransform, sim$studyArea, CRSobj =ageProjection)
  # sim$ageMap <-  Cache(crop, sim$ageMapInit, studyAreaTmp)
  # crs(sim$ageMap) <- ageProjection
  # sim$ageMap <- Cache(mask, sim$ageMap,studyAreaTmp)
  # sim$ageMap <- Cache(projectRaster, sim$ageMap,to=sim$vegMap,method="ngb")
  
  cropMaskReprojectSP <- function(sp, studyArea, toProjection) {
    fireProjection <- CRS(proj4string(sp))
    studyAreaTmp <- Cache(spTransform, studyArea, CRSobj =fireProjection)
    firePoints <- sp[studyAreaTmp,]  #note possibly correct syntax A[B,] rather than A[B]
    #https://cran.r-project.org/web/packages/sp/vignettes/over.pdf page 3
    #crs(sim$firePoints) <- fireProjection
    firePoints <- Cache(spTransform, firePoints,CRSobj=toProjection)
  }
  
  sim$firePoints <- Cache(cropMaskReprojectSP, sim$firePointsInput, 
                          sim$studyArea, simProjection)
  
  #endCluster()
    
  # age will not run with projectRaster directly. Instead, project the vegMap to age, then crop, then project back to vegMap
  #vegMap.crsAge <- projectRaster(sim$vegMap, crs=crs(sim$age))
  #age.crsAge <- crop(sim$age, spTransform(sim$studyArea, CRSobj = crs(sim$age)))
  #age.crsAge <- mask(age.crsAge, spTransform(sim$studyArea, CRSobj = crs(sim$age)))
  #sim$ageMapInit <- projectRaster(age.crsAge, to=sim$vegMap, method="ngb")
  
  
  
  #if(sum(!is.na(getValues(sim$ageMap)))==0) stop("There are no age data provided with input age map")
  #if(sum(!is.na(getValues(sim$vegMap)))==0) stop("There are no vegatation data provided with input vegatation map")
  
  # Faster -- and should be as robust since raster functions update their maxValue to NA if all are NA
  if(is.na(maxValue(sim$ageMap))) stop("There are no age data provided with input age map")
  if(is.na(maxValue(sim$vegMap))) stop("There are no vegatation data provided with input vegatation map")
  
  return(invisible(sim))
}


# scfmCropCacheFunctions <- function(sim) {
#   # for slow functions, add cached versions
#   if(params(sim)$scfmCrop$useCache) {
#     sim$cacheLoc <- file.path(outputPath(sim), "scfmCropCache") 
#     
#     if(!dir.exists(sim$cacheLoc) )
#       createEmptyRepo(file.path(outputPath(sim), "scfmCropCache"))
#     
#     sim$mask <- function(...) {
#       archivist::cache(cacheRepo=sim$cacheLoc, FUN=raster::mask, ...)
#     }
#     sim$crop <- function(...) {
#       archivist::cache(cacheRepo=sim$cacheLoc, FUN=raster::crop, ...)
#     }
#     sim$projectRaster <- function(...) {
#       archivist::cache(cacheRepo=sim$cacheLoc, FUN=raster::projectRaster, ...)
#     }
#     sim$spTransform <- function(...) {
#       archivist::cache(cacheRepo=sim$cacheLoc, FUN=sp::spTransform,  ...)
#     }
#   } else {
#     sim$mask <- raster::mask
#     sim$crop <- raster::crop
#     sim$projectRaster <- raster::projectRaster
#     sim$spTransform <- sp::spTransform
#     
#   }
#   
#   return(invisible(sim))
# }
