stopifnot(packageVersion("SpaDES") >= "0.99.0")

defineModule(sim, list(
  name="scfmLandcoverInit",
  description="Takes the LCC05 classification of 39 land cover classes, and reclassifies it to flammable and inflammable [1,0]",
  keywords=c("fire", "LCC05", "land cover classification 2005", "BEACONs"),
  childModules=character(),
  authors=c(person(c("Eliot", "J", "B"), "McIntire", email="Eliot.McIntire@NRCan.gc.ca", role=c("aut", "cre")),
            person("Steve", "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut"))),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c("2005-01-01", NA)),
  documentation = list("README.txt", "scfmLandCoverInit.Rmd"),
  timeunit="year",
  citation=list(),
  reqdPkgs=list("raster", "purrr"),
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA, desc="Initial time for plotting"),
    defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, desc="Interval between plotting"),
    defineParameter(".saveInitialTime", "numeric", NA_real_,  NA, NA, desc="Initial time for saving"),
    defineParameter(".saveIntervalXXX", "numeric", NA_real_, NA, NA, desc="Interval between save events"),
    defineParameter("useCache", "logical", TRUE, NA, NA, desc="Use cache")),
  inputObjects=data.frame(objectName=c("vegMap","studyArea"),
                          objectClass=c("RasterLayer","SpatialPolygonsDataFrame"),
                          sourceURL="",
                          other=NA_character_, stringsAsFactors=FALSE),
  outputObjects=data.frame(objectName=c("flammableMap", "landscapeAttr", "cellsByZone"), #mapAttr are all things the fir
                           objectClass=c("RasterLayer", "list", "data.frame"),
                           other=rep(NA_character_, 3L), 
                           stringsAsFactors=FALSE)
))

doEvent.scfmLandcoverInit = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    
    sim <- scfmLandcoverInitCacheFunctions(sim)
    sim <- scfmLandcoverInitInit(sim)
    sim <- scheduleEvent(sim, params(sim)$scfmLandcoverInit$.plotInitialTime,
                         "scfmLandcoverInit", "plot")
    sim <- scheduleEvent(sim, params(sim)$scfmLandcoverInit$.saveInitialTime,
                         "scfmLandcoverInit", "save")
  } else if (eventType=="plot") {
    #browser()
    Plot(sim$vegMap, new=TRUE)
    Plot(sim$flammableMap, legend=FALSE) # this is failing probably due to a bug in Plot
    # EJM is working on it 20160224
    # schedule future event(s)
    sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmLandcoverInit$.plotInterval, "scfmLandcoverInit", "plot")
  } else if (eventType=="save") {
    # schedule future event(s)
    sim <- scheduleEvent(sim, time(sim) + params(sim)$scfmLandcoverInit$.saveInterval, "scfmLandcoverInit", "save")
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with=FALSE],
                  "' in module '", events(sim)[1, "moduleName", with=FALSE], "'", sep=""))
  }
  return(invisible(sim))
}

genFireMapAttr<-function(sim){
  #calculate the cell size, total area, and number of flammable cells, etc.
  #
  #All areas in ha
  #
  cellSize<-prod(res(sim$flammableMap))/1e4 #in ha
  
  ecoregionErrors <- as.numeric(sim$studyArea$ECOREGION)==0
  if(any(ecoregionErrors)) {
    warning("ECOREGION in studyArea has invalid values. Removing invalid values")
    sim$studyArea <- sim$studyArea[!ecoregionErrors,]
  }
  
  if (globals(sim)$neighbours==8)
    w<-matrix(c(1,1,1,1,0,1,1,1,1),nrow=3,ncol=3)
  else if (globals(sim)$neighbours==4)
    w<-matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,ncol=3)
  else 
    stop("illegal global neighbours spec")
  #it would be nice to somehow get caching to work on the function argument of focal
  #but I have not been able to make it work.
  neighMap <- focal(1-sim$flammableMap, w, na.rm=TRUE) #default function is sum(...,na.rm)
  neighMapVals <- values(neighMap)
  
  # extract table for each polygon
  valsByPoly <- Cache(extract, neighMap, sim$studyArea, cellnumbers = TRUE)
  
  names(valsByPoly) <- sim$studyArea$ECOREGION
  uniqueZoneNames <- unique(sim$studyArea$ECOREGION)
  valsByZone <- lapply(uniqueZoneNames, function(ecoName) {
    aa <- valsByPoly[names(valsByPoly)==ecoName] 
    if(is.list(aa)) aa <- do.call(rbind, aa)
    aa
  })
  names(valsByZone) <- uniqueZoneNames
  
  # Derive frequency tables of number of flammable cells, per polygon type, currently ECOREGION
  nNbrs <- lapply(valsByZone, function(x) {
    nNbrs<-tabulate(x[,2]+1, 9)#depends on sfcmLandCoverInit
    names(nNbrs)<-0:8
    nNbrs
  })
  
  nFlammable <- lapply(valsByZone, function(x) {
    sum(1-getValues(sim$flammableMap)[x[,1]], na.rm=TRUE)
  })
  
  landscapeAttr <- purrr::transpose(list(cellSize=rep(list(cellSize), length(nFlammable)), 
                                         nFlammable=nFlammable,
                                         nNbrs=nNbrs,
                                         cellsByZone=lapply(valsByZone, function(x) x[,1])))
  
  sim$landscapeAttr <- lapply(landscapeAttr, function(x) {
    append(x, list(burnyArea=x$cellSize*x$nFlammable))
  })
  names(sim$landscapeAttr) <- names(valsByZone)
  sim$cellsByZone <- data.frame(cell=1:ncell(sim$flammableMap), zone=NA)
  for(x in names(sim$landscapeAttr)) {
    sim$cellsByZone[sim$landscapeAttr[[x]]$cellsByZone,"zone"] <- x
  }
  
  
  return(invisible(sim))
}

### template initilization
scfmLandcoverInitInit = function(sim) {
  # these classes are LCC05 specific
  #browser()
  nonFlammClasses<-c(36,37,38,39)
  oldClass <- 0:39
  newClass <- ifelse(oldClass %in% nonFlammClasses,1,0)   #1 codes for non flammable 
  #see mask argument for SpaDES::spread()
  flammableTable <- cbind(oldClass, newClass)
  sim$flammableMap <- ratify(reclassify(sim$vegMap, flammableTable,count=TRUE))
  if ("Mask" %in% names(objs(sim))){
    sim$flammableMap <- sim$flammableMap * sim$Mask
  }
  #the count options should cause that "a column with frequencies is added. 
  
  #setColors(sim$flammableMap, n=2) <- c("blue","red")
  setColors(sim$flammableMap,2) <- colorRampPalette(c("blue", "red"))(2) 
  
  sim<-genFireMapAttr(sim)
  
  return(invisible(sim))
}

testFun<-function(x) {
  sum(na.omit(x)==1)
}

scfmLandcoverInitCacheFunctions <- function(sim) {
  # for slowp functions, add cached versions
  # browser()
  if(params(sim)$scfmLandcoverInit$useCache) {
    sim$cacheLoc <- file.path(outputPath(sim), "scfmLandcoverInitCache") 
    if(!dir.exists(sim$cacheLoc) )
      createEmptyRepo(file.path(outputPath(sim), "scfmLandcoverInitCache"))
    
    sim$focal <- function(...) {
      archivist::cache(cacheRepo=sim$cacheLoc, FUN=raster::focal, ...)
    }
    
    sim$Frabjous <-function(...){
      archivist::cache(cacheRepo=sim$cacheLoc, FUN=sim$testFun, ...)
    }
    
  } else {
    sim$Frabjous <- function(x){sum(na.omit(x)==1)}
    sim$focal <- raster::focal
  }
  
  return(invisible(sim))
}
