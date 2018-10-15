
defineModule(sim, list(
  name="scfmRegime",
  description="estimates fire regime parameters for BEACONs a la Steve's method",
  keywords=c("fire regime", "BEACONs"),
  authors=c(person(c("Steven", "G."), "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut", "cre"))),
  childModules=character(),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c(NA, NA)),
  timeunit="year",
  citation=list(),
  documentation = list("README.txt", "scfmRegime.Rmd"),
  reqdPkgs=list("rgdal"),
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter("fireCause", "character", c("L"), NA, NA,  "subset of c(H,H-PB,L,Re,U)"),
    defineParameter("fireEpoch", "numeric", c(1961,1990), NA, NA, "start of normal period")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "firePoints", objectClass = "SpatialPointsDataFrame", desc = "",
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"),
    expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = ""),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = "")
  ),
  outputObjects = bind_rows(
   createsOutput(objectName = "scfmRegimePars", objectClass = "list", desc =  "")
    )
))


## event types
#   - type `init` is required for initiliazation

doEvent.scfmRegime = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    Init(sim)
  }
   else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with=FALSE],
                    "' in module '", events(sim)[1, "moduleName", with=FALSE], "'", sep=""))
    }
  return(invisible(sim))
}


Init <- function(sim) {
  browser()
  tmp<-sim$firePoints

  #extract and validate fireCause spec
  fc<-params(sim)$scfmRegime$fireCause
  causeSet <- if(is.factor(tmp$CAUSE)) levels(tmp$CAUSE) else unique(tmp$CAUSE)
    
  if(any(!(fc %in% causeSet)))
      stop("illegal fireCause: ", fc)
  tmp<-subset(tmp,CAUSE %in% fc)

  #extract and validate fireEpoch
  epoch<-params(sim)$scfmRegime$fireEpoch
  if (length(epoch)!=2 || !is.numeric(epoch) || any(!is.finite(epoch)) || epoch[1]>epoch[2])
      stop("illegal fireEpoch: ",epoch)
  tmp<-subset(tmp, YEAR>=epoch[1] & YEAR<=epoch[2])
  
  epochLength<-as.numeric(epoch[2]-epoch[1]+1)
  
  
  # Assign polygon label to SpatialPoints of fires object
  tmp[["ECOREGION"]] <- over(tmp, sim$studyArea[, "ECOREGION"])
  
  # Hack to make a study area level cellSize ... TODO -- this should be removed from landscapeAttr
  cellSize <- sim$landscapeAttr[[1]]$cellSize
  
  sim$scfmRegimePars <-lapply(names(sim$landscapeAttr), function(polygonType) {
    tmpA <- tmp[unlist(tmp[["ECOREGION"]])==polygonType,]
    landAttr <- sim$landscapeAttr[[polygonType]]
    
    nFires<-dim(tmpA)[1]
    rate<-nFires/(epochLength * landAttr$burnyArea)   # fires per ha per yr
  
    pEscape <- 0
    maxFireSize <- NA
    xVec <- numeric(0)
    
    if(nFires > 0) {
    #calculate escaped fires
    #careful to subtract cellSize where appropriate
      xVec<-tmpA$SIZE_HA[tmpA$SIZE_HA > cellSize]
      pEscape<-length(xVec)/nFires
      
      zVec<-log(xVec/cellSize)
      if (length(zVec) < 100)
        warning("Less than 100 \"large\" fires. That estimates may be unstable.\nConsider using a larger area and/or longer epoch.")
      #later, this would sim$HannonDayiha
      if(length(zVec) > 0) {
        hdList<-HannonDayiha(zVec) #defined in sourced TEutilsNew.R
        maxFireSize<-exp(hdList$That) * cellSize
      }
      
    } 
    xBar<-mean(xVec)
    lxBar<-mean(log(xVec))
    xMax<-max(xVec)
    #verify estimation results are reasonable. That=-1 indicates convergence failure.
    #
    #need to addd a name or code for basic verification by Driver module, and time field
    #to allow for dynamic regeneration of disturbanceDriver pars.
    list(ignitionRate=rate,
         pEscape=pEscape,
         xBar=xBar,  #mean fire size
         lxBar=lxBar, #mean log(fire size)
         xMax=xMax,  #maximum observed size
         #meanBigFireSize=mean(xVec[xVec>200]),
         emfs=maxFireSize) # Estimated Maximum Fire Size in ha
  })
  
  names(sim$scfmRegimePars) <- names(sim$landscapeAttr)
  
  
  sim$firePoints <- tmp
  
  return(invisible(sim))
}



.inputObjects <- function(sim) {
  
  dPath <- inputPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")
  
  if (!suppliedElsewhere("firePoints", sim)) {
    if (!dir.exists(file.path(dPath, "NFDB_point"))) {
     
      download.file(url = extractURL(objectName = "firePoints"), 
                    destfile = file.path(dPath, "NFDB_point.zip"))
      unzip(zipfile = file.path(dPath, "NFDB_point.zip"), exdir = file.path(dPath, "NFDB_point"))
    }
  
    zipContents <- list.files(file.path(dPath, "NFDB_point"), all.files = TRUE, full.names = TRUE)
    outFile <- grep(pattern = "*.shp$", x = zipContents, value = TRUE)
    
    #browser()
    firePoints <- Cache(rgdal::readOGR, outFile)
    firePoints <- spTransform(firePoints, CRSobj = crs(sim$studyArea))
    firePoints <- firePoints[sim$studyArea,]
    
    sim$firePoints <- firePoints
  }
  
  return(invisible(sim))
}
  
