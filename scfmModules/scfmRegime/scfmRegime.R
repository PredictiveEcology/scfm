
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
  reqdPkgs=list(),
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter("fireCause", "character", c("L"), NA, NA,  "subset of c(H,H-PB,L,Re,U)"),
    defineParameter("fireEpoch", "numeric", c(1961,1990), NA, NA, "start of normal period")
  ),
  inputObjects=data.frame(
    objectName=c("firePoints","flammableMap","landscapeAttr"),
    objectClass=c("SpatialPointsDataFrame", "RasterLayer", "list"),
    sourceURL="",
    other=rep(NA_character_,3),
    stringsAsFactors=FALSE),
  outputObjects=data.frame(
    objectName=c("scfmRegime"),
    objectClass=c("list"),
    other=rep(NA_character_,1),
    stringsAsFactors=FALSE)
))


## event types
#   - type `init` is required for initiliazation

doEvent.scfmRegime = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    # browser()
    #genlandscapeAttr() now called in scfmLandcoverInit.R
    scfmRegimeInit(sim)
  }
   else {
      warning(paste("Undefined event type: '", events(sim)[1, "eventType", with=FALSE],
                    "' in module '", events(sim)[1, "moduleName", with=FALSE], "'", sep=""))
    }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initilization


scfmRegimeInit = function(sim) {

  #browser()

  #curModule <- events(mySim)$module[1] #NOPE When I need it, it's checkpoint
  #in SpaDES 1.4+, place this in a R subdirectory of the module. All such R files are sources into the simlist
  source(file.path(paths$modulePath, "scfmRegime", "TEutilsNew.R"), local=TRUE,echo=FALSE)


  #subset fires by cause and epoch.
  tmp<-sim$firePoints

  #extract and validate fireCause spec
  browser()
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
  
  
  #Poisson!!!

  lapply(seq_len(NROW(sim$studyArea)), function(polyRowIndex) {
    browser()
      
    
    nFires<-dim(tmp)[1]
    rate<-nFires/(epochLength * sim$landscapeAttr$burnyArea)   # fires per ha per yr
  
    #calculate escaped fires
    #careful to subtract cellSize where appropriate
    xVec<-tmp$SIZE_HA[tmp$SIZE_HA > sim$landscapeAttr$cellSize]
    pEscape<-length(xVec)/nFires
    xBar<-mean(xVec)
    lxBar<-mean(log(xVec))
    xMax<-max(xVec)
  
    zVec<-log(xVec/sim$landscapeAttr$cellSize)
    if (length(zVec)<100)
      warning("Less than 100 \"large\" fires. That estimates may be unstable.\nConsider using a larger area and/or longer epoch.")
    #later, this would sim$HannonDayiha
    hdList<-HannonDayiha(zVec) #defined in sourced TEutilsNew.R
    maxFireSize<-exp(hdList$That) * sim$landscapeAttr$cellSize
    #verify estimation results are reasonable. That=-1 indicates convergence failure.
    #
    #need to addd a name or code for basic verification by Driver module, and time field
    #to allow for dynamic regeneration of disturbanceDriver pars.
    sim$scfmRegime<-list(ignitionRate=rate,
                         pEscape=pEscape,
                         xBar=xBar,  #mean fire size
                         lxBar=lxBar, #mean log(fire size)
                         xMax=xMax,  #maximum observed size
                         meanBigFireSize=mean(xVec[xVec>200]),
                         emfs=maxFireSize) # Estimated Maximum Fire Size in ha
  })
  
  return(invisible(sim))
}

### template for save events
scfmRegimeSave = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
scfmRegimePlot = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot("object")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
scfmRegimeEvent1 = function(sim) {
  # ! ----- EDIT BELOW ----- ! #



  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}


