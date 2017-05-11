eben
defineModule(sim, list(
  name="disturbanceDriver",
  description="generate parameters for the generic percolation model",# spades::spread()",
  keywords=c("fire"),
  authors=c(person(c("Steve", "G"), "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut", "cre"))),
  childModules=character(),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c(NA, NA)),
  timeunit="year", 
  citation=list(),
  reqdPkgs=list("stats"),
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur")),
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
  inputObjects=data.frame(objectName=c("spreadCalibrationData","fireRegimeParameters","fireMapAttr"),
                          objectClass=c("data.frame","list","list"), other=rep(NA_character_,3), stringsAsFactors=FALSE),
  #inputObjects=data.frame(objectName=c("fireRegimeParameters","spreadCalibrationData"), objectClass=c("list","data.frame"), other=rep(NA_character_,2), stringsAsFactors=FALSE),
  outputObjects=data.frame(objectName=c("disturbanceGeneratorParameters"),
                          objectClass=c("list"),
                          other=NA_character_, stringsAsFactors=FALSE)
))



## event types
#   - type `init` is required for initiliazation

doEvent.disturbanceDriver = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    sim<-disturbanceDriverInit(sim)

    # schedule future event(s)
    scheduleEvent(sim, params(sim)$disturbanceDriver$.plotInitialTime, "disturbanceDriver", "plot")
    scheduleEvent(sim, params(sim)$disturbanceDriver$.saveInitialTime, "disturbanceDriver", "save")
  } else if (eventType=="templateEvent") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event

    # e.g., call your custom functions/methods here
    # you can define your own methods below this `doEvent` function

    # schedule future event(s)

    # e.g.,
    # scheduleEvent(sim, time(sim) + increment, "disturbanceDriver", "templateEvent")

    # ! ----- STOP EDITING ----- ! #
    } else {
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


# 1 - (1-p0)**N = pEscape
# 1 - pEscape = (1-p0)**N
# (1 - pEscape)**1/N = 1 - p0
# p0 = 1 - (1 - pEscape)**1/N

hatP0<-function(pEscape,n=8){
  1 - (1-pEscape)**(1/n)
}

#a real clever boots would minimise the abs log odds ratio. 
#be my guest.

escapeProbDelta<-function(p0,w,hatPE){
  abs(sum(w*(1-(1-p0)**(0:8)))-hatPE)  
}


disturbanceDriverInit = function(sim) {


 
  
  #this table contains calibration data for several landscape sizes
  #and several min fire sizes (1 or 2 cells), organised by collumn.
  #The data were made by Steve Cumming in June 2013 for a whole other purpose.
  #I chose the one that seems most appropriate to me
  #browser()
  y<-log(sim$spreadCalibrationData[,paste("ls",1e3,"fs",2,sep="")])
  x<-sim$spreadCalibrationData$pjmp
  m.glm<-glm(x~I(log(y)))
  mfs<-sim$fireRegimeParameters$xBar/sim$fireMapAttr$cellSize #mean size escaped fires in cells
  pJmp<-sum(m.glm$coeff*c(1,log(mfs)))
  
  w<-sim$fireMapAttr$nNbrs
  w<-w/sum(w)
  hatPE<-sim$fireRegimeParameters$pEscape
  foo<-optimise(sim$escapeProbDelta,
                interval=c(sim$hatP0(hatPE,globals(sim)$neighbours),sim$hatP0(hatPE,floor(sum(w*0:8)))),
                tol=1e6,
                w=w,
                hatPE=hatPE)
  #do some sanity tests to ensure convergence
  #also, it is almost obvious that the true minimum must occurr within the interval specified in the 
  #call to optimise, but I have not proved it, nor am I certain that the function being minimised is 
  #monotone.
  
  rate<-sim$fireRegimeParameters$rate * sim$fireMapAttr$cellSize #fireRegimeModel and this module must agree on 
                                                                 #an annual time step. How to test / enforce?
  pIgnition <- rate #approximate Poisson arrivals as a Bernoulli process at cell level.
                    #for Poisson rate << 1, the expected values are the same, partially accounting
                    #for multiple arrivals within years. Formerly, I used a poorer approximation
                    #where 1-p = P[x==0 | lambda=rate] (Armstrong and Cumming 2003).
  
  
  sim$disturbanceGeneratorParameters<-list(pJmp=pJmp,p0=foo$minimum,
                                           naiveP0=sim$hatP0(sim$fireRegimeParameters$pEscape,8), 
                                           pIgnition=pIgnition)

  return(invisible(sim))
}
 
### template for save events
disturbanceDriverSave = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
disturbanceDriverPlot = function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot("object")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
disturbanceDriverEvent1 = function(sim) {
  # ! ----- EDIT BELOW ----- ! #



  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}


### add additional events as needed by copy/pasting from above

