
defineModule(sim, list(
  name="scfmDriver",
  description="generate parameters for the generic percolation model",# spades::spread()",
  keywords=c("fire"),
  authors=c(person(c("Steve", "G"), "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut", "cre"))),
  childModules=character(),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c(NA, NA)),
  timeunit="year", 
  citation=list(),
  documentation = list("README.txt", "scfmDriver.Rmd"),
  reqdPkgs=list("stats"),
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur")),
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
  inputObjects=data.frame(objectName=c("scfmRegime","landscapeAttr"),
                          objectClass=c("list","list"), 
                          sourceURL="",
                          other=rep(NA_character_,2), 
                          stringsAsFactors=FALSE),
  
  outputObjects=data.frame(objectName=c("scfmPars"),
                          objectClass=c("list"),
                          other=NA_character_, stringsAsFactors=FALSE)
))



## event types
#   - type `init` is required for initiliazation

doEvent.scfmDriver = function(sim, eventTime, eventType, debug=FALSE) {
  if (eventType=="init") {
    sim<-scfmDriverInit(sim)
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

scfmDriverInit = function(sim) {

  #this table contains calibration data for several landscape sizes
  #and several min fire sizes (1 or 2 cells), organised by collumn.
  #The data were made by Steve Cumming in June 2013 for a whole other purpose.
  #I chose the one that seems most appropriate to me
  cellSize <- sim$landscapeAttr[[1]]$cellSize
  
  sim$scfmPars<- lapply(names(sim$landscapeAttr), function(polygonType) {
  
    regime <- sim$scfmRegime[[polygonType]]
    landAttr <- sim$landscapeAttr[[polygonType]]
      
    if (FALSE) {
      y<-log(sim$spreadCalibrationTable[,paste("ls",1e3,"fs",2,sep="")])
      x<-sim$spreadCalibrationTable$pjmp
      m.glm<-glm(x~y,family=gaussian)
      mfs<-regime$xBar/cellSize #mean size escaped fires in cells
      pJmp<-sum(m.glm$coeff*c(1,log(mfs)))
      #mfs <- sim$scfmRegime$lxBar - log(cellSize)
      #pjmp <- sum(m.glm$coeff*c(1,mfs))
    }
    else {
      #we know this table was produced with MinFireSize=2cells.
      #browser()
      y <- sim$cTable2$y
      x <- sim$cTable2$p
      m.lw <- lowess(y~x,iter=2)
      if (sum(diff(m.lw$y)<0)>0)
        warning("lowess curve non-monotone. Proceed with caution")
      targetSize <- regime$xBar/cellSize - 1 
      pJmp <- approx(m.lw$y,m.lw$x,targetSize,rule=2)$y
    }
    #browser()
    w<-landAttr$nNbrs
    w<-w/sum(w)
    hatPE<-regime$pEscape
    if(hatPE == 0) { # no fires in polygon zone escapted
      foo <- 0 
    } else if (hatPE == 1) { # all fires in polygon zone escaped
      foo <- 1
    } else {
      foo<-optimise(sim$escapeProbDelta,
                    interval=c(sim$hatP0(hatPE,globals(sim)$neighbours),
                               sim$hatP0(hatPE,floor(sum(w*0:8)))),
                    tol=1e6,
                    w=w,
                    hatPE=hatPE)$minimum
    } 
    #do some sanity tests to ensure convergence
    #also, it is almost obvious that the true minimum must occurr within the interval specified in the 
    #call to optimise, but I have not proved it, nor am I certain that the function being minimised is 
    #monotone.
    
    #don't forget to scale by number of years, as well.
    rate<-regime$ignitionRate * cellSize #fireRegimeModel and this module must agree on 
                                                                   #an annual time step. How to test / enforce?
    pIgnition <- rate #approximate Poisson arrivals as a Bernoulli process at cell level.
                      #for Poisson rate << 1, the expected values are the same, partially accounting
                      #for multiple arrivals within years. Formerly, I used a poorer approximation
                      #where 1-p = P[x==0 | lambda=rate] (Armstrong and Cumming 2003).
    
    list(pSpread=pJmp,
         p0=foo,
         naiveP0=sim$hatP0(regime$pEscape,8), 
         pIgnition=pIgnition,
         maxBurnCells=as.integer(round(regime$emfs/cellSize)))
  })
  names(sim$scfmPars) <- names(sim$landscapeAttr)
  
  return(invisible(sim))
}
 

