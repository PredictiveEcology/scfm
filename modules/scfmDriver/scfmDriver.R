
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
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParemater("neighbours", "numeric", 8, 4, 8, "number of cell immediate neighbours")),
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmRegimePars", objectClass = "list", desc = ""),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = ""),
    expectsInput(objectName = "cTable2", objectClass = "data.frame", desc = "A csv containing results of fire experiment")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName="scfmPars", objectClass = "list", desc = "")
  )                      
))

## event types
#   - type `init` is required for initiliazation

doEvent.scfmDriver = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
      init = {
        sim <- Init(sim)
      },
    
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with=FALSE],
                    "' in module '", events(sim)[1, "moduleName", with=FALSE], "'", sep=""))
  )
  return(invisible(sim))
}


# 1 - (1-p0)**N = pEscape
# 1 - pEscape = (1-p0)**N
# (1 - pEscape)**1/N = 1 - p0
# p0 = 1 - (1 - pEscape)**1/N

hatP0 <- function(pEscape,n=8){
  1 - (1-pEscape)**(1/n)
}

#a real clever boots would minimise the abs log odds ratio. 
#be my guest.

escapeProbDelta <- function(p0,w,hatPE){
  abs(sum(w*(1-(1-p0)**(0:8)))-hatPE)  
}

Init <- function(sim) {
  
  #this table contains calibration data for several landscape sizes
  #and several min fire sizes (1 or 2 cells), organised by collumn.
  #The data were made by Steve Cumming in June 2013 for a whole other purpose.
  #I chose the one that seems most appropriate to me
  cellSize <- sim$landscapeAttr[[1]]$cellSize
  
  sim$scfmPars<- lapply(names(sim$landscapeAttr), function(polygonType) {
   
    regime <- sim$scfmRegimePars[[polygonType]]
    landAttr <- sim$landscapeAttr[[polygonType]]
      
    #we know this table was produced with MinFireSize=2cells.
    
     y <- sim$cTable2$y #What are these supposed to be?
     x <- sim$cTable2$p 
     m.lw <- lowess(y~x,iter=2)
     if (sum(diff(m.lw$y)<0)>0)
        warning(sprintf("lowess curve non-monotone in zone %s. Proceed with caution", polygonType))
      targetSize <- regime$xBar/cellSize - 1 
      pJmp <- approx(m.lw$y,m.lw$x,targetSize,rule=2)$y
  
    w <- landAttr$nNbrs
    w <- w/sum(w)
    hatPE<-regime$pEscape
    if(hatPE == 0) { # no fires in polygon zone escapted
      foo <- 0 
    } else if (hatPE == 1) { # all fires in polygon zone escaped
      foo <- 1
    } else {
      res<-optimise(sim$escapeProbDelta,
                    interval=c(hatP0(hatPE,P(sim)$neighbours),
                               hatP0(hatPE,floor(sum(w*0:8)))),
                    tol=1e-4,
                    w=w,
                    hatPE=hatPE)
      foo <-res$minimum
      #It is almost obvious that the true minimum must occurr within the interval specified in the 
      #call to optimise, but I have not proved it, nor am I certain that the function being minimised is 
      #monotone.
    }
    #don't forget to scale by number of years, as well.
    rate <- regime$ignitionRate * cellSize #fireRegimeModel and this module must agree on 
                                                                   #an annual time step. How to test / enforce?
    pIgnition <- rate #approximate Poisson arrivals as a Bernoulli process at cell level.
                      #for Poisson rate << 1, the expected values are the same, partially accounting
                      #for multiple arrivals within years. Formerly, I used a poorer approximation
                      #where 1-p = P[x==0 | lambda=rate] (Armstrong and Cumming 2003).
    
    return(list(pSpread = pJmp,
                p0 = foo,
                naiveP0 = hatP0(regime$pEscape,8), 
                pIgnition = pIgnition,
                maxBurnCells = as.integer(round(regime$emfs/cellSize)))
    )
  })
  
  names(sim$scfmPars) <- names(sim$landscapeAttr)
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  
  if (!suppliedElsewhere("cTable2", sim)) {
    cTable2 <- read.csv(file.path(dataPath(sim), "FiresN1000MinFiresSize2NoLakes.csv"))
    sim$cTable2 <- cTable2
    }

  
  return(invisible(sim))
}

