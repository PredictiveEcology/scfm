defineModule(sim, list(
  name="andisonDriver",
  description="make scfm fire modules simulate Andisons FRIs by fudging scfmRegime parameters.",
  keywords=c("fire"),
  authors=c(person(c("Steven", "G"), "Cumming", email="stevec@sbf.ulaval.ca", role=c("aut", "cre"))),
  childModules=character(),
  version=numeric_version("0.1.0"),
  spatialExtent=raster::extent(rep(NA_real_, 4)),
  timeframe=as.POSIXlt(c(NA, NA)),
  timeunit="year", 
  citation=list(),
  documentation = list("README.txt", "andisonDriver.Rmd"),
  reqdPkgs=list("stats"), 
  parameters=rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter("neighbours", "numeric", 8, 4, 8, "number of cell immediate neighbours")),
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
  inputObjects = bind_rows(
    expectsInput(objectName = "scfmRegimePars", objectClass = "list", desc = ""),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = ""),
    expectsInput(objectName = "cTable2", objectClass = "data.frame", 
                 desc = "A csv containing results of fire experiment", 
                 sourceURL = "https://drive.google.com/open?id=155fOsdEJUJNX0yAO_82YpQeS2-bA1KGd"),
    expectsInput(objectName = "AndisonFRI", objectClass = "spatialPolygonsDataFrame", desc = "Dave's FRI data", 
                 sourceURL = "https://drive.google.com/file/d/1JptU0R7qsHOEAEkxybx5MGg650KC98c6/view?usp=sharing")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName="scfmPars", objectClass = "list", desc = ""),
    createsOutput(objectName = "meanFRI", objectClass = "lsit", desc = "list of mean FRI over studyArea polygons")
  )                      
))

## event types
#   - type `init` is required for initiliazation

doEvent.andisonDriver = function(sim, eventTime, eventType, debug = FALSE) {
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

odds <- function(p)p/(1-p)

oddsRatio <- function(p,q) odds(p)/odds(q)

calcp <-function(q,or){ # given q and or(p,q), solve for p
  x <- or * odds(q)
  return ( 1 / ((1 / x) +1))
}

Init <- function(sim) {
  
  cellSize <- sim$landscapeAttr[[1]]$cellSize
  
  sim$scfmPars<- lapply(names(sim$landscapeAttr), function(polygonType) {
    
    regime <- sim$scfmRegimePars[[polygonType]]
    landAttr <- sim$landscapeAttr[[polygonType]]
    browser()
    fri <- sim$  #will return LTHFC FIX THIS LINE
    fri <- ifelse(fri < 20, 20, fri)
    targetAAB <- landAttr$burnyArea / fri
    
    rate <- regime$ignitionRate
    pEscape <- regime$pEscape
    mfs <- regime$xBar
    
    scfmAAB <- rate * landAttr$burnyArea * pEscape * mfs
    ratio <- targetAAB/scfmAAB
    
    if (ratio < 0.05){
        warning(sprintf("AAB ratio %5.3f in %s: no action taken",ratio, polygonType))
    }
    else{
    if (ratio > 1.05) { #we have work to do
       or <- oddsRatio(0.179,0.111)
       pEscape <- calcP(pEscape,or)
       scfmAAB <- rate * landAttr$burnyArea * pEscape * mfs
       ratio <- targetAAB/scfmAAB
       #caution this step could over-correct
    }
    if (ratio > 1.05){ #we have more work to do}
         or <- oddsRatio(0.319,0.179)
         pEscape <- calcP(pEscape,or)
         scfmAAB <- rate * landAttr$burnyArea * pEscape * mfs
         ratio <- targetAAB/scfmAAB
         #caution this step could over-correct, too
    }
    if (ratio < 0.95){
    #we got this way by dialing up pEscape, so dial it back down again, and we're apples.
      pEscape <- pEscape * ratio
      scfmAAB <- rate * landAttr$burnyArea * pEscape * mfs
      ratio <- targetAAB/scfmAAB
    }
    if (ratio > 1.05){
      mfs <- mfs * max(ratio, 2)
      scfmAAB <- rate * landAttr$burnyArea * pEscape * mfs
      ratio <- targetAAB/scfmAAB
    }
    if (ratio > 1.05){
      rate <- rate * max(ratio, 2)
      scfmAAB <- rate * landAttr$burnyArea * pEscape * mfs
      ratio <- targetAAB/scfmAAB
    }
    if (ratio > 1.05)
      warning(sprintf("Target FRI of %03d in zone %s not achievable from data",fri, polygonType))
    }
    
    
    #this table contains calibration data for several landscape sizes
    #and several min fire sizes (1 or 2 cells), organised by collumn.
    #The data were made by Steve Cumming in June 2013 for a whole other purpose.
    #I chose the one that seems most appropriate to me
    #we know this table was produced with MinFireSize=2cells.
    
     y <- sim$cTable2$y 
     x <- sim$cTable2$p 
     m.lw <- lowess(y~x,iter=2)
     if (sum(diff(m.lw$y)<0)>0)
        warning(sprintf("lowess curve non-monotone in zone %s. Proceed with caution", polygonType))
      targetSize <- mfs/cellSize - 1 
      pJmp <- approx(m.lw$y,m.lw$x,targetSize,rule=2)$y
  
    w <- landAttr$nNbrs
    w <- w/sum(w)
    hatPE <- pEscape
    if(hatPE == 0) { # no fires in polygon zone escapted
      foo <- 0 
    } else if (hatPE == 1) { # all fires in polygon zone escaped
      foo <- 1
    } else {
      res <- optimise(escapeProbDelta,
                    interval=c(hatP0(hatPE,P(sim)$neighbours),
                               hatP0(hatPE,floor(sum(w*0:8)))),
                    tol=1e-4,
                    w=w,
                    hatPE=hatPE)
      foo <-res[["minimum"]]
      #It is almost obvious that the true minimum must occurr within the interval specified in the 
      #call to optimise, but I have not proved it, nor am I certain that the function being minimised is 
      #monotone.
    }
    #don't forget to scale by number of years, as well.
    rate <- rate * cellSize #fireRegimeModel and this module must agree on 
                                                                   #an annual time step. How to test / enforce?
    pIgnition <- rate #approximate Poisson arrivals as a Bernoulli process at cell level.
                      #for Poisson rate << 1, the expected values are the same, partially accounting
                      #for multiple arrivals within years. Formerly, I used a poorer approximation
                      #where 1-p = P[x==0 | lambda=rate] (Armstrong and Cumming 2003).
    
    return(list(pSpread = pJmp,
                p0 = foo,
                naiveP0 = hatP0(pEscape,8), 
                pIgnition = pIgnition,
                maxBurnCells = as.integer(round(regime$emfs/cellSize)))
    )
  })
  
  names(sim$scfmPars) <- names(sim$landscapeAttr)
  
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  
  dPath <- dataPath(sim)
  if (!suppliedElsewhere("cTable2", sim)) {
    cTable2 <- Cache(prepInputs,
                     targetFile = file.path(dPath, "FiresN1000MinFiresSize2NoLakes.csv"),
                     url = extractURL(objectName = "cTable2", sim), 
                     fun = "read.csv",
                     destinationPath = dPath,
                     overwrite = TRUE,
                     filename2 = TRUE)
    
    
    sim$cTable2 <- cTable2
  }
  return(invisible(sim))
}