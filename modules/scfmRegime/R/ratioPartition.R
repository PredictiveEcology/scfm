#library("partitions")
ratioPartition <- function (targetBurnRate, #
                          empiricalBurnRate, #rate estimated by the model
                          pEscape, #escape probability
                          xBar, # mean fire size
                          rate){ # iginition rate
  browser()


    ratio <-  targetBurnRate / empiricalBurnRate 


## How do I partition this 'ratio' among 3 processes: same increase in all processes for now.
pEscapeRatio <- ratio^(1/3)
#ignitionRatio <- ratio^(1/3)
#sizeRatio <- ratio^(1/3)
#ratioPart <-partitions::restrictedparts(pEscapeRatio, 3) ### only works with small numbers 
switch <- FALSE
if (pEscapeRatio <= 2){
  pEscapeTemp <- pEscape * 2
  pEscapeRatioNew <- pEscapeRatio / 2 
  # Evaluate pEscape for being > 1
  if (pEscapeTemp > 1){
    warning("careful, pEscape cannot be greater than one")
    # If > 1, keep pEscape at 1 and switch the function 
    pEscapeTemp <- 1
    # we switch to 2
    switch <- TRUE
  }
 } 
if (any(switch, all(pEscapeRatio > 2, pEscapeRatio <= 4, pEscape <= 1))){
  pEscapeTemp <- pEscape * 4
  pEscapeRatioNew <- pEscapeRatio / 4 
  } 
if (all(pEscapeRatio > 4, pEscape <= 1)){  
  pEscapeTemp <- pEscape * 4 
  xBar <- xBar * 2
  pEscapeRatioNew <- pEscapeRatio / 6
} else {
    warning("pEscape bigger than 1")
  }
if (all(pEscapeRatioNew > 4, pEscape <= 1)){
  pEscapeTemp <- pEscape * 4
  xBar <- xBar * 2
  rate <- rate * 2
  
  counter <- 0
  while (all(pEscapeRatioNew > 4, pEscape <= 1)){
    xBar <- xBar * (2+ counter)
    rate <- rate *(2 + counter)
    pEscapeRatioNew <- pEscapeRatio / (8 + counter)
    counter <- counter + 1
  }
   
}
return(list(
  rate = rate,
  pEscape = pEscapeTemp,
  xBar = xBar  #mean fire size
))
}
