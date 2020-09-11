ratioPartitionSC <- function (targetBurnRate, #
                            empiricalBurnRate, #rate estimated by the model
                            pEscape, #escape probability
                            xBar, # mean fire size
                            rate){ # iginition rate
  browser()
  
  
  ratio <-  targetBurnRate / empiricalBurnRate 
 
  remains <- ratio
  
  step <- min(remains, 2)
  pEscape <- pEscape * step
  remains <- remains / step
  step <- min(remains, 2)
  xBar <- xBar * step
  remains <- remains / step
  step <- min(remains, 2)
  pEscape <- pEscape * step 
  remains <- remains/step
  step <- min(remains, 2)
  xBar <- xBar * step
  remains <- remains / step
  step <- min(remains, 2)
  rate <- rate  * step
  remains <- remains / step
  xBar <- xBar * remains
}