ratioPartition2 <- function (targetBurnRate, #
                             empiricalBurnRate, #rate estimated by the model
                             pEscape, #escape probability
                             xBar, # mean fire size
                             rate){ # iginition rate

  ratio <-  targetBurnRate / empiricalBurnRate 
  
  remains <- ratio
#   
#   step <- function(remains){
#     step = min(remains, 2)
#   } 
#   
#   remains <- function(remains, step){
#     remains <- remains / step
#     print(remains)
#   }
#   
#   while (step > 2){
#     pEscape <- pEscape * step
#     remains(remains, step)
#     step <- step (remains)
#     if (step > 2){
#       xBar <- xBar * step
#       remains(remains, step)}
#   }
# return(list(
#   rate = rate,
#   pEscape = pEscape,
#   xBar = xBar  #mean fire size
# ))   
# }
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

return(list(
  rate = rate,
  pEscape = pEscape,
  xBar = xBar  #mean fire size
))
}