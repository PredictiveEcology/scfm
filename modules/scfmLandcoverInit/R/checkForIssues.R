checkForIssues <- function(fireRegimePolys, studyArea, rasterToMatch, flammableMap, sliverThresh, cacheTag) {
  ## TODO: bug the group for the better function
  compareCRS(rasterToMatch, flammableMap, fireRegimePolys)

  if (is.null(fireRegimePolys$PolyID)) {
    #this is done in .inputObjects but should cover when the object
    stop("please supply fireRegimePolys with a PolyID")
  }

  if (sf::st_is_longlat(fireRegimePolys)) {
    stop("scfm requires projected coordinate systems - lat/long too prone to error")
  }
  ## TODO: fix all this rgeos business
  fireRegimePolys$trueArea <- round(sf::st_area(fireRegimePolys), digits = 0)

  if (any(as.numeric(fireRegimePolys$trueArea) < sliverThresh)) {
    message("sliver polygon(s) detected. Merging to their nearest valid neighbour")
    fireRegimePolys <- Cache(deSliver, fireRegimePolys, threshold = sliverThresh,
                                 userTags = cacheTag)
  }
  # #this is a problem if there is an upstream PROJ bug with gridded shapefiles...
  # if (length(unique(fireRegimePolys$PolyID)) != length(fireRegimePolys)) {
  #   stop("mismatch between PolyID and fireRegimePolys. Must be 1 PolyID value per multipolygon object")
  # }

  return(fireRegimePolys)
}
