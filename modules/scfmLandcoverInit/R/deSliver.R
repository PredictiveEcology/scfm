#' A function for merging sliver polygons into non-sliver neighbours.
#' The threshold is applied to the area of the multipolygon object, not each
#' individual polygon. Non-sliver polygons keep their original attributesl
#' Inteded to be used when it is important to retain the original extent of an
#' area while removing sliver polygons
#' @keywords sliver polygons intersect
#' @param x A spatialPolygonsDataFrame or sf object
#' @param threshold the minimum area below which a polygon is considered a sliver
#' @return an object of class sf or Spatial with sliver polygons
#' merged to their nearest valid neighbour.
#'
#' @export
#' @import sf #you basically need the whole package at thsi point
#' @examples
#'deSliver(x = intersectedPolygons, threshold = 500)

deSliver <- function(x, threshold) {
  x$tempArea <- as.numeric(st_area(x))

  #determine slivers by area
  xSlivers <- x[x$tempArea < threshold, ]
  xNotSlivers <- x[x$tempArea > threshold, ]
  if (nrow(xNotSlivers) < 1) {
    stop("Threshold exceeds the area of every polygon. Please select a smaller number")
  }

  #Split slivers from multipolygon, or nearest feature may be incorrect
  xSlivers <- suppressWarnings(st_cast(xSlivers, "POLYGON"))

  #Find nearest non-sliver
  nearestFeature <- st_nearest_feature(xSlivers, xNotSlivers)

  #Merge each sliver polygon into nearest neighbour
  mergeSlivers <- lapply(
    unique(nearestFeature),
    FUN = function(i,
                   ns = xNotSlivers,
                   s = xSlivers,
                   nf = nearestFeature) {
      featurePolys <- nearestFeature == i
      xMerge <- s[featurePolys, ] %>%
        st_union(.)
      yMerge <- ns[i, ]
      #convert slivers back to multipolygon
      out <- sf::st_union(x = xMerge, y = yMerge)

      #update the geometry
      yMerge$geometry <- out
      return(yMerge)
    }
  )
  otherPolys <- xNotSlivers[!(1:nrow(xNotSlivers) %in% nearestFeature),]
  if (length(mergeSlivers) > 1) {

    #these polygons must be tracked and merged. They may be nrow(0) if every feature was modified in some way
    if (nrow(otherPolys) != 0) {
      mergeSlivers <- do.call(rbind, mergeSlivers)
      m <- rbind(otherPolys, mergeSlivers)
    } else {
      m <- rbind(mergeSlivers)
    }
  } else { #lapply over length 1 is special
    if (nrow(otherPolys) != 0) {
      mergeSlivers <- mergeSlivers[[1]]
      m <- rbind(mergeSlivers, otherPolys)
    } else {
      m <- mergeSlivers[[1]]
    }
  }
  #the geometry will be sfc
  m <- st_cast(m, to = "MULTIPOLYGON")
  #Remove the temporary column
  m$tempArea <- NULL

  #remove self-intersecting geometries
  if (any(!st_is_valid(m))) {
    m <- gBuffer(spgeom = m, byid = TRUE, width = 0)
  }

  return(m)
}
