## scfmEscape's escape event on a toy landscape: 10 x 10 pixels of 250 m, all flammable.

withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                          reproducible.verbose = -2))

## simInit() of scfmEscape alone for years 1..1, with escape probability p0 = 1.
toySim <- function() {
  rtm <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2500, ymin = 0, ymax = 2500,
                     crs = "EPSG:3005", vals = 1)
  sa <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2500, ymax = 2500),
                                                       crs = sf::st_crs(3005))))
  frp <- sa
  frp$PolyID <- 1L
  frp$p0 <- 1
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 1),
    modules = "scfmEscape",
    params = list(scfmEscape = list(.useCache = FALSE)),
    objects = list(rasterToMatch = rtm, studyArea = sa, fireRegimePolys = frp,
                   fireRegimeRas = rtm, flammableMap = rtm, ignitionLoci = integer(0)),
    paths = list(modulePath = normalizePath(testthat::test_path("..", "..", "..")),
                 cachePath = root, inputPath = root, outputPath = root)
  ))
}

escapeYear <- function(sim, year, ignitionLoci) {
  sim$ignitionLoci <- ignitionLoci
  SpaDES.core::end(sim) <- year
  suppressMessages(SpaDES.core::spades(sim, debug = FALSE))
}

test_that("a year with no ignitions leaves no spreadState from the year before", {
  for (noIgnitions in list(integer(0), NULL)) { # scfmIgnition gives integer(0); its init gives NULL
    sim <- escapeYear(toySim(), 1, 45L)
    expect_identical(sort(unique(sim$spreadState$initialPixels)), 45L)
    expect_gt(NROW(sim$spreadState[state == "activeSource"]), 0) # p0 = 1: it escaped

    sim <- escapeYear(sim, 2, noIgnitions)
    expect_null(sim$spreadState) # scfmSpread would otherwise spread year 1's fire again
  }
})
