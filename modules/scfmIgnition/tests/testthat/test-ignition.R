## scfmIgnition on a toy landscape: 10 x 10 pixels of 250 m. Column 5 is not flammable. Fire regime
## polygon 1 is columns 1-5, polygon 2 columns 6-10; `pIgnition` (from scfmDriver's calibration) is the
## per-pixel annual ignition probability of each polygon.

withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                          reproducible.verbose = -2))

toyRTM <- function() {
  terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2500, ymin = 0, ymax = 2500,
              crs = "EPSG:3005", vals = 1)
}
toyCol <- function(r) terra::colFromCell(r, seq_len(terra::ncell(r)))
vals <- function(r) as.vector(terra::values(r))

toySim <- function(pIgnition = c(1, 0)) {
  rtm <- toyRTM()
  flam <- rtm
  flam[toyCol(rtm) == 5] <- 0
  bbox <- function(xmin, xmax) {
    sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = 0, xmax = xmax, ymax = 2500), crs = sf::st_crs(3005)))
  }
  frp <- sf::st_sf(PolyID = 1:2, pIgnition = pIgnition, geometry = c(bbox(0, 1250), bbox(1250, 2500)))
  frr <- terra::rast(rtm, vals = ifelse(toyCol(rtm) <= 5, 1, 2))
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 1),
    modules = "scfmIgnition",
    objects = list(rasterToMatch = rtm, studyArea = sf::st_sf(geometry = bbox(0, 2500)),
                   fireRegimePolys = frp, fireRegimeRas = frr, flammableMap = flam),
    paths = list(modulePath = normalizePath(testthat::test_path("..", "..", "..")),
                 cachePath = root, inputPath = root, outputPath = root)
  ))
}

test_that("the ignition probability map takes each polygon's pIgnition, and NA where not flammable", {
  sim <- suppressMessages(SpaDES.core::spades(toySim(pIgnition = c(0.3, 0.02)), debug = FALSE))
  pIg <- vals(sim$pIg)
  col <- toyCol(toyRTM())
  expect_true(all(is.na(pIg[col == 5])))
  expect_identical(unique(pIg[col <= 4]), 0.3)
  expect_identical(unique(pIg[col >= 6]), 0.02)
})

test_that("ignitions fall only on flammable pixels, at their polygon's probability", {
  sim <- suppressMessages(SpaDES.core::spades(toySim(pIgnition = c(1, 0)), debug = FALSE))
  ## probability 1 in polygon 1 (its flammable pixels are columns 1-4), 0 in polygon 2
  expect_setequal(sim$ignitionLoci, which(toyCol(toyRTM()) <= 4))
  expect_false(anyDuplicated(sim$ignitionLoci) > 0)
})

test_that("each year draws new ignitions rather than adding to last year's", {
  sim <- toySim(pIgnition = c(1, 0))
  SpaDES.core::end(sim) <- 3
  sim <- suppressMessages(SpaDES.core::spades(sim, debug = FALSE))
  expect_length(sim$ignitionLoci, 40L)
})

test_that("a polygon with pIgnition 0 never ignites", {
  sim <- toySim(pIgnition = c(0, 0))
  SpaDES.core::end(sim) <- 5
  sim <- suppressMessages(SpaDES.core::spades(sim, debug = FALSE))
  expect_length(sim$ignitionLoci, 0L)
})
