## scfmDataPrep's land-cover and fire-regime steps on a toy landscape, with every input supplied so
## nothing is downloaded: 20 x 20 pixels of 250 m (6.25 ha). Fire regime polygon 1 is columns 1-10,
## polygon 2 columns 11-20. Column 1 is not flammable, so polygon 1 has 200 - 20 = 180 flammable
## pixels and polygon 2 has 200. The calibration area is the study area.

withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                          reproducible.verbose = -2, reproducible.useCache = FALSE))

toyRTM <- function() {
  terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 5000, ymin = 0, ymax = 5000,
              crs = "EPSG:3005", vals = 1L)
}
toyCol <- function(r) terra::colFromCell(r, seq_len(terra::ncell(r)))
bbox <- function(xmin, xmax) {
  sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = 0, xmax = xmax, ymax = 5000), crs = sf::st_crs(3005)))
}

## Fire points (NFDB-like): in polygon 1, four lightning fires in the epoch (2000-2009) of 1, 2, 50 and
## 100 ha, plus one outside the epoch and one human-caused, which must both be left out.
toyFirePoints <- function() {
  sf::st_sf(
    CAUSE = c("L", "L", "L", "L", "L", "H"),
    YEAR = c(2001, 2003, 2005, 2009, 1990, 2005),
    SIZE_HA = c(1, 2, 50, 100, 500, 300),
    geometry = sf::st_sfc(lapply(c(1000, 1200, 1400, 1600, 1800, 2000), function(y) sf::st_point(c(1500, y))),
                          crs = sf::st_crs(3005))
  )
}

toySim <- function(eventsToPrepare = c("scfmLandcoverInit", "scfmRegime"), ...) {
  rtm <- toyRTM()
  flam <- rtm
  flam[toyCol(rtm) == 1] <- 0L
  sa <- sf::st_sf(geometry = bbox(0, 5000))
  frp <- sf::st_sf(PolyID = 1:2, geometry = c(bbox(0, 2500), bbox(2500, 5000)))
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 1),
    modules = "scfmDataPrep",
    params = list(scfmDataPrep = list(eventsToPrepare = eventsToPrepare, fireEpoch = c(2000, 2009),
                                      sliverThreshold = 1, .plots = NA, .useCache = FALSE, ...)),
    objects = list(studyArea = sa, studyAreaCalibration = sa,
                   rasterToMatch = rtm, rasterToMatchCalibration = rtm,
                   flammableMap = flam, flammableMapCalibration = flam,
                   fireRegimePolys = frp, fireRegimePolysCalibration = frp,
                   firePoints = toyFirePoints()),
    paths = list(modulePath = normalizePath(testthat::test_path("..", "..", "..")),
                 cachePath = root, inputPath = root, outputPath = root)
  ))
}

runInit <- function(sim) suppressWarnings(suppressMessages(SpaDES.core::spades(sim, debug = FALSE)))

test_that("each fire regime polygon gets its flammable area", {
  frp <- runInit(toySim())$fireRegimePolys
  frp <- frp[order(frp$PolyID), ]
  expect_identical(as.numeric(frp$nFlammable), c(180, 200))
  expect_identical(unique(frp$cellSize), 6.25)
  expect_identical(frp$burnyArea, c(180, 200) * 6.25)
})

test_that("the fire regime map holds each polygon's PolyID", {
  sim <- runInit(toySim())
  expect_equal(as.vector(terra::values(sim$fireRegimeRas)), ifelse(toyCol(toyRTM()) <= 10, 1, 2))
})

test_that("only lightning fires inside the epoch are used", {
  pts <- runInit(toySim())$fireRegimePoints
  expect_setequal(pts$SIZE_HA, c(1, 2, 50, 100))
  expect_identical(unique(pts$PolyID), 1L)
})

test_that("ignition rate, escape probability, mean fire size and burn rate follow from those fires", {
  frp <- runInit(toySim())$fireRegimePolys
  p1 <- frp[frp$PolyID == 1, ]
  burnyArea <- 180 * 6.25
  expect_equal(p1$ignitionRate, 4 / (10 * burnyArea))               # fires per ha per year
  expect_equal(p1$pEscape, 2 / 4)                                    # 50 and 100 ha exceed one pixel
  expect_equal(p1$xBar, mean(c(50, 100)))                            # mean size of escaped fires
  expect_equal(p1$xMax, 100)
  expect_equal(p1$empiricalBurnRate, (1 + 2 + 50 + 100) / (10 * burnyArea))
  ## no fires in polygon 2: no estimates
  expect_true(is.na(frp$ignitionRate[frp$PolyID == 2]))
})

test_that("the spread calibration gives each polygon with fires its spread, escape and ignition probabilities", {
  ## a small targetN keeps the calibration to seconds; the bounds, not the fitted value, are under test
  sim <- runInit(toySim(eventsToPrepare = c("scfmLandcoverInit", "scfmRegime", "scfmDriver"),
                        targetN = 100, .useParallelFireRegimePolys = FALSE))
  frp <- sim$fireRegimePolys
  p1 <- frp[frp$PolyID == 1, ]
  expect_gte(p1$pSpread, 0.185) # pMin
  expect_lte(p1$pSpread, 0.253) # pMax
  expect_gt(p1$p0, 0)
  expect_lte(p1$p0, 1)
  ## per-pixel annual ignition probability: ignitions per ha per year times the pixel's area
  expect_equal(p1$pIgnition, p1$ignitionRate * p1$cellSize)
  expect_equal(p1$maxBurnCells, round(p1$emfs_ha / p1$cellSize))
  ## polygon 2 had no fires: it never ignites
  expect_identical(frp$pIgnition[frp$PolyID == 2], 0)
})
