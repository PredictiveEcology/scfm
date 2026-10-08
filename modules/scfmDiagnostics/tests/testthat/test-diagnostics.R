## scfmDiagnostics in "single" mode on a toy landscape: 20 x 20 pixels of 250 m (6.25 ha), all
## flammable, one fire regime polygon (2500 ha). The simulation is years 1-2, and `burnSummary` is what
## scfmSpread would record for three fires: one that did not escape (1 pixel) and two that did.

withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                          reproducible.verbose = -2, reproducible.useCache = FALSE))

toySim <- function() {
  rtm <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 5000, ymin = 0, ymax = 5000,
                     crs = "EPSG:3005", vals = 1L)
  sa <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 5000, ymax = 5000),
                                                       crs = sf::st_crs(3005))))
  frp <- sa
  frp$PolyID <- 1L
  frp$cellSize <- 6.25
  frp$ignitionRate <- 4e-4  # target: 4e-4 * 2500 ha = 1 ignition per year
  frp$pEscape <- 0.5
  frp$xBar <- 30
  frp$p0 <- 0.1
  frp$pSpread <- 0.23
  frp$empiricalBurnRate <- 0.01
  burnSummary <- data.table::data.table(
    igLoc = c(1L, 50L, 100L), grp = 1L, N = c(1, 8, 16), year = c(1, 1, 2),
    areaBurned = c(6.25, 50, 100), PolyID = 1L
  )
  ## historical fires: 1 ha (did not escape), 20 and 40 ha
  pts <- sf::st_sf(PolyID = 1L, SIZE_HA = c(1, 20, 40),
                   geometry = sf::st_sfc(lapply(c(1000, 2000, 3000), function(y) sf::st_point(c(2500, y))),
                                         crs = sf::st_crs(3005)))
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 2),
    modules = "scfmDiagnostics",
    params = list(scfmDiagnostics = list(mode = "single", .plots = NA)),
    objects = list(burnSummary = burnSummary, fireRegimePoints = pts, fireRegimePolys = frp,
                   flammableMap = rtm, studyAreaReporting = sa, burnMap = rtm),
    paths = list(modulePath = normalizePath(testthat::test_path("..", "..", "..")),
                 cachePath = root, inputPath = root, outputPath = root)
  ))
}

test_that("the summary compares the simulated fires with the fire regime's targets", {
  sim <- suppressWarnings(suppressMessages(SpaDES.core::spades(toySim(), debug = FALSE)))
  dt <- sim$scfmSummaryDT
  ## values per year of simulation carry the time unit from times(sim); compare the numbers
  expect_identical(nrow(dt), 1L)
  expect_equal(dt$burnableArea_ha, 2500)
  expect_equal(dt$targetIgnitions, 1)                  # ignitionRate * burnable area
  expect_equal(as.numeric(dt$achievedIgnitions), 3 / 2)            # 3 fires in 2 years
  expect_equal(dt$targetEscapes, 0.5)                  # pEscape * targetIgnitions
  expect_equal(as.numeric(dt$achievedEscapes), 2 / 2)              # 2 fires larger than one pixel in 2 years
  expect_equal(dt$modMeanSize, mean(c(50, 100)))       # mean size of the escaped fires
  expect_equal(dt$histMeanSize, 30)                    # the regime's xBar
  expect_equal(dt$histMedianSize, median(c(20, 40)))   # historical fires larger than one pixel
  expect_equal(as.numeric(dt$achievedFRI), 2 / ((50 + 100) / 2500)) # years / (area burned / burnable area)
  expect_equal(dt$targetFRI, 1 / 0.01)                 # 1 / empiricalBurnRate
})

test_that("the summary is also written to the output folder", {
  sim <- suppressWarnings(suppressMessages(SpaDES.core::spades(toySim(), debug = FALSE)))
  f <- file.path(SpaDES.core::outputPath(sim), "scfmDiagnostics_single_summary_dt.csv")
  expect_true(file.exists(f))
  expect_equal(read.csv(f)$achievedIgnitions, 3 / 2)
})
