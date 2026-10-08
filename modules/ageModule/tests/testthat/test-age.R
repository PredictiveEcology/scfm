## ageModule's age event on a toy landscape: 10 x 10 pixels of 250 m.

withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                          reproducible.verbose = -2))

toyRTM <- function() {
  terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2500, ymin = 0, ymax = 2500,
              crs = "EPSG:3005", vals = 1)
}

## simInit() of ageModule alone, starting at year 1. `ages` is the initial ageMap; `burned` the
## pixels of rstCurrentBurn that are 1 (what scfmSpread gives for this year's fires).
toySim <- function(ages, burned = integer(0), maxAge = 200) {
  rtm <- toyRTM()
  age <- terra::rast(rtm, vals = ages)
  burn <- terra::rast(rtm, vals = 0)
  burn[burned] <- 1
  sa <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2500, ymax = 2500),
                                                       crs = sf::st_crs(3005))))
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 1),
    modules = "ageModule",
    params = list(ageModule = list(maxAge = maxAge, .plots = NA)),
    objects = list(ageMap = age, rasterToMatch = rtm, rstCurrentBurn = burn, studyArea = sa),
    paths = list(modulePath = normalizePath(testthat::test_path("..", "..", "..")),
                 cachePath = root, inputPath = root, outputPath = root)
  ))
}

vals <- function(r) as.vector(terra::values(r))

test_that("each year, unburned pixels age by one year and burned pixels go back to 0", {
  sim <- suppressMessages(SpaDES.core::spades(toySim(ages = 10, burned = 1:5), debug = FALSE))
  age <- vals(sim$ageMap)
  expect_identical(unique(age[1:5]), 0)
  expect_identical(unique(age[-(1:5)]), 11)
})

test_that("ages stop at maxAge", {
  sim <- suppressMessages(SpaDES.core::spades(toySim(ages = c(rep(199, 50), rep(200, 50)), maxAge = 200),
                                              debug = FALSE))
  expect_identical(unique(vals(sim$ageMap)), 200)
})

test_that("ages keep increasing over several years", {
  sim <- toySim(ages = 10)
  SpaDES.core::end(sim) <- 3
  sim <- suppressMessages(SpaDES.core::spades(sim, debug = FALSE))
  expect_identical(unique(vals(sim$ageMap)), 13)
})
