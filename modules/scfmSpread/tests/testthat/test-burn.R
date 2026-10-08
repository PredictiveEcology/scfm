## scfmSpread's burn event, fed by scfmEscape, on a toy landscape. `rstCurrentBurn` is this year's burn:
## Biomass_regeneration regenerates wherever it is > 0 and burn summaries count it, so it must hold this
## year's burned pixels and no others.
##
## The landscape is 10 x 10 pixels of 250 m (6.25 ha). Column 5 is not flammable, so with spread
## probability 1 an escaped fire burns exactly its own block: cols 1-4 ("west") or cols 6-10 ("east").
## Cell 1 is in the west block, cell 10 in the east block.

withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                          reproducible.verbose = -2))

toyCol <- function(r) terra::colFromCell(r, seq_len(terra::ncell(r)))
toyRTM <- function() {
  terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2500, ymin = 0, ymax = 2500,
              crs = "EPSG:3005", vals = 1)
}
westCells <- function() which(toyCol(toyRTM()) <= 4)
eastCells <- function() which(toyCol(toyRTM()) >= 6)
barrierCells <- function() which(toyCol(toyRTM()) == 5)
vals <- function(r) as.vector(terra::values(r))

## simInit() of scfmEscape + scfmSpread for years 1..1; burnYear() runs one year at a time.
## fireRegimePolys gives p0 = 1 (every ignition escapes) and pSpread = 1. With `maxBurnCells` (two
## values) there are two fire regime polygons, PolyID 1 (cols 1-5) and PolyID 2 (cols 6-10).
toySim <- function(maxBurnCells = NULL) {
  rtm <- toyRTM()
  flam <- rtm
  flam[toyCol(rtm) == 5] <- 0
  bbox <- function(xmin, xmax) {
    sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = 0, xmax = xmax, ymax = 2500), crs = sf::st_crs(3005)))
  }
  sa <- sf::st_sf(geometry = bbox(0, 2500))
  frp <- sa
  frp$PolyID <- 1L
  frr <- rtm
  if (!is.null(maxBurnCells)) {
    frp <- sf::st_sf(PolyID = 1:2, maxBurnCells = maxBurnCells, geometry = c(bbox(0, 1250), bbox(1250, 2500)))
    frr[toyCol(rtm) >= 6] <- 2
  }
  frp$cellSize <- 6.25
  frp$p0 <- 1
  frp$pSpread <- 1
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 1, end = 1),
    modules = c("scfmEscape", "scfmSpread"),
    params = list(scfmEscape = list(.useCache = FALSE),
                  scfmSpread = list(.useCache = FALSE, .plots = NA)),
    objects = list(rasterToMatch = rtm, studyArea = sa, studyAreaReporting = sa,
                   fireRegimePolys = frp, fireRegimeRas = frr, flammableMap = flam,
                   ignitionLoci = integer(0)),
    paths = list(modulePath = normalizePath(testthat::test_path("..", "..", "..")),
                 cachePath = root, inputPath = root, outputPath = root)
  ))
}

## Run `year` with these ignitions (what scfmIgnition would give) and, optionally, this escape probability.
burnYear <- function(sim, year, ignitionLoci, p0 = NULL) {
  sim$ignitionLoci <- ignitionLoci
  if (!is.null(p0)) sim$p0 <- p0
  SpaDES.core::end(sim) <- year
  suppressMessages(SpaDES.core::spades(sim, debug = FALSE))
}

test_that("ignitions that do not escape burn their own pixel in a year with escapes", {
  ## the rule the no-escape year follows: spread2() keeps scfmEscape's ignition pixels in burnDT
  sim <- burnYear(toySim(), 1, integer(0)) # year 1 runs the init events, which set p0 from fireRegimePolys
  p0 <- toyRTM()
  p0[] <- ifelse(toyCol(p0) <= 4, 1, 0) # west ignitions escape, east ones do not
  sim <- burnYear(sim, 2, c(1L, 10L), p0 = p0)
  expect_identical(which(vals(sim$rstCurrentBurn) == 1), sort(c(westCells(), 10L)))
  expect_identical(vals(sim$burnMap)[10], 1)
})

test_that("a year whose ignitions do not escape burns only those pixels, not last year's", {
  sim <- burnYear(toySim(), 1, 1L)
  expect_identical(which(vals(sim$rstCurrentBurn) == 1), westCells())

  sim <- burnYear(sim, 2, 10L, p0 = 0) # one ignition in the east block; nothing escapes
  cur <- vals(sim$rstCurrentBurn)
  expect_identical(which(cur == 1), 10L)
  expect_identical(cur[westCells()], rep(0, 40))                # flammable, not burned this year
  expect_true(all(is.na(cur[barrierCells()])))
  bm <- vals(sim$burnMap)
  expect_identical(bm[westCells()], rep(1, 40))                 # year 1 only
  expect_identical(bm[eastCells()], as.numeric(eastCells() == 10L)) # year 2's ignition pixel only
  expect_identical(sim$burnSummary[year == 2]$igLoc, 10L)
  expect_equal(sim$burnSummary[year == 2]$N, 1)
  ## the ignition pixel burned, so its time since fire restarts, as in a year with escapes
  expect_identical(vals(sim$timeSinceFire)[10], 0)
})

test_that("a year with no ignitions burns nothing and does not spread last year's fires again", {
  for (noIgnitions in list(integer(0), NULL)) { # scfmIgnition gives integer(0); its init gives NULL
    sim <- burnYear(toySim(), 1, 1L)
    burnMap1 <- vals(sim$burnMap)
    nSummary1 <- NROW(sim$burnSummary)

    sim <- burnYear(sim, 2, noIgnitions)
    expect_null(sim$spreadState)
    expect_false(any(vals(sim$rstCurrentBurn) == 1, na.rm = TRUE))
    expect_identical(vals(sim$rstCurrentBurn)[westCells()], rep(0, 40))
    expect_identical(vals(sim$burnMap), burnMap1)
    expect_identical(NROW(sim$burnSummary), nSummary1)
  }
})

test_that("time since fire ages by 1 in a year with no ignitions", {
  sim <- burnYear(toySim(), 1, 1L)
  expect_identical(vals(sim$timeSinceFire)[westCells()], rep(0, 40))
  sim <- burnYear(sim, 2, integer(0))
  tsf <- vals(sim$timeSinceFire)
  expect_identical(tsf[westCells()], rep(1, 40))
  expect_true(all(is.na(tsf[eastCells()]))) # never burned, and not supplied: stays NA
})

test_that("with maxBurnCells no fire burns more pixels than its fire regime polygon's cap", {
  ## Caps well below the block sizes (40 west, 50 east), so with spread probability 1 each fire burns
  ## exactly its cap. The east ignition comes first: caps must follow the fire, not the ignition order.
  ## Each corner ignition escapes to its 3 neighbours: a cap of 2 is reached during the escape.
  for (caps in list(c(6L, 12L), c(2L, 12L))) {
    sim <- burnYear(toySim(maxBurnCells = caps), 1, c(10L, 1L))
    bs <- sim$burnSummary[year == 1]
    expect_equal(bs$N[match(c(1L, 10L), bs$igLoc)], caps) # west fire (cell 1), east fire (cell 10)
    expect_identical(sum(vals(sim$rstCurrentBurn) == 1, na.rm = TRUE), sum(caps))
  }
})
