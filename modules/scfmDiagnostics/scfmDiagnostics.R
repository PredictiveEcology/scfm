defineModule(sim, list(
  name = "scfmDiagnostics",
  description = paste("Diagnostic plots for scfm. Can be run for a single simulation,",
                      "as part of the main scfm run, or as part of postprocessing.",
                      "Inputs objects will be loaded from saved simulation files when in 'multi' mode."),
  keywords = c("diagnostics", "scfm"),
  authors = c(
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre"))
  ),
  childModules = character(0),
  version = list(scfmDiagnostics = "2.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmDiagnostics.Rmd"), ## same file
  loadOrder = list(after = c("scfmLandcoverInit", "scfmRegime", "scfmDriver",
                             "scfmEscape", "scfmIgnition", "scfmSpread")),
  reqdPkgs = list(
    "ggplot2", "gridExtra",
    "PredictiveEcology/reproducible@development (>= 2.1.0)",
    "PredictiveEcology/scfmutils (>= 2.0.1)",
    "PredictiveEcology/SpaDES.core@development (>= 2.1.0.9005)"
  ),
  parameters = bindrows(
    defineParameter("mode", "character", "single", NA, NA,
                    paste("use 'single' to run part of an scfm simulation (i.e., along with other scfm modules);",
                          "use 'multi' to run as part of postprocessing multiple scfm runs.")),
    defineParameter("reps", "integer", NA_integer_, 1L, NA_integer_,
                    paste("number of replicates/runs per study area when running in 'multi' mode.")),
    defineParameter("simOutPrefix", "character", "mySimOut", NA_character_, NA_character_,
                    "saved simList file prefix"),
    defineParameter("simTimes", "numeric", c(NA, NA), NA, NA,
                    "Simulation start and end times when running in 'multi' mode."),
    defineParameter(".plots", "character", c("screen", "png"), NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                    "area obtained using `reproducible::studyAreaName()`"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput("burnSummary", "data.table",
                 "describes details of all burned pixels. Required in single mode.",
                 sourceURL = NA),
    expectsInput("burnMap", "SpatRaster",
                 "cumulative burn map from simulation",
                 sourceURL = NA),
    expectsInput("fireRegimePoints", "sf",
                 "Fire locations. Points outside `studyArea` are removed. Required in single mode.",
                 sourceURL = NA),
    expectsInput("fireRegimePolys", "sf",
                 paste("Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
                       "Must have numeric field 'PolyID' or it will be created for individual polygons.",
                       "Required in single mode."),
                 sourceURL = NA),
    expectsInput("flammableMap", "SpatRaster",
                 "binary flammability map. Required in single mode.",
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "sf",
                 paste("multipolygon (typically smaller/unbuffered than `studyArea`) to use for plotting/reporting.",
                       "Required in single mode."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("scfmSummaryDT", "data.table",
                  desc = paste("Summary data.table containing diagnostic plot data.",
                               "Can be used to create customized diagnostic plots;",
                               "see `?scfmutils::comparePredictions`."))
  )
))

doEvent.scfmDiagnostics = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (P(sim)$mode == "single") {
        sim <- scheduleEvent(sim, end(sim), "scfmDiagnostics", "diagnosticPlotsSingle", eventPriority = .last())
      } else if (P(sim)$mode == "multi") {
        sim <- scheduleEvent(sim, start(sim), "scfmDiagnostics", "diagnosticPlotsMulti")
      }
    },
    diagnosticPlotsSingle = {

      dt <- diagnosticPlotsDT(sim)

      write.csv(dt, file.path(outputPath(sim), "scfmDiagnostics_single_summary_dt.csv"))

      ## Some useful plots
      gg_fri <- scfmutils::comparePredictions_fireReturnInterval(dt, times = times(sim))
      gg_frp <- scfmutils::plot_fireRegimePolys(sim$fireRegimePolys)
      gg_ign <- scfmutils::comparePredictions_annualIgnitions(dt)
      gg_mfs <- scfmutils::comparePredictions_meanFireSize(dt)
      gg_esc <- scfmutils::comparePredictions_annualEscapes(dt)
      ## NOTE: historical distribution is derived purely from historical data
      gg_histDist <- scfmutils::comparePredictions_fireDistribution(
        sim$fireRegimePoints,
        size = min(sim$fireRegimePolys$cellSize),
        burnSummary = sim$burnSummary
      )
      ## note that fireRegimePoints may include SAL but this figure only compares distribution
      ## so total area is irrelevant

      ## removed MAAB as diagnostic plot because it was derived from fire points incorrectly when
      ## studyAreaCalibration is supplied; MAAB can still be calculated manually if needed by user ## TODO

      if ("png" %in% P(sim)$.plots) {
        ggsave(file.path(figurePath(sim), "FRI.png"), gg_fri, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "FRP.png"), gg_frp, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "IGN.png"), gg_ign, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "MFS.png"), gg_mfs, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "ESC.png"), gg_esc, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "histDist.png"), gg_histDist, height = 8, width = 8)
      }

      if ("screen" %in% P(sim)$.plots) {
        clearPlot()
        gridExtra::grid.arrange(gg_frp,  gg_mfs,  gg_fri,  gg_ign, gg_esc, gg_histDist,
                                nrow = 2, ncol = 3)
      }

      sim$scfmSummaryDT <- dt
    },
    diagnosticPlotsMulti = {
      allReps <- P(sim)$reps
      gg_frp <- scfmutils::plot_fireRegimePolys(sim$fireRegimePolys)

      summaryDT <- rbindlist(lapply(allReps, function(r) {
        message("Loading saved simulation rep ", r, "/", max(allReps), " ...")
        fsim <- file.path(outputPath(sim), sprintf("rep%02d", r),
                          sprintf("%s_%04d.qs", P(sim)$simOutPrefix, P(sim)$simTimes[2]))
        if (!file.exists(fsim)) {
          fsim <- paste0(tools::file_path_sans_ext(fsim), ".rds") ## fallback to rds if qs not used
        }
        tmpSimPaths <- paths(sim)
        tmpSimPaths$outputPath <- dirname(fsim)

        tmp <- suppressMessages({
          loadSimList(fsim, paths = tmpSimPaths)
        })

        dt <- diagnosticPlotsDT(tmp)
        dt[, rep := r]

        message("  done")

        return(dt)
      }))

      write.csv(summaryDT, file.path(outputPath(sim), "scfmDiagnostics_multi_summary_dt.csv"))

      gg_fri <- scfmutils::comparePredictions_fireReturnInterval(
        summaryDT, list(start = P(sim)$simTimes[1], end = P(sim)$simTimes[2])) +
        geom_smooth(method = lm)

      gg_ign <- scfmutils::comparePredictions_annualIgnitions(summaryDT) +
        geom_smooth(method = lm)

      gg_mfs <- scfmutils::comparePredictions_meanFireSize(summaryDT) +
        geom_smooth(method = lm)

      gg_esc <- scfmutils::comparePredictions_annualEscapes(summaryDT) +
        geom_smooth(method = lm)

      ## note historical distribution is derived purely from historical data
      gg_histDist <- comparePredictions_fireDistribution(sim$fireRegimePoints,
                                                         size = min(sim$fireRegimePolys$cellSize),
                                                         burnSummary = sim$burnSummary)
      ## note that fireRegimePoints may include SAL but this figure only compares distribution
      ## so total area is irrelevant

      ## removed MAAB as diagnostic plot because it was derived from fire points incorrectly when
      ## studyAreaCalibration is supplied; MAAB can still be calculated manually if needed by user ## TODO

      if ("png" %in% P(sim)$.plots) {
        ggsave(file.path(figurePath(sim), "multi_FRI.png"), gg_fri, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "multi_FRP.png"), gg_frp, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "multi_IGN.png"), gg_ign, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "multi_MFS.png"), gg_mfs, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "multi_ESC.png"), gg_esc, height = 8, width = 8)
        ggsave(file.path(figurePath(sim), "multi_histDist.png"), gg_histDist, height = 8, width = 8)
      }

      if ("screen" %in% P(sim)$.plots) {
        clearPlot()
        gridExtra::grid.arrange(gg_frp,  gg_mfs,  gg_fri,  gg_ign, gg_esc,  gg_histDist,
                                nrow = 2, ncol = 3)
      }

      sim$scfmSummaryDT <- summaryDT
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

diagnosticPlotsDT <- function(sim) {
  sAR <- sim$studyAreaReporting |>
    sf::st_as_sf() |> #in case it is terra
    sf::st_union() |>
    sf::st_make_valid() |>
    terra::vect()

  fireRegimePointsReporting <- postProcess(sim$fireRegimePoints, to = sAR)
  #avoid geometry objects
  frpr <- postProcess(terra::vect(sim$fireRegimePolys),
                                          to = sAR)
  #drop true slivers, not ecological slivers
  frpr <- frpr[expanse(frpr) > res(sim$flammableMap)[1],]
  fireRegimePolysReporting <- sf::st_as_sf(frpr)

  # fireRegimePolysReporting <- sf::st_cast(fireRegimePolysReporting, "MULTIPOLYGON")

  #do not collection extract - it breaks.

  colsToDrop <- c("burnyArea", "nFlammable", "cellSize", paste0("nNbr_", 0:8))
  colsToKeep <- setdiff(names(fireRegimePolysReporting), colsToDrop)
  fireRegimePolysReporting <- fireRegimePolysReporting[colsToKeep]

  fireRegimePolysReporting <- genFireMapAttr(
    flammableMap = postProcessTo(sim$flammableMap, to = sAR),
    fireRegimePolys = fireRegimePolysReporting,
    neighbours = 8 ## TODO: use the param from the sim rather than hardcoding here
  )

  polyNames <- as.character(unique(fireRegimePolysReporting$PolyID))

  dt <- scfmutils::comparePredictions_summaryDT(
    fireRegimePoints = fireRegimePointsReporting,
    fireRegimePolys = fireRegimePolysReporting,
    burnSummary = sim$burnSummary, ## already summarized for studyAreaReporting
    times = times(sim)
  )

  return(dt)
}

## older version of SpaDES.core used here doesn't have this function
if (packageVersion("SpaDES.core") < "2.0.2.9001") {
  figurePath <- function(sim) {
    file.path(outputPath(sim), "figures", current(sim)[["moduleName"]]) |>
      checkPath(create = TRUE)
  }
}
