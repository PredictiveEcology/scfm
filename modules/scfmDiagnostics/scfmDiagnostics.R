defineModule(sim, list(
  name = "scfmDiagnostics",
  description = paste("Diagnostic plots for scfm. Can be run for a single simulation,",
                      "as part of the main scfm run, or as part of postprocessing.",
                      "Inputs objects will be loaded from saved simulation files when in 'multi' mode."),
  keywords = c("diagnostics", "scfm"),
  authors = c(
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(scfmDiagnostics = "0.0.4"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmDiagnostics.Rmd"), ## same file
  reqdPkgs = list("ggplot2", "gridExtra",
                  "PredictiveEcology/scfmutils (>= 0.0.7)",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.0.9001)"),
  parameters = bindrows(
    defineParameter("mode", "character", "single", NA, NA,
                    paste("use 'single' to run part of an scfm simulation (i.e., along with other scfm modules);",
                          "use 'multi' to run as part of postprocessing multiple scfm runs.")),
    defineParameter("reps", "integer", NA_integer_, 1L, NA_integer_,
                    paste("number of replicates/runs per study area when running in 'multi' mode.")),
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
                 "describes details of all burned pixels. Required in single mode.", sourceURL = NA),
    expectsInput("fireRegimePoints", "SpatialPointsDataFrame",
                 "Fire locations. Points outside studyArea are removed. Required in single mode.", sourceURL = NA),
    expectsInput("fireRegimePolys", "sf",
                 paste("Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
                       "Must have numeric field 'PolyID' or it will be created for individual polygons.",
                       "Required in single mode."),
                 sourceURL = NA),
    expectsInput("flammableMap", "RasterLayer",
                 desc = "binary flammability map. Required in single mode.", sourceURL = NA),
    expectsInput("landscapeAttr", "list",
                 "contains landscape attributes for each polygon. Required in single mode.", sourceURL = NA),
    expectsInput("scfmDriverPars", "list",
                 "burn parameters for each polygon in `fireRegimePolys`.  Required in single mode.", sourceURL = NA),
    expectsInput("scfmRegimePars", "list",
                 "list of fire regime parameters for each polygon. Required in single mode.", sourceURL = NA),
    expectsInput("studyAreaReporting", "sf",
                 paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                       "Required in single mode."))
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
      # ! ----- EDIT BELOW ----- ! #

      dt <- diagnosticPlotsDT(sim)

      ## Some useful plots
      gg_fri <- scfmutils::comparePredictions_fireReturnInterval(dt, times = times(sim))
      gg_frp <- scfmutils::plot_fireRegimePolys(sim$fireRegimePolys)
      gg_ign <- scfmutils::comparePredictions_annualIgnitions(dt)
      gg_mfs <- scfmutils::comparePredictions_meanFireSize(dt)

      # removed MAAB as diagnostic plot because it was derived from fire points incorrectly when SAL is supplied
      # MAAB can still be calculated manually if a user desires ## TODO

      if ("png" %in% P(sim)$.plots) {
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_FRI.png"), gg_fri, height = 8, width = 8)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_FRP.png"), gg_frp, height = 8, width = 8)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_IGN.png"), gg_ign, height = 8, width = 8)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_MFS.png"), gg_mfs, height = 8, width = 8)
      }

      if ("screen" %in% P(sim)$.plots) {
        clearPlot()
        gridExtra::grid.arrange(gg_frp,  gg_mfs,  gg_fri,  gg_ign, nrow = 2, ncol = 2)
      }

      sim$scfmSummaryDT <- dt

      # ! ----- STOP EDITING ----- ! #
    },
    diagnosticPlotsMulti = {
      # ! ----- EDIT BELOW ----- ! #

      allReps <- P(sim)$reps
      gg_frp <- NULL

      summaryDT <- rbindlist(lapply(allReps, function(r) {
        fsim <- file.path(outputPath(sim), sprintf("rep%02d", r), sprintf("mySimOut_%04d.qs", P(sim)$simTimes[2]))
        tmp <- suppressMessages(loadSimList(fsim))
        if (r == allReps[1]) {
          gg_frp <<- scfmutils::plot_fireRegimePolys(tmp$fireRegimePolys)
        }

        dt <- diagnosticPlotsDT(tmp)
        dt[, rep := r]
      }))

      gg_fri <- scfmutils::comparePredictions_fireReturnInterval(
        summaryDT, list(start = P(sim)$simTimes[1], end = P(sim)$simTimes[2])) +
        geom_smooth(method = lm)

      gg_ign <- scfmutils::comparePredictions_annualIgnitions(summaryDT) +
        geom_smooth(method = lm)

      gg_mfs <- scfmutils::comparePredictions_meanFireSize(summaryDT) +
        geom_smooth(method = lm)

      # removed MAAB as diagnostic plot because it was derived from fire points incorrectly when SAL is supplied
      # MAAB can still be calculated manually if a user desires ## TODO

      if ("png" %in% P(sim)$.plots) {
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_multi_FRI.png"), gg_fri, height = 8, width = 8)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_multi_FRP.png"), gg_frp, height = 8, width = 8)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_multi_IGN.png"), gg_ign, height = 8, width = 8)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnostics_multi_MFS.png"), gg_mfs, height = 8, width = 8)
      }

      if ("screen" %in% P(sim)$.plots) {
        clearPlot()
        gridExtra::grid.arrange(gg_frp,  gg_mfs,  gg_fri,  gg_ign, nrow = 2, ncol = 2)
      }

      sim$scfmSummaryDT <- summaryDT

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

diagnosticPlotsDT <- function(sim) {
  fireRegimePointsReporting <- sf::st_intersection(sim$fireRegimePoints, sim$studyAreaReporting)

  fireRegimePolysReporting <- sf::st_intersection(sim$fireRegimePolys, sim$studyAreaReporting)

  landscapeAttrReporting <- genFireMapAttr(
    flammableMap = postProcess(sim$flammableMap, studyArea = sim$studyAreaReporting),
    fireRegimePolys = fireRegimePolysReporting,
    neighbours = 8 ## TODO: use the param from the sim rather than hardcoding here
  )

  polyNames <- as.character(unique(fireRegimePolysReporting$PolyID))

  stopifnot(all(polyNames %in% names(landscapeAttrReporting)))

  scfmDriverParsReporting <- subset(sim$scfmDriverPars, names(sim$scfmDriverPars) %in% polyNames)

  scfmRegimeParsReporting <- subset(sim$scfmRegimePars, names(sim$scfmRegimePars) %in% polyNames)

  dt <- scfmutils::comparePredictions_summaryDT(
    scfmDriverPars = scfmDriverParsReporting,
    scfmRegimePars = scfmRegimeParsReporting,
    landscapeAttr = landscapeAttrReporting,
    fireRegimePoints = fireRegimePointsReporting,
    burnSummary = sim$burnSummary, ## already summarized for studyAreaReporting
    times = times(sim)
  )

  return(dt)
}
