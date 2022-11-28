defineModule(sim, list(
  name = "scfmDiagnostics",
  description = "diagnostic plots for scfm",
  keywords = c("diagnostics", "scfm"),
  authors = c(
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(scfmDiagnostics = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfmDiagnostics.Rmd"), ## same file
  reqdPkgs = list("ggplot2", "gridExtra",
                  "PredictiveEcology/scfmutils (>= 0.0.1)",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.0.9001)"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
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
                 "describes details of all burned pixels", sourceURL = NA),
    expectsInput("fireRegimePoints", "SpatialPointsDataFrame",
                 "Fire locations. Points outside studyArea are removed", sourceURL = NA),
    expectsInput("fireRegimePolys", "sf",
                 paste("Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada.",
                       "Must have numeric field 'PolyID' or it will be created for individual polygons."),
                 sourceURL = NA),
    expectsInput("landscapeAttr", "list",
                 "contains landscape attributes for each polygon.", sourceURL = NA),
    expectsInput("scfmDriverPars", "list",
                 "burn parameters for each polygon in `fireRegimePolys`", sourceURL = NA),
    expectsInput("scfmRegimePars", "list",
                 "list of fire regime parameters for each polygon.", sourceURL = NA)
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
      sim <- scheduleEvent(sim, end(sim), "scfmDiagnostics", "diagnosticPlots", eventPriority = .last())
    },
    diagnosticPlots = {
      # ! ----- EDIT BELOW ----- ! #

      dt <- scfmutils::comparePredictions_summaryDT(scfmDriverPars = sim$scfmDriverPars,
                                                    scfmRegimePars = sim$scfmRegimePars,
                                                    landscapeAttr = sim$landscapeAttr,
                                                    fireRegimePoints = sim$fireRegimePoints,
                                                    burnSummary = sim$burnSummary,
                                                    times = times(sim))

      ## Some useful plots
      gg_fri <- scfmutils::comparePredictions_fireReturnInterval(dt, times = times(sim))
      gg_frp <- scfmutils::plot_fireRegimePolys(sim$fireRegimePolys)
      gg_ign <- scfmutils::comparePredictions_annualIgnitions(dt)
      gg_mfs <- scfmutils::comparePredictions_meanFireSize(dt)

      # removed MAAB as diagnostic plot because it was derived from fire points incorrectly when SAL is supplied
      # MAAB can still be calculated manually if a user desires ## TODO

      if ("png" %in% P(sim)$.plots) {
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnstics_FRI.png"), gg_fri)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnstics_FRP.png"), gg_frp)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnstics_IGN.png"), gg_ign)
        ggsave(file.path(outputPath(sim), "figures", "scfmDiagnstics_MFS.png"), gg_mfs)
      }

      if ("screen" %in% P(sim)$.plots) {
        clearPlot()
        gridExtra::grid.arrange(gg_frp,  gg_mfs,  gg_fri,  gg_ign, nrow = 2, ncol = 2)
      }

      sim$scfmSummaryDT <- dt

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}
