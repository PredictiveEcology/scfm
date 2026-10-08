defineModule(sim, list(
  name = "group_scfm",
  description = "Deprecated: use the `scfm` parent module at the root of this repository (PredictiveEcology/scfm@v2.1.0)",
  keywords = "fire",
  authors = c(
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = c("ageModule", "scfmDataPrep", "scfmDiagnostics",
                   "scfmEscape", "scfmIgnition", "scfmSpread"),
  version = list(group_scfm = "2.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "group_scfm.Rmd") ## README rendered from Rmd
))
