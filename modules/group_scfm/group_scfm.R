defineModule(sim, list(
  name = "group_scfm",
  description = "parent module for scfm family of modules",
  keywords = "fire",
  authors = c(
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = c("ageModule", "scfmDiagnostics", "scfmDriver", "scfmEscape", "scfmIgnition",
                   "scfmLandcoverInit", "scfmRegime", "scfmSpread"),
  version = list(group_scfm = "0.0.1"),

  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "group_scfm.Rmd") ## README rendered from Rmd
))
