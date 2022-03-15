defineModule(sim, list(
  name = "group_scfm",
  description = "parent module for scfm family of modules",
  keywords = "fire",
  authors = c(
    person(c("Steve", "Cumming"), "Last", email = "email@example.com", role = c("aut", "cre"))
  ),
  childModules = c("ageModule", "scfmDriver", "scfmEscape", "scfmIgnition", "scfmLandcoverInit",
                   "scfmRegime", "scfmSpread"),
  version = list(SpaDES.core = "0.2.3", group_scfm = "0.0.1"),

  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "group_scfm.Rmd")
))
