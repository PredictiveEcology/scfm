## The scfm family as one module: this repository is the parent, and its children are the
## modules in `modules/`. `PredictiveEcology/scfm@v2.1.0` loads the parent and these children
## from the same release (SpaDES.project >= 1.2.0.9015). Every module in the repository
## carries the family version.
defineModule(sim, list(
  name = "scfm",
  description = "Parent module for the scfm family of fire modules (ignition, escape, spread)",
  keywords = "fire",
  authors = c(
    person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = c("scfmDataPrep", "scfmIgnition", "scfmEscape", "scfmSpread",
                   "scfmDiagnostics"),
  version = list(scfm = "2.1.0.9000",
                 scfmDataPrep = "2.1.0", scfmIgnition = "2.1.0", scfmEscape = "2.1.0",
                 scfmSpread = "2.1.0", scfmDiagnostics = "2.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfm.Rmd")
))
