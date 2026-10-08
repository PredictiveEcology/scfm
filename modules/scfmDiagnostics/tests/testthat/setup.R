## Template `tests/testthat/setup.R` for a SpaDES module.
##
## The module is converted to a package rendition before the tests run (see the
## testthat-module workflow), so the module's own functions are already in a real
## namespace: a test can call `myHelper(x)` directly rather than reaching into
## `sim@.xData$.mods$<module>$myHelper`.
##
## This file only sets options and makes a scratch directory. It deliberately does
## NOT install packages or download anything: dependency resolution happens once,
## before the tests, from the module's own `reqdPkgs` metadata.

withr::local_options(
  list(
    reproducible.useMemoise = TRUE,
    reproducible.verbose    = -2,
    Require.verbose         = -2,
    spades.moduleCodeChecks = FALSE,
    spades.moduleDocument   = FALSE,
    spades.useRequire       = FALSE
  ),
  .local_envir = testthat::teardown_env()
)

## Where the module is, from testthat's working directory.
##
## testthat runs with the working directory set to `tests/testthat`, so the module
## directory is two levels up and the `modulePath` that SpaDES.core wants -- the
## directory *containing* modules -- is three. Getting this wrong is silent: the
## module is simply not found, and every test that needs it skips or errors on
## something unrelated.
moduleRoot <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
moduleName <- basename(moduleRoot)
modulePath <- dirname(moduleRoot)

## A scratch tree removed when the suite finishes. Write here, never beside the
## module: the tests run against a throwaway copy, but the habit matters when they
## are run by hand.
testPaths <- local({
  root <- withr::local_tempdir(.local_envir = testthat::teardown_env())
  paths <- list(
    cachePath  = file.path(root, "cache"),
    inputPath  = file.path(root, "inputs"),
    modulePath = modulePath,
    outputPath = file.path(root, "outputs")
  )
  for (p in paths[c("cachePath", "inputPath", "outputPath")]) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  paths
})
