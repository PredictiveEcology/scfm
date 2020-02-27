#' Get fire points data from the Canadian National Fire Database
#'
#' @param url URL from which to download the fire points data. Default \code{NULL} fetches data from
#'            \url{"http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"}.
#' @studyArea description needed
#' @param redownloadIn time in years that we tolerate the data to be "old", and require redownload.
#'                     I.e. 0.5 would mean "redownload data older than 6 months". Default 1.
#' @param NFDB_pointPath file path to save the download data. Must be provided.
getFirePoints_NFDB <- function(url = NULL,
                               studyArea = NULL, rasterToMatch = NULL, redownloadIn = 1,
                               NFDB_pointPath = NULL) {
  if (is.null(NFDB_pointPath)) stop("NFDB_pointPath cannot be null. Specify a file path.")

  if (is.null(url))
    url <- "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"

  check <- Checksums(NFDB_pointPath, checksumFile = file.path(NFDB_pointPath, "CHECKSUMS.txt"), write = TRUE)
  whRowIsShp <- grep("NFDB_point.*shp$", check$expectedFile)
  whIsOK <- which(check$result[whRowIsShp] == "OK")
  needNewDownload <- TRUE
  if (any(whIsOK)) {
    filesToCheck <- tools::file_path_sans_ext(unlist(lapply(check[whRowIsShp[whIsOK], "expectedFile"], as.character)))
    dateOfFile <- substr(x = filesToCheck, start = nchar(filesToCheck) - 8 +
             1, nchar(filesToCheck))
    if (any((as.Date(dateOfFile, format = "%Y%m%d") + dyear(redownloadIn)) > Sys.Date())) {
      needNewDownload <- FALSE
    }
  }
  if (needNewDownload) {
    print("downloading NFDB")
    firePoints <- Cache(prepInputs, url = url, studyArea = studyArea, fun = "shapefile",
                        destination = NFDB_pointPath, useCache = "overwrite",
                        useSAcrs = TRUE, omitArgs = c("NFDB_pointPath", "overwrite"))
  } else {
    NFDBs <- grep(list.files(NFDB_pointPath), pattern = "^NFDB", value = TRUE)
    shps <- grep(list.files(NFDB_pointPath), pattern = ".shp$", value = TRUE)
    aFile <- NFDBs[NFDBs %in% shps][1] #in case there are multiple files
    firePoints <- Cache(shapefile, file.path(NFDB_pointPath, aFile))
    firePoints <- Cache(postProcess, x = firePoints,
                            studyArea = studyArea, filename2 = NULL,
                            rasterToMatch = rasterToMatch,
                            userTags = c("cacheTags", "NFDB"))
  }
  return(firePoints)
}
