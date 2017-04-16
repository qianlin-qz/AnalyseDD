#' shinny
#' @export
runExemple <- function() {
  appDir <- system.file("shiny-examples","afc", package = "AnalyseDD")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `AnalyseDD`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

