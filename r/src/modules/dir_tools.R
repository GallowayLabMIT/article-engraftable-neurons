#' Sets the current working directory to the data directory. When filename is not specified,
#' uses a git-ignored datadir_path.txt file relative to the root Rproj directory to set the
#' working directory. Returns the original working directory
#'
#' @param data_path Optionally, a path to the data directory
#' @return The old working directory
set_datadir <- function(data_path = NA) {
  old_wd <- getwd()
  if (is.na(data_path)) {
    setwd(readLines('datadir.txt', n = 1))
  } else {
    setwd(data_path)
  }
  return(old_wd)
}
