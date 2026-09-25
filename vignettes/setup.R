options(rgl.useNULL=FALSE)
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(SurfaceMesh))
options(rgl.useNULL=TRUE)
options(rgl.printRglwidget=FALSE)
open3d()

if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !rmarkdown::pandoc_available("1.14")) {
  warning(call. = FALSE, "These vignettes assume rmarkdown and Pandoc
          version 1.14. These were not found. Older versions will not work.")
  knitr::knit_exit()
}

# If Pandoc is not installed, the output format won't be set.
# knitr uses it to determine whether to do
# screenshots; we don't want those. see https://github.com/rstudio/markdown/issues/115
knitr::opts_chunk$set(screenshot.force = FALSE, snapshot = TRUE)
# snapshot = TRUE for snapshots instead of dynamic
