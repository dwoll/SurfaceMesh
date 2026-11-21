# static PNG snapshots need a real OpenGL device. Check if one is available
# (not on, e.g., CRAN macOS machines: "no conforming visual").
options(rgl.useNULL=FALSE)
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(SurfaceMesh))
options(rgl.printRglwidget=FALSE)

hasGL <- tryCatch({
    !rgl.useNULL() && suppressWarnings(open3d()) > 0 && rgl.cur() > 0
}, error=function(e) FALSE)

if(!hasGL) {
    try(while(rgl.cur() > 0) close3d(), silent=TRUE)
    options(rgl.useNULL=TRUE)
}

# no display: skip chunks that call rgl plotting functions (eval=FALSE)
# keyword 'rgl3d' is a chunk option (TRUE by default) that triggers the hook
rglPattern <- "open3d|view3d|shade3d|wire3d|points3d|lines3d|segments3d|spheres3d|text3d|plot3d|next3d|mfrow3d|bg3d|rglwidget|snapshot3d|par3d|axes3d|title3d|decorate3d|bbox3d|aspect3d|triangles3d|quads3d|planes3d)"
knitr::opts_chunk$set(rgl3d=TRUE)
knitr::opts_hooks$set(rgl3d=function(options) {
    if(!hasGL && any(grepl(rglPattern, options$code))) { options$eval <- FALSE }
    options
})

if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !rmarkdown::pandoc_available("1.14")) {
  warning(call. = FALSE, "These vignettes assume rmarkdown and Pandoc
          version 1.14. These were not found. Older versions will not work.")
  knitr::knit_exit()
}

# If Pandoc is not installed, the output format won't be set.
# knitr uses it to determine whether to do
# screenshots; we don't want those.
# see https://github.com/rstudio/markdown/issues/115
knitr::opts_chunk$set(screenshot.force = FALSE, snapshot = TRUE)
# snapshot = TRUE for snapshots (PNG) instead of dynamic widgets
