options(rgl.useNULL=FALSE)
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(SurfaceMesh))
options(rgl.useNULL=TRUE)
options(rgl.printRglwidget=FALSE)
open3d()

if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !rmarkdown::pandoc_available("1.14")) {
  warning(call. = FALSE, "These vignettes assume rmarkdown and Pandoc
          version 1.14.  These were not found. Older versions will not work.")
  knitr::knit_exit()
}

# If Pandoc is not installed, the output format won't be set.
# knitr uses it to determine whether to do
# screenshots; we don't want those. see https://github.com/rstudio/markdown/issues/115
knitr::opts_chunk$set(screenshot.force = FALSE, snapshot = TRUE)
# snapshot = TRUE for snapshots instead of dynamic

backticked <- function(s) { paste0("`", s, "`") }

# Write this once at the start of the document.
cat('<style>
    .nostripes tr.even {background-color: white;}
    table {border-style: none;}
    table th {border-style: none;}
    table td {border-style: none;}
    a[href^=".."] {text-decoration: underline;}
    </style>
    ')

# This displays the string code as `r code` when entered
# as `r rinline(code)`.  Due to Stephane Laurent
rinline <- function(code, script = FALSE){
  if (script)
    html <- "`r CODE`"
  else
    html <- '<code  class="r">``` `r CODE` ```</code>'
  sub("CODE", code, html)
}

# This sets up default "alt text" for screen readers.
defaultAltText <- function() {
  paste(knitr::opts_current$get("label"), "example.")
}

knitr::opts_chunk$set(fig.alt = quote(defaultAltText()))
