## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of an Entzensberger star
#' @description An object of class \code{\link[rgl]{mesh3d}} from package
#'    \strong{rgl} representing an Entzensberger star.
#' @format A list with components \code{vb}, \code{it}, and \code{normals}.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "EntzenStar.R",        package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(-55, -20, zoom=0.8)
#' shade3d(dataEntzenStar, color="steelblue2")
"dataEntzenStar"
