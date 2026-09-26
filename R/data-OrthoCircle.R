## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of an Ortho Circle
#' @description An object of class \code{\link[rgl]{mesh3d}} from package
#'    \strong{rgl} representing an Ortho Circle.
#' @format A list with components \code{vb}, \code{it}, and \code{normals}.
#' @details \url{https://paulbourke.net/geometry/orthocircle/}
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "OrthoCircle.R",       package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(-15, -10, zoom=0.8)
#' shade3d(dataOrthoCircle, color="darkolivegreen4")
"dataOrthoCircle"
