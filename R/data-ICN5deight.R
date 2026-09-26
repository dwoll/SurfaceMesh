## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of an ICN5d eight
#' @description An object of class \code{\link[rgl]{mesh3d}} from package
#'   \strong{rgl} representing an ICN5D eight, a shape from from user ICN5D
#'   on forum \url{http://hi.gher.space/forum/}.
#' @format A list with components \code{vb}, \code{it}, and \code{normals}.
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "ICN5Deight.R",        package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(35, -20, zoom=0.8)
#' shade3d(dataICN5Deight, color="violetred")
"dataICN5Deight"
