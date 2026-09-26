## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of an oloid
#' @description A \code{CGALmesh} object, i.e., the output of
#'   \code{\link[SurfaceMesh]{makeMesh}} representing an oloid.
#' @format A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "Oloid.R",             package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(20, 15, zoom=0.8)
#' shade3d(toRGL(dataOloid), color="darkviolet")
"dataOloid"
