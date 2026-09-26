## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of a septuaginta
#' @description An object of class \code{\link[rgl]{mesh3d}} from package
#'    \strong{rgl} representing a septuaginta, i.e. a polyhedron with 72
#'    faces like a sphere based on the drawing by Leonardo da Vinci.
#' @format A list with components \code{vb}, \code{it}, and \code{normals}.
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "Septuaginta.R",       package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(-10, -20, zoom=0.8)
#' shade3d(dataSeptuaginta, color="orangered")
"dataSeptuaginta"
