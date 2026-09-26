## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of a pentagrammic prism
#'
#' @description A list representing a pentagrammic prism, giving the vertices
#'   and the faces; it has 20 vertices, 10 triangular faces, 10 rectangular
#'   faces and two pentagonal faces.
#' @format A list with components \code{vertices}, \code{faces}.
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "PentaPrism.R",        package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' mPentaPrism <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' view3d(-15, -15, zoom=0.8)
#' shade3d(toRGL(mPentaPrism), color="forestgreen")
"dataPentaPrism"
