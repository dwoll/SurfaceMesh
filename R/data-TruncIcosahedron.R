## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of the truncated icosahedron
#'
#' @description A list giving the vertices and the faces of a truncated
#'   icosahedron. There are some hexagonal faces and some pentagonal faces.
#' @format A list with two components: \code{vertices} and \code{faces}.
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "TruncIcosahedron.R",  package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' mTruncIcosahedron <- makeMesh(dataTruncIcosahedron, triangulate=TRUE)
#' view3d(-15, -10, zoom=0.8)
#' shade3d(toRGL(mTruncIcosahedron), color="wheat2")
"dataTruncIcosahedron"
