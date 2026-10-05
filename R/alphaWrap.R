## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title 3D alpha wrapping
#' @description Reconstruction of a 3D surface mesh from a cloud of 3D
#'   points by alpha wrapping.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of
#'   \code{\link[SurfaceMesh]{makeMesh}},
#'   a \code{\link[rgl]{mesh3d}} object from package \strong{rgl},
#'   or a numeric matrix with 3 columns which stores the point coordinates,
#'   one point per row, and at least 4 points.
#' @param alphaRel Relative alpha parameter. The actual alpha parameter (see
#'   details) is defined as the length of the diagonal of the bounding box of
#'   the point cloud divided by the relative alpha parameter. Increase for
#'   a more detailed mesh.
#' @param offsetRel Relative offset. The actual offset parameter (see details)
#'   is defined as the length of the diagonal of the bounding box of the
#'   point cloud divided by the relative offset parameter. Increase for
#'   output that is closer to input mesh.
#' @param normals Boolean. Return vertex normals?
#'
#' @returns A \code{CGALmesh} object.
#'
#' @details See \url{https://doc.cgal.org/latest/Alpha_wrap_3/} for details.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso See also \code{\link[SurfaceMesh]{reconstructAFS}},
#'    \code{\link[SurfaceMesh]{reconstructSSS}},
#'    \code{\link[SurfaceMesh]{reconstructPoisson}} for other surface reconstruction
#'    methods.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#' mesh        <- dataHeart1
#' mesh_rgl    <- toRGL(mesh)
#' mesh_aw     <- alphaWrap(mesh,
#'                          alphaRel =10,
#'                          offsetRel=300)
#' mesh_aw_rgl <- toRGL(mesh_aw)
#'
#' mfrow3d(1, 2)
#' view3d(0, -90, zoom=0.7)
#' shade3d(mesh_rgl, col="gray")
#' next3d()
#' view3d(0, -90, zoom=0.7)
#' shade3d(mesh_aw_rgl, col="gray")
#'
#' @export
#' @importFrom rgl tmesh3d
alphaWrap <- function(x, alphaRel, offsetRel, normals=FALSE) {
  stopifnot(isPositiveNumber(alphaRel))
  stopifnot(isPositiveNumber(offsetRel))
  stopifnot(isBoolean(normals))
  xIn     <- getVertsMat(x, nPtsMin=4L)
  meshOut <- alphaWrapPoints_cpp(t(xIn), alphaRel, offsetRel, normals)
  # meshCPP <- fromR(x)
  # alphaWrapMesh_cpp(meshCPP, alphaRel, offsetRel, normals)
  fromCPP(meshOut)
}
