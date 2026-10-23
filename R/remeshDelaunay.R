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

#' @title Surface Delaunay remeshing
#' @description Remeshing algorithm based on Delaunay refinement of a
#'   restricted Delaunay triangulation.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param edgeSize Positive number: Upper bound for the lengths of curve
#'   edges. Has to be set to a positive value when \code{protectConstraints=TRUE}.
#' @param facetAngle Positive number: Lower bound for the angles (in degrees)
#'   of the surface mesh facets.
#' @param facetDist Positive number. Upper bound for the distance between the
#'   facet circumcenter and the center of its surface Delaunay ball.
#' @param protectConstraints Boolean. Protect sharp edges.
#'   See \code{dihedralAngle}.
#' @param dihedralAngle Positive number. Constrain edges with a dihedral
#'   angle over given value to protect sharp edges.
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/} for details.
#' @seealso See \code{\link[SurfaceMesh]{remeshACVD}},
#'   \code{\link[SurfaceMesh]{remeshIsotropic}},
#'   \code{\link[SurfaceMesh]{remeshSimplify}},
#'   \code{\link[SurfaceMesh]{remeshSmoothShape}},
#'   \code{\link[SurfaceMesh]{remeshSmoothAngle}},
#'   \code{\link[SurfaceMesh]{remeshSmoothTangentRelax}}
#'   for other remeshing operations.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh       <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_rgl   <- toRGL(mesh)
#' mesh_d     <- remeshDelaunay(mesh,
#'                              edgeSize=0.2,
#'                              facetDist=0.2,
#'                              protectConstraints=TRUE)
#' mesh_d_rgl <- toRGL(mesh_d)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 30, zoom=0.8)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 30, zoom=0.8)
#' wire3d(mesh_d_rgl)
#'
#' @export
remeshDelaunay <- function(x,
                           edgeSize,
                           facetAngle = 0.0,
                           facetDist = 0.5,
                           protectConstraints = TRUE,
                           dihedralAngle = 60.0,
                           normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
           " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(protectConstraints))
  if(protectConstraints) {
    stopifnot(isPositiveNumber(edgeSize))
  } else {
    edgeSize = -1.0
  }
  stopifnot(isNonNegativeNumber(facetAngle))
  stopifnot(isPositiveNumber(facetDist))
  stopifnot(isNonNegativeNumber(dihedralAngle))
  stopifnot(isBoolean(normals))
  storage.mode(edgeSize)       <- "double"
  storage.mode(facetAngle)     <- "double"
  storage.mode(facetDist)      <- "double"
  storage.mode(dihedralAngle)  <- "double"
  meshCPP <- fromR(x)
  meshOut <- remeshDelaunay_cpp(
    meshCPP, edgeSize, facetAngle, facetDist, protectConstraints, dihedralAngle, normals)
  fromCPP(meshOut)
}
