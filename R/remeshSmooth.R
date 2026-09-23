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

#' @title Smooth shape
#' @description Smooths the overall shape of a 3D surface mesh by using
#'   the mean curvature flow.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param indices Indices of the faces to be smoothed.
#'   If missing, the whole mesh is smoothed.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param time Positive number: A time step that corresponds to the speed by
#'   which the surface is smoothed (the larger the faster).
#'   typical values lie between \code{1e-6} and \code{1}. See details.
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh       <- dataHeart1
#' mesh_rgl   <- toRGL(mesh)
#' mesh_s     <- remeshSmoothShape(mesh, nIter=5, time=1)
#' mesh_s_rgl <- toRGL(mesh_s)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_s_rgl)
#'
#' @export
remeshSmoothShape <- function(x, indices, nIter = 1L, time = 0.001, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
           " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isPositiveNumber(time))
  stopifnot(isBoolean(normals))
  if(missing(indices)) {
    indices <- integer(0L)
  } else {
    stopifnot(isAtomicVector(indices))
    stopifnot(is.numeric(indices))
    integers <- isTRUE(all.equal(indices, floor(indices)))
    if(!integers) {
      stop("The indices must be positive integers.")
    }
    positive <- all(indices >= 1L)
    if(!positive) {
      stop("The indices must be positive integers.")
    }
    indices <- unique(as.integer(indices)) - 1L
  }
  meshCPP <- fromR(x)
  meshOut <- remeshSmoothShape_cpp(
    meshCPP, indices, as.integer(nIter), time, normals)
  fromCPP(meshOut)
}

#' @title Smooth angle
#' @description Smoothing by angle optimization.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param dihedralAngle Positive number. Constrain edges with a dihedral
#'   angle over given value.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param useSafeConstr Boolean. Use safety constraints for moving vertices?
#' @param doProject Boolean. Project points onto initial surface after
#'   each iteration?
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/} for details.
#'   Note that CGAL `angle_and_area_smoothing()` supports smoothing by area
#'   optimization as well. However, this feature requires the external
#'   Ceres library (\url{http://ceres-solver.org/}) and is thus currently
#'   not enabled here.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh       <- dataHeart1
#' mesh_rgl   <- toRGL(mesh)
#' mesh_s     <- remeshSmoothAngle(mesh, nIter=5)
#' mesh_s_rgl <- toRGL(mesh_s)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_s_rgl)
#'
#' @export
remeshSmoothAngle <- function(x,
                              dihedralAngle = 60,
                              nIter = 1L,
                              useSafeConstr = FALSE,
                              doProject = TRUE,
                              normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
           " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isPositiveNumber(dihedralAngle))
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(useSafeConstr))
  stopifnot(isBoolean(normals))
  storage.model(dihedralAngle) <- "double"
  meshCPP <- fromR(x)
  meshOut <- remeshSmoothAA_cpp(meshCPP,
                                dihedralAngle,
                                as.integer(nIter),
                                useSafeConstr,
                                doProject,
                                normals)
  fromCPP(meshOut)
}

#' @title Smooting by tangential relaxation
#' @description Iterative area-based smoothing by tangential relaxation.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param relaxConstr Boolean. If \code{TRUE}, the end vertices of the edges
#'   set as constrained and boundary edges move along the constrained
#'   polylines they belong to.
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/} for details.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh       <- dataHeart1
#' mesh_rgl   <- toRGL(mesh)
#' mesh_s     <- remeshSmoothTangentRelax(mesh, nIter=5)
#' mesh_s_rgl <- toRGL(mesh_s)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_s_rgl)
#'
#' @export
remeshSmoothTangentRelax <- function(x,
                                     nIter = 1L,
                                     relaxConstr = FALSE,
                                     normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
           " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(relaxConstr))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  meshOut <- remeshSmoothTR_cpp(
    meshCPP, as.integer(nIter), relaxConstr, normals)
  fromCPP(meshOut)
}
