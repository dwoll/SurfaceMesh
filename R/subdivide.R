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

#' @title Subdivision and deformation
#' @description Performs the Catmull-Clark, Doo-Sabin, Loop, and Sqrt3
#'   subdivision and deformation of a 3D surface mesh.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param method \code{character}. One of \code{"CC"} (Catmull-Clark),
#'   \code{"DS"} (Doo-Sabin), \code{"L"} (Loop), \code{"S3"} (Sqrt3).
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param triangulate Boolean. For \code{method="DS"}. Triangulate resulting mesh?
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/Subdivision_method_3/} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh        <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_cc     <- subdivideCatmullClark(mesh, nIter=2L)
#' mesh_cc_rgl <- toRGL(mesh_cc)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_cc_rgl)
#'
#' @export
subdivide <- function(x,
                      method = c("CC", "DS", "Loop", "Sqrt3"),
                      nIter = 1L,
                      triangulate = TRUE,
                      normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			     " (i.e., the output of the `makeMesh()` function).")
  }
  method <- match.arg(method)
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(triangulate))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  meshOut <- if(method == "CC") {
    subdivideCatmullClark_cpp(meshCPP, as.integer(nIter), normals)
  } else if(method == "DS") {
    subdivideDooSabin_cpp(meshCPP, as.integer(nIter), triangulate, normals)
  } else if(method == "Loop") {
    subdivideLoop_cpp(meshCPP, as.integer(nIter), normals)
  } else if(method == "Sqrt3") {
    subdivideSqrt3_cpp(meshCPP, as.integer(nIter), normals)
  } else {
    stop("Wrong method.")
  }
  fromCPP(meshOut)
}
