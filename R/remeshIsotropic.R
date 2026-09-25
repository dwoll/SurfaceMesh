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

#' @title Isotropic remeshing
#' @description Isotropic remeshing of a triangular 3D surface mesh.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param method \code{character}. Either \code{"uniform"} for uniform sizing field
#'   (all triangle edges are targeted to have equal lengths)
#'   or \code{"adaptive"} for adaptive sizing field
#'   (triangle edge lengths depend on the local curvature).
#' @param targetEdgeLen Positive number for \code{method="uniform"}:
#'   The target edge length of the remeshed mesh.
#' @param tol Positive number for \code{method="adaptive"}.
#'   Error tolerance. See details.
#' @param edgeMin Positive number for \code{method="adaptive"}.
#'   Minimum edge length. See details.
#' @param edgeMax Positive number for \code{method="adaptive"}.
#'   Maximum edge length. See details.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param nRelaxSteps Positive \code{integer}: Number of relaxation steps.
#' @param protectConstraints Boolean. Protect constraints? See details.
#' @param dihedralAngle Positive number. Constrain edges with a dihedral
#'   angle over given value to preserve sharp edges.
#' @param normals Boolean. Return vertex normals?
#' @return A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/} for details.
#' @seealso See \code{\link[SurfaceMesh]{remeshSimplify}},
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
#' mesh         <- makeMesh(dataTruncIcosahedron, triangulate=TRUE)
#' mesh_rgl     <- toRGL(mesh)
#' mesh_rem     <- remeshIsotropic(mesh, targetEdgeLen=0.7)
#' mesh_rem_rgl <- toRGL(mesh_rem)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_rem_rgl)
#'
#' @export
remeshIsotropic <- function(
        x,
        method = c("uniform", "adaptive"),
        tol = 0.001,
        edgeMin = 0.001,
        edgeMax,
        targetEdgeLen,
        nIter = 1L,
        nRelaxSteps = 1L,
        protectConstraints = TRUE,
        dihedralAngle = 60,
        normals = FALSE) {
    method <- match.arg(method)
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    stopifnot(isStrictPositiveInteger(nIter))
    stopifnot(isStrictPositiveInteger(nRelaxSteps))
    stopifnot(isBoolean(protectConstraints))
    stopifnot(isNonNegativeNumber(dihedralAngle))
    stopifnot(isBoolean(normals))
    meshCPP <- fromR(x)
    meshRem <- if(method == "uniform") {
      stopifnot(isPositiveNumber(targetEdgeLen))
      remeshIsoUniform_cpp(meshCPP,
                           targetEdgeLen,
                           as.integer(nIter),
                           as.integer(nRelaxSteps),
                           protectConstraints,
                           dihedralAngle,
                           normals)
    } else if(method == "adaptive") {
      stopifnot(isPositiveNumber(tol))
      stopifnot(isPositiveNumber(edgeMin))
      stopifnot(isPositiveNumber(edgeMax))
      remeshIsoAdapt_cpp(meshCPP,
                         tol,
                         edgeMin,
                         edgeMax,
                         as.integer(nIter),
                         as.integer(nRelaxSteps),
                         protectConstraints,
                         dihedralAngle,
                         normals)
    }

    fromCPP(meshRem)
}
