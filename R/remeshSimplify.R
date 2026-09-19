## ----------------------------------------------------------------------- //
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title Triangulated Surface Mesh Simplification
#' @description Simplification of a triangular 3D surface mesh.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param repairSoup Boolean. Whether to clean the mesh (merging duplicated
#'   vertices and duplicated faces, removing isolated vertices).
#' @param repairMesh Boolean. Whether to try to remove self-intersections from output mesh.
#' @param method \code{character}. Cost and placement strategy. One of
#'   \code{"LT-R"} (Lindstrom-Turk using undirected edge ratio as stop predicate),
#'   \code{"LT-C"} (Lindstrom-Turk using undirected edge count as stop predicate),
#'   \code{"LT-BNCF"} (Lindstrom-Turk, applying bounded normal change filter), and
#'   \code{"GH"} (Garland-Heckbert).
#' @param stopUeRatio Positive number. Stop predicate. For \code{method="LT-R"}, \code{"GH"}.
#'   Ratio of undirected edges relative to initial count. Lower -> coarser.
#' @param stopUeCount \code{integer}. Stop predicate. For \code{method="LT-C"}, \code{"LT-BNCF"}.
#'   Count of undirected edges left. If missing, set to (half of all edges) - 1.
#'   Lower -> coarser.
#' @param policy \code{character}. For \code{method="GH"}. Simplification strategy.
#'   One of \code{"CP"} (classic plane), \code{"CT"} (classic tri),
#'   \code{"PP"} (prob plane), \code{"PT"} (prob tri), and
#'   \code{"PL"} (plane and line).
#' @param normals Boolean. Whether to return vertex normals.
#' @param verbose Boolean. Whether to print out messages about mesh processing.
#' @return A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/Surface_mesh_simplification/} for details.
#' @seealso \code{\link[SurfaceMesh]{remeshIsotropic}}
#' @author Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh       <- makeMesh(dataHopfTorus, triangulate=TRUE)
#' mesh_rgl   <- toRGL(mesh)
#' mesh_s     <- remeshSimplify(mesh, method="LT-R", stopUeRatio=0.1)
#' mesh_s_rgl <- toRGL(mesh_s)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_s_rgl)
#'
#' @export
remeshSimplify <- function(x,
                           repairSoup = TRUE,
                           repairMesh = TRUE,
                           method = c("LT-R", "LT-C", "LT-BNCF", "GH"),
                           stopUeRatio = 0.1,
                           stopUeCount,
                           policy = c("CP", "CT", "PP", "PT", "PL"),
                           normals = FALSE,
                           verbose = FALSE) {
    if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			     " (i.e., the output of the `makeMesh()` function).")
    }
    stopifnot(isBoolean(repairSoup))
    stopifnot(isBoolean(repairMesh))
    method <- match.arg(method)
    stopifnot(isPositiveNumber(stopUeRatio), stopUeRatio < 1.0)
    stopifnot(isBoolean(normals))
    stopifnot(isBoolean(verbose))
    meshCPP    <- fromR(x)
    meshSimple <- if(method %in% c("LT-R", "LT-C", "LT-BNCF")) {
      nEdges <- nrow(x[["edges"]])
      if(missing(stopUeCount)) {
        stopUeCount <- ceiling((nEdges / 2) - 1)
      } else {
        stopifnot(isStrictPositiveInteger(stopUeCount),
                  stopUeCount < nEdges)
      }

      simplifyLT_cpp(meshCPP,
                     repairSoup,
                     repairMesh,
                     as.character(method),
                     stopUeRatio,
                     as.integer(stopUeCount),
                     normals,
                     verbose)
    } else if(method == "GH") {
      policy  <- match.arg(policy)
      simplifyGH_cpp(meshCPP,
                     repairSoup,
                     repairMesh,
                     stopUeRatio,
                     as.character(policy),
                     normals,
                     verbose)
    }
    fromCPP(meshSimple)
}
