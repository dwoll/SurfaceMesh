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

#' @title Approximated centroidal Voronoi diagram remeshing
#' @description Approximated centroidal Voronoi diagram remeshing (ACVD) uses
#'   clustering on polygonal meshes as to approximate a centroidal Voronoi
#'   diagram construction.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#'   The mesh must be a single connected component.
#'   The mesh must be triangle or be able to be made triangle.
#' @param method \code{character}. \code{"uniform_iso"} for uniform isotropic ACVD,
#'   \code{"uniform_iso_qem"} for uniform isotropic ACVD with post-processing QEM optimization,
#'   \code{"uniform_qem"} for uniform ACVD QEM energy minimization,
#'   \code{"adapt_iso"} for adaptive isotropic ACVD.
#' @param nVerts Positive \code{integer}. Lower bound on the target number
#'   of vertices in the output mesh
#' @param vertCountRatio Positive number. Ratio used when input mesh does not
#'   have enough vertices compared to \code{nVerts}. The number of vertices
#'   of the input mesh should at least be \code{vertCountRatio} times
#'   \code{nVerts} If not, the mesh will first be subdivided appropriately.
#' @param gradationFactor Positive number. For \code{method="adapt_iso"}.
#'   A factor used to gradate the weights of the vertices based on their
#'   curvature values. The larger the value is, the more the curvature will
#'   impact the distribution of output vertices on the surface.
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/} for details.
#' @seealso See \code{\link[SurfaceMesh]{remeshDelaunay}},
#'   \code{\link[SurfaceMesh]{remeshIsotropic}},
#'   \code{\link[SurfaceMesh]{remeshSimplify}},
#'   \code{\link[SurfaceMesh]{remeshSmoothShape}},
#'   \code{\link[SurfaceMesh]{remeshSmoothAngle}},
#'   \code{\link[SurfaceMesh]{remeshSmoothTangentRelax}}
#'   for other remeshing operations.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh       <- makeMeshValid(dataOloid)
#' mesh_rgl   <- toRGL(mesh)
#' mesh_d     <- remeshACVD(mesh, method="uniform_iso_qem", nVerts=100L)
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
remeshACVD <- function(x,
                       method = c("uniform_iso", "uniform_iso_qem", "uniform_qem", "adapt_iso"),
                       nVerts = 500L,
                       vertCountRatio = 0.1,
                       gradationFactor = 2.0,
                       normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
           " (i.e., the output of the `makeMesh()` function).")
  }

  method_choices <- c("uniform_iso", "uniform_iso_qem", "uniform_qem", "adapt_iso")
  method         <- match.arg(tolower(method), choices=method_choices)
  stopifnot(isPositiveInteger(nVerts))
  stopifnot(isPositiveNumber(vertCountRatio))
  if(method == "adapt_iso") {
    stopifnot(isPositiveNumber(gradationFactor))
  } else {
    gradationFactor = -1.0
  }
  stopifnot(isBoolean(normals))
  storage.mode(vertCountRatio)  <- "double"
  storage.mode(gradationFactor) <- "double"
  meshCPP <- fromR(x)
  meshOut <- remeshACVD_cpp(
    meshCPP, method, as.integer(nVerts), vertCountRatio, gradationFactor, normals)
  fromCPP(meshOut)
}
