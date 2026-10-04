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

#' @title Hausdorff distance between two meshes (full or quantile)
#' @description Hausdorff distance between two 3D surface meshes. The
#'   approximate distance, the distance estimate with a given error bound,
#'   or the quantile Hausdorff distance, e.g., HD95.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param mesh2 A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param symmetric Boolean. Report the symmetric Hausdorff distance?
#' @param errBound A positive number. If provided, upper bound on the
#'   error of the Hausdorff distance estimate. If missing, and \code{p}
#'   is also missing, the approximate distance is returned.
#' @param p A number in \eqn{[0, 1]}. The quantile probability, defaults
#'   to \code{0.95}. If provided and \code{errBound} is missing,
#'   the quantile Hausdorff distance is returned.
#' @param ... Options for random sampling passed to
#'   \code{\link[SurfaceMesh]{getSurfaceDist}}. See there fore details.
#'   Relevant for the approximate and quantile distance.
#' @returns A number. For the apprixmate and quantile distance,
#'   the algorithm uses random sampling, and thus the result can vary.
#' @details The quantile Hausdorff distance approximates the Hausdorff distance
#'   using a quantile (often the 95th percentile, i.e. "HD95") of the
#'   distances from a random sample of points on one mesh to the other mesh,
#'   instead of their maximum. This is less sensitive to small, isolated
#'   outlying regions than the full Hausdorff distance.
#'
#'   See \url{https://metrics-reloaded.dkfz.de/metric-library/hd} and
#'   \url{https://metrics-reloaded.dkfz.de/metric-library/xhd} for more information
#'   on the (quantile) Hausdorff distance.
#'   See \url{https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#PMPDistance}
#'   for implementation details.
#' @seealso See \code{\link[SurfaceMesh]{getSurfaceDist}} for other surface distance metrics, and
#'   \code{\link[SurfaceMesh]{getJSCDSC}} for volume-overlap based mesh similarity metrics (JSC, DSC).
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' ## approximate symmetric Hausdorff distance
#' getHausdorff(dataHeart1, dataHeart2, symmetric=TRUE)
#'
#' ## estimate with error bound
#' getHausdorff(dataHeart1, dataHeart2, symmetric=TRUE, errBound=0.01)
#'
#' ## quantile Hausdorff distance
#' getHausdorff(dataHeart1, dataHeart2, symmetric=TRUE, p=0.95)
#' @export
getHausdorff <- function(mesh1, mesh2, symmetric = TRUE, errBound, p, ...) {
  stopifnot(inherits(mesh1, "CGALmesh"))
  stopifnot(inherits(mesh2, "CGALmesh"))
  stopifnot(isBoolean(symmetric))
  meshCPP1 <- fromR(mesh1)
  meshCPP2 <- fromR(mesh2)
  if(!missing(errBound)) {
    stopifnot(isPositiveNumber(errBound))
    getHausdorffEst_cpp(meshCPP1, meshCPP2, symmetric, errBound)
  } else if(!missing(p)) {
    metroL <- getSurfaceDist(mesh1,
                             mesh2,
                             symmetric=symmetric,
                             p=p,
                             ...)
    metroL[["HDq"]]
  } else {
    sampleOptL <- checkSampleOpts(list(...))
    getHausdorffApprox_cpp(meshCPP1, meshCPP2, symmetric, sampleOptL)
  }
}

#' @title Several distance metrics between two meshes
#' @description Quantile Hausdorff distance (HD), average symmetric surface distance (ASSD),
#'   and root mean squared error (RMSE) for the surface distance between two 3D surface meshes.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param mesh2 A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param returnDists Boolean. Return the sampled distances from mesh 1 to 2, and from 2 to 1?
#' @param symmetric Boolean. Report the symmetric quantile Hausdorff distance?
#' @param p A number in \eqn{[0, 1]}. The quantile probability, defaults
#'   to \code{0.95}.
#' @param method \code{character}. \code{"random"} for random uniform sampling,
#'   \code{"grid"} for grid sampling, \code{"mc"} for Monte Carlo sampling.
#'   See details.
#' @param sampleVerts Boolean. Do sample vertices?
#' @param sampleEdges Boolean. Do sample edges?
#' @param sampleFaces Boolean. Do sample faces?
#' @param gridSpacing \code{numeric}. The grid spacing for \code{method="grid"}.
#' @param ptsOnEdges \code{integer}. For the random sampling method as the
#'   number of points to pick exclusively on edges.
#'   If missing, the number of edges is used.
#' @param ptsOnFaces \code{integer}. For the random sampling method as the
#'   number of points to pick on the surface.
#'   If missing, the number of vertices is used.
#' @param ptsPerDist \code{numeric}. Points per distance unit. Used for
#'   \code{method="random"} and \code{"mc"} to respectively
#'   determine the total number of points on edges and the number of points per edge.
#' @param ptsPerEdge \code{integer}. Points per edge. Used for \code{method="mc"}
#'   as the number of points per edge to pick.
#' @param ptsPerArea \code{numeric}. Points per area unit. Used for
#'   \code{method="random"} and \code{"mc"} to respectively determine
#'   the total number of points inside faces and the number of points per face.
#' @param ptsPerFace \code{integer}. Points per face. Used for \code{method="mc"}
#'   as the number of points per face to pick.
#' @returns A list with components \code{"HDq"} (quantile Hausdorff distance),
#'   \code{"ASSD"} (average symmetric surface distance),
#'   \code{"RMSE"} (root mean squared error).
#'   For \code{returnDists=TRUE} also \code{"dists12"} and \code{"dists21"}.
#'   The algorithm uses random sampling, so the result can vary.
#' @details See \code{\link[SurfaceMesh]{getHausdorff}} for details on the quantile
#'   Hausdorff distance, \url{https://metrics-reloaded.dkfz.de/metric-library/assd} for ASSD.
#'   Inspired by \url{http://vcglib.net/metro.html}.
#'   For details on sampling options, see
#'   \url{https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__distance__grp.html}.
#' @seealso See \code{\link[SurfaceMesh]{getHausdorff}} for the (approximate) Hausdorff distance, and
#'   \code{\link[SurfaceMesh]{getJSCDSC}} for volume-overlap based mesh similarity metrics (JSC, DSC).
#'
#' @examples
#' library(SurfaceMesh)
#' getSurfaceDist(dataHeart1, dataHeart2, p=0.95, ptsOnFaces=1000L)
#'
#' @export
getSurfaceDist <- function(mesh1,
                           mesh2,
                           returnDists = FALSE,
                           symmetric = TRUE,
                           p = 0.95,
                           method = c("random", "grid", "mc"),
                           sampleVerts = TRUE,
                           sampleEdges = TRUE,
                           sampleFaces = TRUE,
                           gridSpacing = NULL,
                           ptsOnEdges  = NULL,
                           ptsOnFaces  = NULL,
                           ptsPerDist  = NULL,
                           ptsPerEdge  = NULL,
                           ptsPerArea  = NULL,
                           ptsPerFace  = NULL) {
  stopifnot(inherits(mesh1, "CGALmesh"))
  stopifnot(inherits(mesh2, "CGALmesh"))
  stopifnot(isBoolean(returnDists))
  stopifnot(isBoolean(symmetric))
  stopifnot(is.numeric(p), length(p) == 1L, p > 0, p < 1)
  sampleOptL <- checkSampleOpts(list(method     =method,
                                     sampleVerts=sampleVerts,
                                     sampleEdges=sampleEdges,
                                     sampleFaces=sampleFaces,
                                     gridSpacing=gridSpacing,
                                     ptsOnEdges =ptsOnEdges,
                                     ptsOnFaces =ptsOnFaces,
                                     ptsPerDist =ptsPerDist,
                                     ptsPerEdge =ptsPerEdge,
                                     ptsPerArea =ptsPerArea,
                                     ptsPerFace =ptsPerFace))
  meshCPP1 <- fromR(mesh1)
  meshCPP2 <- fromR(mesh2)
  getSurfaceDist_cpp(meshCPP1,
                     meshCPP2,
                     returnDists,
                     symmetric,
                     p,
                     sampleOptL)
}
